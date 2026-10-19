"""Pilot: four representations of the same cells, one clustering method, Weber-style scoring.

Arms
  A1  arcsinh(x / cofactor), no scaling (CATALYST and diffcyt default)
  A2  A1, then z-score per marker
  A3  A1, then 1st to 99th percentile scaling to [0, 1]
  B   arcsinh((x / brightness) / cofactor): divide by the conjugate's relative brightness on the
      linear scale, before the transform. Only where the fluorochromes are known.
  C   lineage markers as 0/1 (two-component Gaussian mixture on A1), functional markers as A3
  C2  lineage markers as 0/1, functional markers as A2 (z-score), because percentile clipping
      at the 99th percentile flattens a population that is 0.03% of cells
  D   every marker as 0/1
  C_gate, D_gate  the same two binary arms with the supervised gate thresholds from
      results/gate_thresholds.json (code/gate_thresholds.py). Functional markers in C_gate are
      z-scored, as in C2. These arms run only where that file exists for the dataset, and they
      are the upper bound the plan names, not the recommended method.
  D_gate_clean  only the markers whose gate cut separates its positive and negative populations
      with balanced accuracy of at least 0.9 (15 of 32 on Levine_32dim). The other 17 bits are
      noise under any threshold, so this is the upper bound for a binary panel.

Usage: python3 pilot.py <project_dir> [--datasets NAME,NAME] [--arms ARM,ARM] [--tag TAG]
  --datasets  run only these jobs (Nilsson_rare, Mosmann_rare, Levine_32dim, Levine_32dim_k40,
              Levine_13dim, Levine_13dim_k40, Samusik_01, Samusik_01_k40, Samusik_all, Samusik_all_k40). The default is the
              four pilot jobs; the phase 2 jobs run with --datasets and --tag phase2.
  --arms      run only these arms
  --tag       write results/<TAG>_metrics.csv, <TAG>_perpop.csv, <TAG>_meta.json and
              labels_<TAG>_<dataset>.npz instead of the pilot files
Without options it writes results/pilot_metrics.csv, results/pilot_perpop.csv, results/pilot_meta.json.
"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from flowsom.models import FlowSOMEstimator

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(20261018)
SUB_N = 50_000          # subsample for binarization fits, stability and dimensionality
SEEDS = list(range(1, 11))

# Relative brightness, ordinal 1 (dim) to 5 (bright). Source column says where each came from.
# fluor-ai = a curated table of vendor brightness ratings (BioLegend, BD, Cytek).
BRIGHTNESS = {
    "FITC": (3, "fluor-ai"), "PerCP-Cy5.5": (3, "fluor-ai"), "PE": (5, "fluor-ai"),
    "PE-Cy5": (5, "fluor-ai"), "PE-Cy7": (5, "fluor-ai"), "APC": (5, "fluor-ai"),
    "Alexa Fluor 700": (3, "fluor-ai"), "APC-Cy7": (3, "fluor-ai"), "Pacific Blue": (3, "fluor-ai"),
    "Qdot 605": (4, "assumed, not in fluor-ai"), "Qdot 655": (4, "assumed, not in fluor-ai"),
    "Qdot 705": (3, "assumed, not in fluor-ai"), "Qdot 800": (3, "assumed, not in fluor-ai"),
    "Alexa Fluor 488": (3, "fluor-ai"), "PE-Texas Red": (4, "fluor-ai"), "AmCyan": (2, "fluor-ai, as V500"),
    "Qdot 585": (3, "assumed, not in fluor-ai"),
    # DURAClone IM Basic plus drop-ins, the local LSRFortessa set (code/stain_index.py)
    "BUV661": (4, "fluor-ai"), "BUV737": (4, "fluor-ai"), "BV786": (4, "fluor-ai"),
    "APC-Alexa Fluor 750": (3, "assumed, not in fluor-ai; rated with the other APC tandems"),
    "Krome Orange": (3, "assumed, not in fluor-ai"),
    "Viability dye, violet": (3, "assumed; the FCS keywords name no product, fluor-ai rates every violet viability dye 3"),
}
# Two provisional mappings from the ordinal rank to a multiplicative divisor.
BRIGHT_MAPS = {"B_geom": lambda b: 2.0 ** (b - 3), "B_lin": lambda b: b / 3.0}


# ----------------------------------------------------------------------------- loaders
def load_levine32(root):
    import flowio
    f = flowio.FlowData(str(root / "data/FR-FCM-ZZPH/Levine_32dim_notransform.fcs"))
    ch = [f.channels[i]["pnn"] for i in sorted(f.channels)]
    d = np.reshape(f.events, (-1, f.channel_count)).astype(np.float32)
    skip = {"Time", "Cell_length", "DNA1", "DNA2", "Viability", "file_number", "event_number", "label", "individual"}
    cols = [i for i, c in enumerate(ch) if c not in skip]
    markers = [ch[i] for i in cols]
    lab = d[:, ch.index("label")]
    lab = np.where(np.isnan(lab), 0, lab).astype(int)
    names = pd.read_csv(root / "data/FR-FCM-ZZPH/population_names_Levine_32dim.txt", sep="\t")
    pop_names = dict(zip(names["label"].astype(int), names["population"]))
    continuous = {"CXCR4", "CD44", "CD47", "CD49d", "CD321", "CD38", "HLA-DR"}
    return dict(name="Levine_32dim", kind="mass", X=d[:, cols], markers=markers, labels=lab,
                pop_names=pop_names, cofactor=5.0, rare=None,
                lineage=[m for m in markers if m not in continuous], fluor=None, k=14)


def load_npz(root, name, cofactor, keep, lineage, fluor, rare_label=2):
    z = np.load(root / f"data/HDCytoData/{name}.npz", allow_pickle=True)
    markers_all = list(z["marker"])
    cols = [markers_all.index(m) for m in keep]
    lab = (z["pop"] == rare_label).astype(int)   # 1 = the rare population, 0 = everything else
    return dict(name=name, kind="flow", X=z["expr"][:, cols].astype(np.float32), markers=keep,
                labels=lab, pop_names={1: "rare"}, cofactor=cofactor, rare=1,
                lineage=lineage, fluor=fluor, k=40)


def load_nilsson(root):
    keep = ["CD38", "CD34", "CD123", "CD19", "CD10", "CD110", "CD45", "CD3", "CD45RA", "CD11b", "CD4", "CD49fpur", "CD90bio"]
    fluor = {"CD38": "FITC", "CD34": "PerCP-Cy5.5", "CD123": "PE", "CD19": "PE-Cy5", "CD10": "PE-Cy7",
             "CD110": "APC", "CD45": "Alexa Fluor 700", "CD3": "APC-Cy7", "CD45RA": "Pacific Blue",
             "CD11b": "Qdot 605", "CD4": "Qdot 655", "CD49fpur": "Qdot 705", "CD90bio": "Qdot 800"}
    lineage = [m for m in keep if m != "CD38"]
    return load_npz(root, "Nilsson_rare", 150.0, keep, lineage, fluor)


def load_mosmann(root):
    """Panel from the JMW034 MIFlowCyt checklist, Supporting Information Table 1 of Mosmann et al.
    2014 (PMC4238823, file cyto0085-0422-sd1.docx, sections 2.4 and 3.3). LSR-II, PE and its
    tandems on the 532 nm laser. The Live/Dead Pacific Orange channel is not in the HDCytoData set."""
    typ = ["CXCR5", "CD3", "CD45RA", "CD8a", "CD14", "CD69", "CD4"]
    state = ["IL17A", "TNFa", "GZB-SA", "IL5", "IL2", "CCL4", "IFNg"]
    fluor = {"CXCR5": "Alexa Fluor 488", "IL17A": "PerCP-Cy5.5", "TNFa": "Pacific Blue", "GZB-SA": "Qdot 585",
             "CD3": "Qdot 605", "CD45RA": "Qdot 655", "CD8a": "Qdot 705", "CD14": "Qdot 800", "IL5": "APC",
             "IL2": "Alexa Fluor 700", "CD69": "APC-Cy7", "CCL4": "PE", "CD4": "PE-Texas Red", "IFNg": "PE-Cy7"}
    return load_npz(root, "Mosmann_rare", 150.0, typ + state, typ, fluor)


def load_levine13(root):
    """Levine_13dim: 13 surface markers, 24 gated populations, 1 donor, mass cytometry.
    The label set defines CD11b at three levels (CD11b-, mid, hi monocytes) and CD38 at two
    (CD38lo, CD38mid mature B cells), which are the hypothesis 3 cases the plan names.
    CD38 stays continuous as in Levine_32dim. CD11b stays in the lineage set on purpose, so the
    binary arms show what a single cut does to a three-level marker."""
    import flowio
    f = flowio.FlowData(str(root / "data/FR-FCM-ZZPH/Levine_13dim_notransform.fcs"))
    ch = [f.channels[i]["pnn"] for i in sorted(f.channels)]
    d = np.reshape(f.events, (-1, f.channel_count)).astype(np.float32)
    markers = [c for c in ch if c != "label"]
    cols = [ch.index(m) for m in markers]
    lab = d[:, ch.index("label")]
    lab = np.where(np.isnan(lab), 0, lab).astype(int)
    names = pd.read_csv(root / "data/FR-FCM-ZZPH/population_names_Levine_13dim.txt", sep="\t")
    pop_names = dict(zip(names["label"].astype(int), names["population"]))
    continuous = {"CD38"}
    return dict(name="Levine_13dim", kind="mass", X=d[:, cols], markers=markers, labels=lab,
                pop_names=pop_names, cofactor=5.0, rare=None,
                lineage=[m for m in markers if m not in continuous], fluor=None, k=24)


def load_samusik01(root, name="Samusik_01"):
    """Samusik_01: 39 surface markers, 24 gated populations, mouse bone marrow, mass cytometry.
    Population 25 is the unassigned class and scores as unlabeled. The continuous set holds the
    activation, maturation and stem markers that the label set uses at more than one level
    (CD44, CD25, CD27, Ly6C, Sca1, CD150, CD43, CD103, CCR7, CD16_32, Foxp3).
    Samusik_all is the same panel on all 10 mice (841,644 cells); rowvec0 is the mouse id and the
    dataset dict keeps it as "sample" for the between-mouse variance analysis."""
    z = np.load(root / f"data/HDCytoData/{name}.npz", allow_pickle=True)
    markers_all = [str(m) for m in z["marker"]]
    skip = {"Time", "Cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "DNA1", "DNA2", "Cisplatin", "beadDist"}
    markers = [m for m in markers_all if m not in skip]
    cols = [markers_all.index(m) for m in markers]
    lab = z["pop"].astype(int)
    lab = np.where(lab == 25, 0, lab)
    pop_names = {int(k): v for k, v in json.load(open(root / "data/HDCytoData/population_names_Samusik.json")).items() if int(k) != 25}
    continuous = {"CD44", "CD25", "CD27", "Ly6C", "Sca1", "CD150", "CD43", "CD103", "CCR7", "CD16_32", "Foxp3"}
    return dict(name=name, kind="mass", X=z["expr"][:, cols].astype(np.float32), markers=markers,
                labels=lab, pop_names=pop_names, cofactor=5.0, rare=None,
                lineage=[m for m in markers if m not in continuous], fluor=None, k=24,
                sample=z["rowvec0"].astype(int))


# HIPC Lyoplate T cell panel: fixed BD conjugates, the same at every center (Finak 2016 Table 1 and
# the single-stain control names in data/HIPC/COMP). Never read these from $PnS, which is stale in
# the deposited files; see data/README.md.
HIPC_TCELL_FLUOR = {"Live": "FITC", "CCR7": "PE", "CD4": "PerCP-Cy5.5", "CD45RA": "PE-Cy7", "CD3": "Pacific Blue",
                    "HLADR": "AmCyan", "CD38": "APC", "CD8": "APC-Cy7"}
# The transform the depositor applied is FlowJo 9's, as flowWorkspace reproduced it from the
# central-analysis workspaces: a logicle with T 262144, M 4.5 decades, A 0, scaled to 4096,
# and a width W that FlowJo chose per parameter. The widths below were fitted on the Stanford
# raw files against the deposited ones (median error 1 to 6 units of 4096; code/hipc_transform_check.py)
# and are applied to every center, because only Stanford deposited raw files. The README's
# iplogicle (T 16777215, M 7.22, W 1.16) does not reproduce the deposited values.
HIPC_LOGICLE = dict(t=262144.0, m=4.5, a=0.0, scale=4096.0,
                    w={"Live": 0.40, "CCR7": 0.90, "CD4": 0.80, "CD45RA": 0.80, "CD3": 0.15, "HLADR": 0.75, "CD38": 0.60, "CD8": 1.10})


def load_hipc_tcell(root, per_file=40_000, centers=("Baylor", "CIMR", "Miami", "NHLBI", "Stanford", "UCLA", "Yale")):
    """HIPC Lyoplate T cell panel, the seven centers with central manual gates. The deposited
    files are compensated and iplogicle-transformed by the HIPC pipeline, so the dataset carries
    the transformed matrix as X and kind "pretransformed": build_arms uses X as A1 directly and
    inverts the logicle for arm B. Scatter, Time and Live are not features; Live stays in the
    dict for the between-site script. Each file contributes up to per_file events, drawn with
    a fixed seed, so the pooled clustering weighs every center and lot the same.

    Labels come from the central manual gates applied in Python (code/flowjo9_gates.py, files
    data/HIPC/labels/<Center>/<file>.npz, one entry per event of the deposited file): 1..8 are
    the CD4 and CD8 naive, central memory, effector memory and effector subsets, 0 is everything
    the analyst did not put in one of those (debris, doublets where a singlets gate could be
    applied, dead, non-T, double positive and double negative T). `activated` is the CD38+
    HLA-DR+ flag within CD4 or CD8, kept beside the labels, not in them."""
    import flowio
    hdir = root / "data/HIPC"
    markers = ["CCR7", "CD4", "CD45RA", "CD3", "HLADR", "CD38", "CD8"]
    Xs, live, site, lot, rep, fname, labs, act = [], [], [], [], [], [], [], []
    rng = np.random.default_rng(20261019)
    for c in centers:
        pd_ = pd.read_csv(hdir / f"pData/{c}.csv")
        pd_.columns = ["file", "institution", "panel", "replicate", "sample"]
        for _, r in pd_[pd_.panel == "T cell"].iterrows():
            f = hdir / f"FCS/{c}/{r.file}"
            if not f.exists(): continue
            fd = flowio.FlowData(str(f)); t = fd.text; n = fd.channel_count
            pnn = [t.get(f"p{i}n", t.get(f"$P{i}N")) for i in range(1, n + 1)]
            d = np.reshape(fd.events, (-1, n)).astype(np.float32)
            cols = [pnn.index(m) for m in markers]
            # the depositor's transform wrote absurd negatives (down to -6e7) on 0.6 percent of
            # events; the logicle floor on this scale is about -1100, so anything below -4096 is broken
            idx = np.flatnonzero((d[:, cols] > -4096).all(1))
            if len(idx) > per_file: idx = idx[rng.choice(len(idx), per_file, replace=False)]   # same draw as before the labels existed
            d = d[idx]
            lf = hdir / "labels" / c / f"{r.file}.npz"
            if not lf.exists(): raise FileNotFoundError(f"{lf}: run code/flowjo9_gates.py --write first")
            with np.load(lf) as z:
                labs.append(z["label"][idx].astype(int)); act.append(z["activated"][idx]); pop_names = [str(n) for n in z["pop_names"]]
            Xs.append(d[:, cols]); live.append(d[:, pnn.index("Live")])
            site += [c] * len(d); lot += [str(r["sample"])] * len(d); rep += [int(r.replicate)] * len(d); fname += [r.file] * len(d)
    X = np.vstack(Xs)
    return dict(name="HIPC_tcell", kind="pretransformed", X=X, markers=markers, labels=np.concatenate(labs),
                pop_names={i: n for i, n in enumerate(pop_names) if i > 0}, cofactor=None, rare=None,
                lineage=["CCR7", "CD4", "CD45RA", "CD3", "CD8"],
                fluor={m: HIPC_TCELL_FLUOR[m] for m in markers}, k=12, logicle=HIPC_LOGICLE, activated=np.concatenate(act),
                live=np.concatenate(live), site=np.array(site), lot=np.array(lot), replicate=np.array(rep), file=np.array(fname),
                site_divisors=hipc_site_divisors(root, markers, centers))


def hipc_site_divisors(root, markers, centers):
    """Measured brightness divisors from results/hipc_stain_index.csv (code/hipc_stain_index.py),
    the single-stain controls of each center. Two divisor sources:
      B_site: per center, the stain index of each marker's dye on that center's own instrument,
              normalized to geometric mean 1 over the markers within the center, so the cofactor
              keeps its meaning and only the relative brightness of the channels changes. Where a
              center has no control for a dye the CSV carries the across-center mean in
              stain_index_used, and that cell takes the pooled value.
      B_pool: the across-center mean stain index per dye, one divisor set for every center,
              normalized the same way. B_pool against B_geom is the measured-versus-rank effect;
              B_site against B_pool is the per-instrument effect.
    None when the CSV is absent, and build_arms then makes neither arm."""
    f = Path(root) / "results/hipc_stain_index.csv"
    if not f.exists(): return None
    si = pd.read_csv(f)
    def geo1(v):
        v = np.asarray(v, float); return v / np.exp(np.mean(np.log(v)))
    dyes = [HIPC_TCELL_FLUOR[m] for m in markers]
    pooled = si.drop_duplicates("dye").set_index("dye").pooled_stain_index
    out = {"B_site": {}, "B_pool": dict(zip(markers, geo1([pooled[d] for d in dyes]).round(4).tolist()))}
    for c in centers:
        s = si[si.center == c].set_index("dye").stain_index_used
        out["B_site"][c] = dict(zip(markers, geo1([s[d] for d in dyes]).round(4).tolist()))
    return out


# Mass cytometry has no dye, but the detector's sensitivity varies with isotope mass. This table
# holds the relative transmission per mass channel of the first-generation CyTOF that acquired the
# Levine sets (Bendall 2011), as a divisor for arm B_meas; values and source are set in
# results/mass_vs_flow.md. The Levine_13dim metals come from the HDCytoData channel names.
# Relative sensitivity of the first-generation CyTOF by mass, percent of each instrument's peak,
# mean of the three CyTOF 1 instruments in Tricot et al. 2015 Cytometry A 87:357 Figure 2B (values
# read off the figure, about 3 units of uncertainty). 145 is interpolated between 144 and 146.
# Indium 115 is below the measured range (139 to 176); the curve falls toward low mass, so it
# takes 0.8 of the 139 value, an assumption. The summed cadmium channel (CD3, 110 to 114) has
# no published value and takes the median of the others, which makes it neutral in the divisor.
CYTOF_SENSITIVITY = {"139": 26.0, "142": 29.0, "144": 29.7, "145": 32.0, "146": 34.3, "147": 37.5, "148": 40.7,
                     "158": 84.7, "160": 92.7, "167": 75.7, "170": 74.7, "115": 20.8}
CYTOF_SENSITIVITY["110_114"] = float(sorted(CYTOF_SENSITIVITY.values())[len(CYTOF_SENSITIVITY) // 2])
LEVINE13_METAL = {"CD45": "115", "CD45RA": "139", "CD19": "142", "CD11b": "144", "CD4": "145", "CD8": "146",
                  "CD34": "148", "CD20": "147", "CD33": "158", "CD123": "160", "CD38": "167", "CD90": "170", "CD3": "110_114"}
SHARED_MARKERS = ["CD45", "CD45RA", "CD19", "CD11b", "CD4", "CD34", "CD123", "CD38", "CD90", "CD3"]


def load_levine13_shared(root):
    """Levine_13dim on the 10 markers it shares with Nilsson_rare, scored on the HSC (label 7,
    261 cells) as a rare population so the two sets use the same metric. CD38 stays continuous.
    Arm B_meas divides by the CyTOF sensitivity of each marker's mass channel."""
    ds = load_levine13(root)
    cols = [ds["markers"].index(m) for m in SHARED_MARKERS]
    lab = (ds["labels"] == 7).astype(int)
    div = {m: CYTOF_SENSITIVITY.get(LEVINE13_METAL[m]) for m in SHARED_MARKERS}
    if any(v is None for v in div.values()): div = None
    return dict(ds, name="Levine_13dim_shared", X=ds["X"][:, cols], markers=list(SHARED_MARKERS), labels=lab,
                pop_names={1: "HSC"}, rare=1, lineage=[m for m in SHARED_MARKERS if m != "CD38"], fluor=None, k=40,
                divisor_values=div, metal={m: LEVINE13_METAL[m] for m in SHARED_MARKERS})


def load_nilsson_shared(root):
    """Nilsson_rare on the same 10 markers (CD90bio is CD90 on Qdot 800), HSC as the rare population."""
    ds = load_nilsson(root)
    names = {"CD90bio": "CD90"}
    have = [names.get(m, m) for m in ds["markers"]]
    cols = [have.index(m) for m in SHARED_MARKERS]
    fluor = {SHARED_MARKERS[i]: ds["fluor"][ds["markers"][c]] for i, c in enumerate(cols)}
    return dict(ds, name="Nilsson_shared", X=ds["X"][:, cols], markers=list(SHARED_MARKERS),
                lineage=[m for m in SHARED_MARKERS if m != "CD38"], fluor=fluor, k=40)


# ----------------------------------------------------------------------------- representations
def gmm_fit(v):
    """Two-component Gaussian mixture on a subsample. Returns the equal-posterior cut, the low and
    high means, the low and high variances, the fitted model and the index of the high component."""
    s = v if len(v) <= SUB_N else RNG.choice(v, SUB_N, replace=False)
    g = GaussianMixture(2, random_state=0).fit(s.reshape(-1, 1))
    lo, hi = np.argsort(g.means_.ravel())
    grid = np.linspace(g.means_.ravel()[lo], g.means_.ravel()[hi], 2000).reshape(-1, 1)
    post = g.predict_proba(grid)[:, hi]
    i = np.argmax(post >= 0.5)
    mu, var = g.means_.ravel(), g.covariances_.ravel()
    return float(grid[i, 0]), float(mu[lo]), float(mu[hi]), float(var[lo]), float(var[hi]), g, int(hi)


def gmm_threshold(v):
    """Equal-posterior cut between the two components of a 2-GMM fitted on a subsample."""
    return gmm_fit(v)[0]


def resolution_weights(fits, cap=3.0):
    """Arm W: a weight per marker from the separation of its two mixture modes in pooled-SD units,
    (mu_high - mu_low) / sqrt((var_high + var_low) / 2), a stain index measured on the cells.
    Normalized to geometric mean 1, then capped at `cap` and floored at 1/cap so a unimodal
    continuous marker stays in."""
    w = np.array([(f[2] - f[1]) / np.sqrt((f[3] + f[4]) / 2) for f in fits], dtype=np.float64)
    w = np.maximum(w, 1e-3)
    w = w / np.exp(np.mean(np.log(w)))
    return np.clip(w, 1 / cap, cap).astype(np.float32)


def percentile01(a):
    lo, hi = np.percentile(a, [1, 99], axis=0)
    return np.clip((a - lo) / np.maximum(hi - lo, 1e-9), 0, 1)


def load_gate_thresholds(root, ds, min_balanced_accuracy=0.9):
    """Supervised thresholds from code/gate_thresholds.py, or None when the file is absent.
    Returns the threshold per marker and a mask of the markers whose cut reaches the balanced
    accuracy floor; the mask defines arm D_gate_clean."""
    f = Path(root) / "results/gate_thresholds.json"
    if not f.exists(): return None
    g = json.load(open(f))
    if g["meta"]["dataset"] != ds["name"].replace("_k40", ""): return None
    thr = np.array([g["markers"][m]["threshold"] for m in ds["markers"]], dtype=np.float32)
    clean = np.array([g["markers"][m].get("balanced_accuracy", 0) >= min_balanced_accuracy for m in ds["markers"]])
    return thr, clean


def build_arms(ds, gate=None):
    X, c = ds["X"], ds["cofactor"]
    if ds.get("kind") == "pretransformed":
        A1 = X.copy()     # already compensated and logicle-transformed by the depositor
    else:
        A1 = np.arcsinh(X / c)
    A2 = (A1 - A1.mean(0)) / np.maximum(A1.std(0), 1e-9)
    A3 = percentile01(A1)
    fits = [gmm_fit(A1[:, j]) for j in range(A1.shape[1])]
    thr = np.array([f[0] for f in fits], dtype=np.float32)
    binary = (A1 > thr).astype(np.float32)
    lin = np.array([m in ds["lineage"] for m in ds["markers"]])
    C = np.where(lin, binary, A3)
    C2 = np.where(lin, binary, A2)
    # W: intensity weighted by each marker's mode separation; W_z: the same weights on z-scores;
    # P: the mixture posterior of the high component in place of the bit, functional markers z-scored
    w = resolution_weights(fits)
    post = np.column_stack([f[5].predict_proba(A1[:, [j]].astype(np.float64))[:, f[6]] for j, f in enumerate(fits)]).astype(np.float32)
    # W2: weight by the sharpness of the posterior, the share of cells scored below 0.1 or above 0.9.
    # It rewards a marker whose two modes are clean wherever they sit and ignores how many cells
    # are in the high mode. W2s squares the sharpness before normalizing, for a wider spread.
    sharp = ((post < 0.1) | (post > 0.9)).mean(0).astype(np.float64)
    def norm(v, cap=3.0):
        v = np.maximum(v, 1e-3); v = v / np.exp(np.mean(np.log(v))); return np.clip(v, 1 / cap, cap).astype(np.float32)
    w2, w2s = norm(sharp), norm(sharp ** 2)
    arms = {"A1": A1, "A2": A2, "A3": A3, "C": C, "C2": C2, "D": binary,
            "W": (A1 * w).astype(np.float32), "W_z": (A2 * w).astype(np.float32), "P": np.where(lin, post, A2).astype(np.float32),
            "W2": (A1 * w2).astype(np.float32), "W2s": (A1 * w2s).astype(np.float32)}
    if gate is not None:
        thr_gate, clean = gate
        gbin = (A1 > thr_gate).astype(np.float32)
        arms["C_gate"] = np.where(lin, gbin, A2)
        arms["D_gate"] = gbin
        arms["D_gate_clean"] = gbin[:, clean]
    divisors = {"W_weights": dict(zip(ds["markers"], w.round(3).tolist())),
                "W2_sharpness": dict(zip(ds["markers"], sharp.round(3).tolist())),
                "W2_weights": dict(zip(ds["markers"], w2.round(3).tolist())), "W2s_weights": dict(zip(ds["markers"], w2s.round(3).tolist()))}
    if ds.get("divisor_values"):
        # a measured per-channel sensitivity or stain index, used directly as the divisor
        div = np.array([ds["divisor_values"][m] for m in ds["markers"]], dtype=np.float32)
        div = div / np.exp(np.mean(np.log(div)))
        arms["B_meas"] = np.arcsinh((X / div) / c)
        divisors["B_meas"] = dict(zip(ds["markers"], div.round(3).tolist()))
    if ds["fluor"]:
        for mname, fmap in BRIGHT_MAPS.items():
            div = np.array([fmap(BRIGHTNESS[ds["fluor"][m]][0]) for m in ds["markers"]], dtype=np.float32)
            div = div / np.exp(np.mean(np.log(div)))     # geometric mean 1, so the cofactor keeps its meaning
            arms[mname] = divide_channels(X, ds, div)
            divisors[mname] = dict(zip(ds["markers"], div.round(3).tolist()))
    # measured divisors from the dataset's own single-stain controls (hipc_site_divisors)
    for mname, src in (ds.get("site_divisors") or {}).items():
        if isinstance(next(iter(src.values())), dict):          # one divisor set per center
            table = {s_: np.array([src[s_][m] for m in ds["markers"]], np.float32) for s_ in src}
            rows = np.stack([table[s_] for s_ in table]); idx = {s_: i for i, s_ in enumerate(table)}
            arms[mname] = divide_channels(X, ds, rows[np.array([idx[s_] for s_ in ds["site"]])])
            divisors[mname] = {s_: dict(zip(ds["markers"], v.round(3).tolist())) for s_, v in table.items()}
        else:                                                    # one divisor set for every cell
            div = np.array([src[m] for m in ds["markers"]], np.float32)
            arms[mname] = divide_channels(X, ds, div)
            divisors[mname] = dict(zip(ds["markers"], div.round(3).tolist()))
    return arms, thr, lin, divisors


def divide_channels(X, ds, div):
    """Arm B's division on the linear scale. div is one divisor per marker, or one per cell and
    marker as a (cells, markers) array. A pretransformed dataset goes back to the linear scale
    through the depositor's logicle, per channel because the width differs per channel, and
    forward again through the same logicle; any other dataset is linear already and takes the
    arcsinh with the dataset's cofactor."""
    div = np.asarray(div, dtype=np.float64)
    if ds.get("kind") != "pretransformed":
        return np.arcsinh((X / div) / ds["cofactor"]).astype(np.float32)
    from flowutils import transforms as FT
    L = ds["logicle"]; out = np.empty_like(X, dtype=np.float32)
    for j, m in enumerate(ds["markers"]):
        w = L["w"][m] if isinstance(L["w"], dict) else L["w"]
        d = div[j] if div.ndim == 1 else div[:, [j]]
        linear = FT.logicle_inverse(X[:, [j]].astype(np.float64) / L["scale"], [0], t=L["t"], m=L["m"], w=w, a=L["a"])
        out[:, j] = (FT.logicle(linear / d, [0], t=L["t"], m=L["m"], w=w, a=L["a"]) * L["scale"]).ravel()
    return out


# ----------------------------------------------------------------------------- scoring
def f1_matrix(true, pred, pops, clusters):
    M = np.zeros((len(pops), len(clusters)))
    for i, p in enumerate(pops):
        t = true == p
        for j, q in enumerate(clusters):
            c = pred == q
            tp = np.sum(t & c)
            if tp == 0: continue
            prec, rec = tp / c.sum(), tp / t.sum()
            M[i, j] = 2 * prec * rec / (prec + rec)
    return M


def hungarian_f1(true, pred):
    """Weber and Robinson: one-to-one matching that maximises summed F1, unweighted mean."""
    pops = sorted(set(true)); clusters = sorted(set(pred))
    M = f1_matrix(true, pred, pops, clusters)
    r, c = linear_sum_assignment(-M)
    per = {pops[i]: M[i, j] for i, j in zip(r, c)}
    for p in pops: per.setdefault(p, 0.0)
    return float(np.mean(list(per.values()))), per


def rare_f1(true, pred, rare):
    t = true == rare; best = (0.0, 0.0, 0.0)
    for q in set(pred):
        c = pred == q; tp = np.sum(t & c)
        if tp == 0: continue
        prec, rec = tp / c.sum(), tp / t.sum(); f1 = 2 * prec * rec / (prec + rec)
        if f1 > best[0]: best = (f1, prec, rec)
    return best


def cluster(X, k, seed):
    return FlowSOMEstimator(n_clusters=k, xdim=10, ydim=10, seed=seed).fit_predict(X)


def two_nn_id(X, n=10_000):
    s = X[RNG.choice(len(X), min(n, len(X)), replace=False)]
    s = s + RNG.normal(0, 1e-6, s.shape)   # break exact ties in binary data
    d, _ = NearestNeighbors(n_neighbors=3).fit(s).kneighbors(s)
    mu = d[:, 2] / np.maximum(d[:, 1], 1e-12)
    mu = np.sort(mu)[: int(0.9 * len(mu))]
    return float(len(mu) / np.sum(np.log(mu)))


def dimensionality(Xa, binary_cols):
    s = Xa[RNG.choice(len(Xa), min(SUB_N, len(Xa)), replace=False)]
    pca = PCA(n_components=min(5, s.shape[1])).fit(s)
    out = dict(pca5=float(pca.explained_variance_ratio_.sum()), pca1=float(pca.explained_variance_ratio_[0]),
               two_nn_id=two_nn_id(s))
    if binary_cols is not None and binary_cols.any():
        B = Xa[:, binary_cols]
        out["occupied_phenotypes"] = int(len(np.unique(B, axis=0)))
        out["possible_phenotypes"] = int(2 ** binary_cols.sum())
        counts = np.unique(B, axis=0, return_counts=True)[1]
        out["phenotypes_covering_95pct"] = int(np.searchsorted(np.cumsum(np.sort(counts)[::-1]) / counts.sum(), 0.95) + 1)
    return out


# ----------------------------------------------------------------------------- main
def run(root, datasets=None, only_arms=None, tag="pilot"):
    root = Path(root); res = root / "results"; res.mkdir(exist_ok=True)
    rows, perpop, meta = [], [], {}
    loaders = {"Nilsson_rare": load_nilsson, "Mosmann_rare": load_mosmann, "Levine_32dim": load_levine32,
               "Levine_32dim_k40": lambda r: dict(load_levine32(r), name="Levine_32dim_k40", k=40),
               "Levine_13dim": load_levine13,
               "Levine_13dim_k40": lambda r: dict(load_levine13(r), name="Levine_13dim_k40", k=40),
               "Samusik_01": load_samusik01,
               "Samusik_01_k40": lambda r: dict(load_samusik01(r), name="Samusik_01_k40", k=40),
               "Samusik_all": lambda r: load_samusik01(r, "Samusik_all"),
               "Samusik_all_k40": lambda r: dict(load_samusik01(r, "Samusik_all"), name="Samusik_all_k40", k=40),
               "Levine_13dim_shared": load_levine13_shared, "Nilsson_shared": load_nilsson_shared,
               "HIPC_tcell": load_hipc_tcell,
               "HIPC_tcell_k30": lambda r: dict(load_hipc_tcell(r), name="HIPC_tcell_k30", k=30),
               # one center at a time, same 40,000 events per file and k 12, so the site effect of the
               # pooled job separates from the representation effect (results/hipc.md, F1 section)
               **{f"HIPC_tcell_{c}": (lambda r, c=c: dict(load_hipc_tcell(r, centers=(c,)), name=f"HIPC_tcell_{c}"))
                  for c in ("Baylor", "CIMR", "Miami", "NHLBI", "Stanford", "UCLA", "Yale")}}
    names = datasets or ["Nilsson_rare", "Mosmann_rare", "Levine_32dim", "Levine_32dim_k40"]
    # a job name may carry "@cf=<cofactor>" to run the same cells with another arcsinh cofactor
    def make(n):
        base, _, cf = n.partition("@cf=")
        ds = loaders[base](root)
        if cf: ds = dict(ds, name=f"{base}@cf={cf}", cofactor=float(cf))
        return ds
    jobs = [make(n) for n in names]   # the k40 job is the same cells with a free number of metaclusters
    label_prefix = "labels" if tag == "pilot" else f"labels_{tag}"
    for ds in jobs:
        t0 = time.time()
        gate = load_gate_thresholds(root, ds)
        arms, thr, lin, divisors = build_arms(ds, gate)
        if only_arms: arms = {a: arms[a] for a in only_arms if a in arms}
        labeled = ds["labels"] > 0
        meta[ds["name"]] = dict(cells=int(len(ds["X"])), labeled=int(labeled.sum()), markers=ds["markers"],
                                lineage=ds["lineage"], thresholds=dict(zip(ds["markers"], thr.round(3).tolist())),
                                divisors=divisors, k=ds["k"], cofactor=ds["cofactor"], fluor=ds["fluor"],
                                gate_thresholds=None if gate is None else dict(zip(ds["markers"], gate[0].round(3).tolist())),
                                gate_clean_markers=None if gate is None else [m for m, c in zip(ds["markers"], gate[1]) if c])
        sub = RNG.choice(len(ds["X"]), min(SUB_N, len(ds["X"])), replace=False)
        full_labels = {}
        for arm, Xa in arms.items():
            pred = cluster(Xa, ds["k"], seed=1); full_labels[arm] = pred
            row = dict(dataset=ds["name"], arm=arm, k=ds["k"])
            if ds["rare"] is None and labeled.any():
                f1, per = hungarian_f1(ds["labels"][labeled], pred[labeled])
                row.update(mean_f1=f1, ari=adjusted_rand_score(ds["labels"][labeled], pred[labeled]),
                           nmi=normalized_mutual_info_score(ds["labels"][labeled], pred[labeled]))
                for p, v in per.items():
                    perpop.append(dict(dataset=ds["name"], arm=arm, population=ds["pop_names"][p], f1=v,
                                       n=int(np.sum(ds["labels"] == p))))
            elif ds["rare"] is not None:
                f1, prec, rec = rare_f1(ds["labels"], pred, ds["rare"])
                row.update(rare_f1=f1, rare_precision=prec, rare_recall=rec)
            # stability: same 50k cells, ten seeds, mean pairwise ARI
            labs = [cluster(Xa[sub], ds["k"], seed=s) for s in SEEDS]
            pair = [adjusted_rand_score(labs[i], labs[j]) for i in range(len(labs)) for j in range(i + 1, len(labs))]
            row.update(stability_ari=float(np.mean(pair)), n_clusters_found=int(len(np.unique(pred))))
            bcols = lin if arm in ("C", "C2", "C_gate") else (np.ones(Xa.shape[1], bool) if arm in ("D", "D_gate", "D_gate_clean") else None)
            row.update(dimensionality(Xa, bcols))
            rows.append(row); print(ds["name"], arm, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items() if k not in ("dataset", "arm")}, flush=True)
        # agreement between arms on the full clustering (hypothesis 1)
        for a in full_labels:
            for b in full_labels:
                if a < b:
                    rows.append(dict(dataset=ds["name"], arm=f"{a}~{b}", k=ds["k"],
                                     arm_agreement_ari=adjusted_rand_score(full_labels[a], full_labels[b])))
        extra = {("fcs_file" if k == "file" else k): ds[k] for k in ("sample", "site", "lot", "replicate", "file", "live") if k in ds}
        np.savez_compressed(res / f"{label_prefix}_{ds['name']}.npz", **full_labels, truth=ds["labels"], **extra)
        print(f"{ds['name']} done in {time.time() - t0:.0f} s", flush=True)
    pd.DataFrame(rows).to_csv(res / f"{tag}_metrics.csv", index=False)
    pd.DataFrame(perpop).to_csv(res / f"{tag}_perpop.csv", index=False)
    json.dump(meta, open(res / f"{tag}_meta.json", "w"), indent=1)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--datasets", default=None, help="comma-separated job names")
    ap.add_argument("--arms", default=None, help="comma-separated arm names")
    ap.add_argument("--tag", default="pilot", help="output file prefix")
    a = ap.parse_args()
    run(a.root, a.datasets.split(",") if a.datasets else None, a.arms.split(",") if a.arms else None, a.tag)
