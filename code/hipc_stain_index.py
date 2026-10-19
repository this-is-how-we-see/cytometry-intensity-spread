"""Stain index per center and dye from the HIPC Lyoplate single-stain compensation controls.

Arm B divides each channel by the relative brightness of its fluorochrome. The phase 3 run
(results/hipc.md) used vendor ranks, which are the same at every center, and the between-site
spread went up. This script measures the brightness on each center's own instrument, so that
pilot.py can build arm B_site (one divisor set per center) and arm B_pool (the mean over
centers, one divisor set for all).

Data: data/HIPC/COMP/<Center>/*.fcs, the raw FACSDiva exports of each center's compensation
controls, linear scale, with a SPILL keyword. The antibody controls are BD CompBeads, a tube of
capture beads with about half positive and half negative beads. The Live/Dead FITC controls are
cells. Yale has two sets (Compensation Controls_* and Specimen_001_*); the Specimen_001 set is
beads and the Compensation Controls FITC and PE tubes are cells. NHLBI names some tubes by
detector (B515, B710, V450, V545).

Method, per center and per T cell panel dye, on the uncompensated primary detector:
  1. the file: the control whose name carries the dye, with the panel's own reagent first
     (Live FITC, CCR7 PE, CD4 PerCP-Cy5.5, CD45RA PE-Cy7, CD3 V450, HLA-DR V500, CD38 APC,
     CD8 APC-H7), then the other antibodies on that dye in alphabetical order;
  2. the detector: the channel the control labels in $PnS, else the channel whose name carries
     the dye (NHLBI detectors are named by laser and filter);
  3. the split: a two-component Gaussian mixture on arcsinh(x / 150), as code/stain_index.py and
     the thresholds in pilot.py do. A cell tube holds a small positive fraction and the
     two-component cut sits inside the negative, so on cells the positive is the brightest of
     three components and the negative stays the lower of two. When the smaller component holds
     under 10 percent of events the tube has no negative population and the unstained tube of
     the center is the negative;
  4. on the linear scale:
       stain index = (median positive - median negative) / (2 * robust SD negative),
       robust SD   = 1.4826 * median absolute deviation,
       separation  = (median positive - median negative) / median negative.

The index is a property of one conjugate on one instrument at one setting. On beads it carries
the conjugate (its fluorochrome-to-protein ratio and the amount the bead captures) and the
detector gain, and not the antigen. The control antibodies differ by center (Baylor has no CD8
control and takes CD20 APC-H7; Stanford has no CD45RA control and takes CD194 PE-Cy7), so part
of the between-center difference on a dye is the antibody and not the instrument. The FITC
value is a Live/Dead stain on cells and is in the table for completeness; Live is not a feature
of the clustering.

Usage: python3 -I code/hipc_stain_index.py <project_dir>
Writes results/hipc_stain_index.csv (one row per center and dye, with the pooled value and the
fallback flag) and results/hipc_stain_index.md.
"""
import re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import flowio
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")
COFACTOR = 150.0
MIN_COMPONENT = 0.10
CENTERS = ("Baylor", "CIMR", "Miami", "NHLBI", "Stanford", "UCLA", "Yale")

# Dye names as pilot.HIPC_TCELL_FLUOR spells them, so the loader joins on them directly. The
# patterns match file names and channel names; tandems are tested before their parent dyes.
DYES = [
    ("PE-Cy7", r"PE[- ]?Cy7|G780"),
    ("APC-Cy7", r"APC[- ]?(?:H7|Cy7)|R780"),
    ("PerCP-Cy5.5", r"PerCP|B710"),
    ("APC", r"\bAPC\b|AF647|Alexa[- ]647|R660"),
    ("PE", r"\bPE\b|G560"),
    ("FITC", r"FITC|Alexa Fluor 488|Live[ _,2f]*Green|LIVE[ _]GREEN|Live Green|B515"),
    ("Pacific Blue", r"V450|Pacific Blue"),
    ("AmCyan", r"V500|Am ?Cyan|V545"),
]
ALIASES = {"PE-Cy7": "PE-Cy7", "APC-Cy7": "APC-H7 / APC-Cy7", "PerCP-Cy5.5": "PerCP-Cy5.5", "APC": "APC", "PE": "PE",
           "FITC": "FITC", "Pacific Blue": "V450 / Pacific Blue", "AmCyan": "V500 / AmCyan"}
# the panel's own reagent per dye, preferred when the center has it
PANEL_MARKER = {"FITC": "Live", "PE": "CCR7", "PerCP-Cy5.5": "CD4", "PE-Cy7": "CD45RA", "Pacific Blue": "CD3",
                "AmCyan": "HLA-DR", "APC": "CD38", "APC-Cy7": "CD8"}
MARKER_RE = re.compile(r"^(CD[0-9A-Za-z+_]+|CCR\d|CXCR\d|HLA[- ]?DR|Lineage|LIVE[ _]?(?:DEAD|GREEN)?|Live\S*(?: Dead| Green)?|IGD)", re.I)
MARKER_NORM = {"CD197": "CCR7", "HLA DR": "HLA-DR", "HLADR": "HLA-DR", "CD3+19+20": "CD3+CD19+CD20", "CD3_CD19_CD20": "CD3+CD19+CD20",
               "CD3+CD19+CD20+": "CD3+CD19+CD20", "Lineage": "CD3+CD19+CD20"}


def dye_of(text):
    for dye, pat in DYES:
        if re.search(pat, text, re.I): return dye
    return None


def marker_of(stem):
    m = MARKER_RE.match(stem)
    if not m: return ""
    mk = m.group(1)
    if re.match(r"live", mk, re.I): return "Live"
    return MARKER_NORM.get(mk, mk)


def read_fcs(path):
    fd = flowio.FlowData(str(path), ignore_offset_error=True)
    keys = sorted(fd.channels, key=int)
    names = [fd.channels[k]["pnn"] for k in keys]
    stains = [fd.channels[k].get("pns", "") or "" for k in keys]
    volts = [fd.text.get(f"p{k}v", "") for k in keys]
    X = np.reshape(fd.events, (-1, fd.channel_count)).astype(np.float64)
    return dict(path=Path(path), names=names, stains=stains, volts=volts, X=X, text=fd.text)


def robust_sd(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def split_pos(x, k=2, seed=0):
    """Events in the brightest of k Gaussian components on arcsinh(x / COFACTOR)."""
    v = np.arcsinh(x / COFACTOR).reshape(-1, 1)
    g = GaussianMixture(k, random_state=seed, n_init=3).fit(v)
    hi = int(np.argmax(g.means_.ravel()))
    return g.predict(v) == hi


def primary_channel(f, dye):
    labeled = [n for n, s in zip(f["names"], f["stains"]) if s]
    if labeled:
        if dye_of(labeled[0]) != dye:
            raise SystemExit(f"{f['path'].name}: labeled channel {labeled[0]} is not the {dye} detector")
        return labeled[0]
    cand = [n for n in f["names"] if n not in ("Time",) and not n.startswith(("FSC", "SSC")) and dye_of(n) == dye]
    if len(cand) != 1:
        raise SystemExit(f"{f['path'].name}: {len(cand)} channels match {dye}: {cand}")
    return cand[0]


def control_stem(path):
    s = path.stem
    s = re.sub(r"^(Compensation Controls_|Specimen_001_)", "", s)
    s = re.sub(r"\s*Stained Control(_[A-H]\d\d)?$", "", s)
    return s


def candidates(cdir):
    out = []
    for p in sorted(cdir.glob("*.fcs")):
        stem = control_stem(p)
        if re.search(r"unstained|Tube_001|^\d+$", stem, re.I): continue
        dye = dye_of(stem)
        if dye is None: continue
        out.append(dict(file=p, stem=stem, dye=dye, marker=marker_of(stem)))
    return out


def unstained_file(cdir):
    u = sorted(cdir.glob("Compensation Controls_Unstained Control*.fcs"))
    return u[0] if u else None


def measure(f, ch, unst):
    j = f["names"].index(ch)
    x = f["X"][:, j]
    pos = split_pos(x)
    frac_small = min(pos.mean(), 1 - pos.mean())
    fsc = f["X"][:, f["names"].index("FSC-A")]
    kind = "beads" if robust_sd(fsc) / max(np.median(fsc), 1) < 0.2 else "cells"
    if kind == "cells" and frac_small >= MIN_COMPONENT:
        # a stained cell tube holds a small positive fraction, and the two-component cut sits
        # inside the negative; the positive is the brightest of three components, the negative
        # stays the lower of two
        top = split_pos(x, k=3)
        if top.mean() >= 0.01 and np.median(x[top]) > np.median(x[pos]):
            neg = x[~pos]; pos = top; source = "in-tube, 3 components"
        else:
            neg = x[~pos]; source = "in-tube"
    elif frac_small >= MIN_COMPONENT:
        neg = x[~pos]; source = "in-tube"
    elif unst is not None:
        neg = unst["X"][:, unst["names"].index(ch)]; source = "unstained tube"
        pos = np.arcsinh(x / COFACTOR) > np.median(np.arcsinh(neg / COFACTOR)) + 3 * robust_sd(np.arcsinh(neg / COFACTOR))
    else:
        return None
    xp = x[pos]
    med_p, med_n, sd_n = np.median(xp), np.median(neg), robust_sd(neg)
    return dict(detector=ch, volts=f["volts"][j], control_type=kind, n_events=len(x), n_positive=int(pos.sum()),
                positive_fraction=float(pos.mean()), negative=source,
                median_negative=med_n, robust_sd_negative=sd_n, median_positive=med_p,
                stain_index=(med_p - med_n) / (2 * sd_n) if sd_n > 0 else np.nan,
                separation=(med_p - med_n) / med_n if med_n > 0 else np.nan)


def main(root):
    root = Path(root); res = root / "results"; res.mkdir(exist_ok=True)
    comp = root / "data/HIPC/COMP"
    dyes = [d for d, _ in DYES]
    rows, skipped = [], []
    for c in CENTERS:
        cdir = comp / c
        cands = candidates(cdir)
        uf = unstained_file(cdir); unst = read_fcs(uf) if uf else None
        for dye in dyes:
            cs = [k for k in cands if k["dye"] == dye]
            cs.sort(key=lambda k: (k["marker"] != PANEL_MARKER[dye], k["marker"] == "", k["marker"], k["file"].name))
            for k in cs[1:]:
                skipped.append(dict(center=c, dye=dye, file=k["file"].name, marker=k["marker"]))
            if not cs:
                rows.append(dict(center=c, dye=dye, marker="", file="", n_events=0, stain_index=np.nan, separation=np.nan,
                                 note="no control for this dye at this center"))
                continue
            k = cs[0]; f = read_fcs(k["file"])
            ch = primary_channel(f, dye)
            m = measure(f, ch, unst)
            if m is None:
                rows.append(dict(center=c, dye=dye, marker=k["marker"], file=k["file"].name, n_events=len(f["X"]), stain_index=np.nan,
                                 separation=np.nan, note="one population in the tube and no unstained tube"))
                continue
            rows.append(dict(center=c, dye=dye, marker=k["marker"], file=k["file"].name, **m, note=""))
            print(f"{c:9s} {dye:12s} {k['marker']:16s} {ch:18s} n={m['n_events']:5d} pos={m['positive_fraction']:.2f} "
                  f"{m['control_type']:5s} {m['negative']:14s} SI={m['stain_index']:7.1f}", flush=True)
    df = pd.DataFrame(rows)
    # pooled value per dye: the arithmetic mean over the centers that have a usable control
    pooled = df.groupby("dye").stain_index.mean()
    df["pooled_stain_index"] = df.dye.map(pooled)
    df["uses_pooled"] = df.stain_index.isna()
    df["stain_index_used"] = np.where(df.uses_pooled, df.pooled_stain_index, df.stain_index)
    df["cytometer"] = df.center.map({c: cyt_of(comp / c) for c in CENTERS})
    cols = ["center", "cytometer", "dye", "marker", "file", "detector", "volts", "control_type", "n_events", "n_positive", "positive_fraction",
            "negative", "median_negative", "robust_sd_negative", "median_positive", "stain_index", "separation",
            "pooled_stain_index", "uses_pooled", "stain_index_used", "note"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(res / "hipc_stain_index.csv", index=False)
    pd.DataFrame(skipped).to_csv(res / "hipc_stain_index_skipped.csv", index=False)
    write_md(df, pd.DataFrame(skipped), res / "hipc_stain_index.md")
    print(matrix_md(df))


def cyt_of(cdir):
    uf = unstained_file(cdir) or next(iter(sorted(cdir.glob("*.fcs"))))
    t = read_fcs(uf)["text"]
    return f"{t.get('cyt', '')} {t.get('cytnum', '')}".strip()


def f0(v): return "" if pd.isna(v) else f"{v:,.0f}"


def matrix_md(df):
    """Dye by center: stain index, with the control marker, bold when the value is pooled."""
    dyes = [d for d, _ in DYES]
    order = ["FITC", "PE", "PerCP-Cy5.5", "PE-Cy7", "Pacific Blue", "AmCyan", "APC", "APC-Cy7"]
    L = ["| Dye | " + " | ".join(CENTERS) + " | Mean over centers | Brightest over dimmest |", "|---|" + "---|" * (len(CENTERS) + 2)]
    for dye in order:
        cells = []
        for c in CENTERS:
            r = df[(df.center == c) & (df.dye == dye)].iloc[0]
            if pd.isna(r.stain_index): cells.append("none, pooled " + f0(r.pooled_stain_index))
            else: cells.append(f"{f0(r.stain_index)} ({r.marker or 'no marker named'}{', cells' if r.control_type == 'cells' else ''})")
        si = df[(df.dye == dye)].stain_index.dropna()
        L.append(f"| {ALIASES[dye]} | " + " | ".join(cells) + f" | {f0(si.mean())} | {si.max() / si.min():.1f} |")
    return "\n".join(L)


def write_md(df, skipped, path):
    L = ["# Stain index per center from the HIPC Lyoplate single-stain controls\n",
         "**Run:** 2026-10-19 · **Code:** `code/hipc_stain_index.py` · **Data:** `data/HIPC/COMP/<Center>/`, the raw FACSDiva "
         "compensation controls of the seven centers with central gates · **Output:** `results/hipc_stain_index.csv` "
         "(one row per center and dye), `results/hipc_stain_index_skipped.csv` (the other controls on the same dye)\n",
         "## What ran\n",
         "The antibody controls are BD CompBeads, about half positive and half negative beads in one tube, 5,000 to 11,000 events. "
         "The Live/Dead FITC controls are cells. For each center and each of the eight T cell panel dyes the script takes the "
         "control that carries the panel's own reagent when the center has it, else another antibody on the same dye, and "
         "splits the tube on its primary detector with a two-component Gaussian mixture on arcsinh(x / 150). The stain index "
         "is the median positive minus the median negative over twice the robust SD of the negative, on the linear scale, "
         "uncompensated. On the cell tubes the positive is the brightest of three components, because the two-component cut "
         "sits inside the negative when the positive fraction is small. When a tube holds one population, the unstained tube "
         "of the center is the negative.\n",
         "The index is a property of one conjugate on one instrument at one setting. On beads it carries the conjugate and the "
         "detector gain, and not the antigen. The control antibodies differ by center, so part of a between-center difference "
         "on a dye is the antibody and not the instrument. The FITC value is a Live/Dead stain on cells; Live is not a feature "
         "of the clustering and the value is here for completeness.\n",
         "## Stain index by dye and center\n",
         "The control marker is in parentheses. \"none, pooled\" marks a center without a control for that dye; arm B_site takes "
         "the mean over the other centers there.\n",
         matrix_md(df), "",
         "## Every control used\n",
         "| Center | Cytometer | Dye | Marker | File | Detector | Volts | Control | Events | Positive fraction | Negative | Median neg | Robust SD neg | Median pos | Stain index | Separation | Note |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        L.append(f"| {r.center} | {r.cytometer} | {ALIASES[r.dye]} | {r.marker} | {r.file} | {r.get('detector', '') or ''} | {r.get('volts', '') or ''} | "
                 f"{r.get('control_type', '') or ''} | {r.n_events} | {'' if pd.isna(r.get('positive_fraction')) else f'{r.positive_fraction:.2f}'} | "
                 f"{r.get('negative', '') or ''} | {f0(r.get('median_negative'))} | {f0(r.get('robust_sd_negative'))} | {f0(r.get('median_positive'))} | "
                 f"{f0(r.stain_index)} | {'' if pd.isna(r.separation) else f'{r.separation:,.0f}'} | {r.note} |")
    L += ["", "## Controls on the same dye that the script did not use\n",
          "| Center | Dye | Marker | File |", "|---|---|---|---|"]
    for _, r in skipped.iterrows():
        L.append(f"| {r.center} | {ALIASES[r.dye]} | {r.marker} | {r.file} |")
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
