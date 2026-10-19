"""Measured stain index per fluorochrome from single-stain compensation controls.

The pilot divides each channel by a vendor brightness rank (pilot.py BRIGHTNESS, ordinal 1 to 5).
Section 4 of plan-and-datasets.md names a stain index measured from the dataset's own single
stains as the next source up. This script computes it and compares the two.

Method, per single-stain control, on the primary channel:
  1. the positive beads are the events of the control tube above the unstained level,
     threshold = unstained median + 3 robust SD on arcsinh(x / cofactor), cofactor 150;
  2. the negative is the unstained control tube on the same channel (every control tube in
     this run holds 95 to 99 percent positive beads, so the tube has no negative population;
     the script still fits a two-component Gaussian mixture per tube and reports the in-tube
     stain index when the smaller component holds at least 10 percent of events);
  3. on the linear scale, uncompensated:
       stain index = (median positive - median negative) / (2 * robust SD of the negative),
       robust SD = 1.4826 * median absolute deviation,
       separation  = (median positive - median negative) / median negative.
  A single-stain control carries one dye, so its own detector needs no compensation; the
  uncompensated value is the measurement. The compensated value (raw @ inv(S), S from --spill,
  rows fluorochromes, columns detectors, FCS 3.1 $SPILLOVER convention) is a check column.

Consistency check: FACSDiva writes an identity SPILL into each control file and the experiment
matrix into the stained sample files. For each control the script compares the observed spill
pattern (median of the tube minus median of the unstained tube, in every detector, over the same
in the primary detector) with the stored row of S and reports the largest difference. In the
December 2024 run six rows differ by more than 0.5, so the stored matrix did not come from these
tubes as exported, and compensating them with it is wrong. The report flags those tubes.

Usage:
  python3 -I stain_index.py <project_dir> --spill <stained sample.fcs> --unstained <unstained.fcs>
                            [--cofactor 150] <control.fcs> ...
Writes results/stain_index.csv and results/stain_index.md. Reads the FCS files in place and
copies nothing.
"""
import argparse, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import flowio
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pilot import BRIGHTNESS, BRIGHT_MAPS   # noqa: E402

warnings.filterwarnings("ignore")

# $PnS fluorochrome label in the sample files -> key in pilot.py BRIGHTNESS
FLUOR_KEYS = {"APC-A750": "APC-Alexa Fluor 750", "FITC": "FITC", "PE": "PE", "APC": "APC",
              "A700": "Alexa Fluor 700", "ECD": "PE-Texas Red", "KrO": "Krome Orange", "PE-Cy7": "PE-Cy7",
              "BUV737": "BUV737", "BV786": "BV786", "BUV661": "BUV661", "Violet": "Viability dye, violet"}
VIABILITY = {"Viability dye, violet"}


def read_fcs(path):
    fd = flowio.FlowData(str(path), ignore_offset_error=True)
    names = [fd.channels[k]["pnn"] for k in sorted(fd.channels, key=int)]
    stains = [fd.channels[k].get("pns", "") or "" for k in sorted(fd.channels, key=int)]
    X = np.reshape(fd.events, (-1, fd.channel_count)).astype(np.float64)
    sp = fd.text.get("spill") or fd.text.get("$spillover") or fd.text.get("spillover")
    S, order = None, None
    if sp:
        parts = sp.split(","); n = int(parts[0])
        order = parts[1:1 + n]
        S = np.array(parts[1 + n:], dtype=float).reshape(n, n)
    return dict(path=Path(path), text=fd.text, names=names, stains=stains, X=X, spill=S, spill_order=order)


def compensate(f, S, order):
    cols = [f["names"].index(c) for c in order]
    f["comp"] = f["X"][:, cols] @ np.linalg.inv(S)
    return f


def robust_sd(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def split_pos_neg(x, cofactor, seed=0):
    v = np.arcsinh(x / cofactor).reshape(-1, 1)
    g = GaussianMixture(2, random_state=seed, n_init=3).fit(v)
    hi = int(np.argmax(g.means_.ravel()))
    pos = g.predict(v) == hi
    return pos


def primary_channel(f, order):
    tube = f["text"].get("tube name", "")
    base = re.sub(r"\s*Stained Control$", "", tube).strip()
    cand = base + "-A"
    if cand not in order:
        raise SystemExit(f"{f['path'].name}: tube '{tube}' matches no SPILL channel in {order}")
    return cand


def marker_map(sample):
    """{channel $PnN: (marker, brightness key)} from the stained sample's $PnS labels."""
    out = {}
    for n, s in zip(sample["names"], sample["stains"]):
        m = re.match(r"^(.*?)\s*\((.+)\)\s*$", s or "")
        if m:
            out[n] = (m.group(1).strip(), FLUOR_KEYS.get(m.group(2).strip(), m.group(2).strip()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project_dir")
    ap.add_argument("--spill", required=True, help="stained sample FCS that carries the experiment SPILL and $PnS")
    ap.add_argument("--unstained", required=True)
    ap.add_argument("--cofactor", type=float, default=150.0)
    ap.add_argument("controls", nargs="+")
    a = ap.parse_args()
    root = Path(a.project_dir); res = root / "results"; res.mkdir(exist_ok=True)

    sample = read_fcs(a.spill)
    S, order = sample["spill"], sample["spill_order"]
    if S is None or np.allclose(S, np.eye(len(order))):
        raise SystemExit("--spill file carries no experiment matrix")
    mm = marker_map(sample)
    unst = compensate(read_fcs(a.unstained), S, order)
    ucols = [unst["names"].index(o) for o in order]
    uraw = unst["X"][:, ucols]
    umed = np.median(uraw, axis=0)
    rows = []
    for c in sorted(a.controls):
        f = read_fcs(c)
        own_identity = f["spill"] is not None and np.allclose(f["spill"], np.eye(len(f["spill_order"])))
        compensate(f, S, order)
        ch = primary_channel(f, order); j = order.index(ch)
        raw = f["X"][:, [f["names"].index(o) for o in order]]
        x = raw[:, j]; u = uraw[:, j]
        vu = np.arcsinh(u / a.cofactor); thr = np.median(vu) + 3 * robust_sd(vu)
        pos = np.arcsinh(x / a.cofactor) > thr
        xp = x[pos]
        med_p, med_u, sd_u = np.median(xp), np.median(u), robust_sd(u)
        # in-tube split, reported only when the tube really holds two populations
        g = split_pos_neg(x, a.cofactor)
        frac_low = 1 - g.mean()
        if min(frac_low, 1 - frac_low) >= 0.10:
            xn = x[~g]; med_n, sd_n = np.median(xn), robust_sd(xn)
            si_tube = (np.median(x[g]) - med_n) / (2 * sd_n)
        else:
            med_n = sd_n = si_tube = np.nan
        # observed spill pattern of this tube against the stored matrix row
        obs = (np.median(raw[pos], axis=0) - umed); obs = obs / obs[j]
        dev = np.abs(obs - S[j]); dev[j] = 0
        k = int(np.argmax(dev))
        xc = f["comp"][:, j]; uc = unst["comp"][:, j]
        si_comp = (np.median(xc[pos]) - np.median(uc)) / (2 * robust_sd(uc))
        marker, key = mm.get(ch, ("", ch))
        rank, src = BRIGHTNESS.get(key, (np.nan, "missing from BRIGHTNESS"))
        rows.append(dict(
            file=f["path"].name, channel=ch, fluorochrome=key, marker=marker,
            n_events=len(x), n_positive=int(pos.sum()),
            median_negative=med_u, robust_sd_negative=sd_u, median_positive=med_p,
            stain_index=(med_p - med_u) / (2 * sd_u),
            separation=(med_p - med_u) / med_u if med_u > 0 else np.nan,
            stain_index_compensated=si_comp,
            in_tube_low_fraction=frac_low, in_tube_median_negative=med_n, in_tube_sd_negative=sd_n, in_tube_stain_index=si_tube,
            max_spill_observed=float(np.delete(obs, j).max()), max_spill_into=order[int(np.argmax(np.where(np.arange(12) == j, -np.inf, obs)))],
            spill_deviation=float(dev[k]), spill_deviation_detector=order[k], stored_spill_there=float(S[j, k]), observed_spill_there=float(obs[k]),
            own_spill_identity=own_identity, voltage=f["text"].get(f"p{f['names'].index(ch) + 1}v"),
            vendor_rank=rank, vendor_source=src))
    df = pd.DataFrame(rows).sort_values("stain_index", ascending=False).reset_index(drop=True)
    df.to_csv(res / "stain_index.csv", index=False)
    write_md(df, res / "stain_index.md", a, sample)
    pd.set_option("display.width", 250)
    print(df[["fluorochrome", "marker", "n_positive", "median_negative", "robust_sd_negative", "median_positive",
              "stain_index", "stain_index_compensated", "in_tube_low_fraction", "in_tube_stain_index",
              "max_spill_observed", "max_spill_into", "spill_deviation", "spill_deviation_detector", "vendor_rank"]].to_string())


def compare(df):
    """Rank agreement and spread. Four sets: all tubes, the antibody conjugates, the tubes whose
    observed spill pattern matches the stored matrix (largest deviation under 0.5), and the
    antibodies among those."""
    out = {}
    clean = df[df.spill_deviation < 0.5]
    sets = (("all 12 tubes", df), ("11 antibody conjugates", df[~df.fluorochrome.isin(VIABILITY)]),
            ("tubes consistent with the matrix", clean), ("antibodies consistent with the matrix", clean[~clean.fluorochrome.isin(VIABILITY)]))
    for label, sub in sets:
        sub = sub.dropna(subset=["vendor_rank"])
        rho, p = spearmanr(sub.stain_index, sub.vendor_rank)
        si = sub.stain_index
        out[label] = dict(n=len(sub), rho=rho, p=p, members=", ".join(sub.fluorochrome),
                          si_ratio=si.max() / si.min(), si_max=sub.loc[si.idxmax(), "fluorochrome"], si_min=sub.loc[si.idxmin(), "fluorochrome"],
                          rank_max=int(sub.vendor_rank.max()), rank_min=int(sub.vendor_rank.min()),
                          geom_ratio=BRIGHT_MAPS["B_geom"](sub.vendor_rank.max()) / BRIGHT_MAPS["B_geom"](sub.vendor_rank.min()),
                          lin_ratio=BRIGHT_MAPS["B_lin"](sub.vendor_rank.max()) / BRIGHT_MAPS["B_lin"](sub.vendor_rank.min()))
    return out


def geo1(v):
    v = np.asarray(v, float)
    return v / np.exp(np.mean(np.log(v)))


def write_md(df, path, a, sample):
    cmp = compare(df)
    cyt = sample["text"].get("cyt", ""); date = sample["text"].get("date", "")
    L = []
    L.append("# Measured stain index against the vendor brightness rank\n")
    L.append(f"**Run:** 2026-10-18 · **Code:** `code/stain_index.py` · **Data:** the DURAClone IM Basic run of {date} on the {cyt}, "
             "12 compensation controls and one unstained control, read in place from `<DURAClone folder>/` "
             "(permission to publish that set stays open, and the script copied nothing) · **Inventory:** `data/README.md`\n")
    L.append("## What ran\n")
    L.append("Each compensation control is a tube of capture beads, 5,000 events, one antibody conjugate per tube. "
             "The tubes hold 95 to 99 percent positive beads, so a tube has no negative population of its own. "
             "The negative for every channel is the unstained bead tube. The positive beads are the events above the unstained median "
             f"plus 3 robust SD on arcsinh(x / {a.cofactor:.0f}). The stain index is the median positive minus the median negative, "
             "over twice the robust SD of the negative. The script computes it on the linear scale, on the uncompensated primary detector.\n")
    L.append("A single stain carries one dye, so its own detector needs no compensation. "
             "The column \"SI compensated\" applies the stored matrix of the stained sample files as a check. "
             "\"Separation\" is the median positive minus the median negative, over the median negative. "
             "The unstained median sits between 0 and 40 on nine channels, so that ratio is unstable there. "
             "The table keeps it only because the plan asked for it.\n")
    L.append("The two-component Gaussian mixture that the plan named splits the positive mode in two on these tubes, because the tubes have no "
             "negative beads. The script reports the in-tube index only when the smaller component holds at least 10 percent of events. "
             "That happens on three tubes (BUV661, Krome Orange, APC). Those three values are splits inside the positive mode and "
             "are no stain index. They are in `stain_index.csv` and out of the table.\n")
    L.append("## Table\n")
    L.append("Sorted by stain index. The spill pattern of a tube is the median of the tube minus the median of the unstained tube, "
             "in every detector. The script divides that by the same difference in the tube's own detector. \"Spill check\" is the largest difference between that "
             "pattern and the stored SPILL row of the stained sample files, with the detector where it occurs. "
             "A difference above 0.5 means the stored matrix did not come from this tube as exported.\n")
    L.append("| Fluorochrome | Marker | Detector | Volts | Positive beads | Median neg | Robust SD neg | Median pos | Stain index | SI compensated | Separation | Spill check | Vendor rank | Rank source |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in df.iterrows():
        sep = "" if pd.isna(r.separation) else f"{r.separation:,.0f}"
        flag = "**" if r.spill_deviation >= 0.5 else ""
        L.append(f"| {r.fluorochrome} | {r.marker} | {r.channel} | {r.voltage} | {r.n_positive:,} | {r.median_negative:,.0f} | {r.robust_sd_negative:,.0f} | "
                 f"{r.median_positive:,.0f} | {r.stain_index:,.0f} | {r.stain_index_compensated:,.0f} | {sep} | {flag}{r.spill_deviation:.2f} in {r.spill_deviation_detector}{flag} | "
                 f"{r.vendor_rank:.0f} | {r.vendor_source} |")
    L.append("")
    L.append("## Comparison with the vendor rank\n")
    L.append("| Set | n | Spearman rho | p | Brightest / dimmest, measured | Brightest / dimmest, 2^rank | Brightest / dimmest, linear rank |")
    L.append("|---|---|---|---|---|---|---|")
    for label, c in cmp.items():
        L.append(f"| {label} | {c['n']} | {c['rho']:.2f} | {c['p']:.2f} | {c['si_ratio']:,.0f} ({c['si_max']} / {c['si_min']}) | "
                 f"{c['geom_ratio']:.0f} (rank {c['rank_max']} / {c['rank_min']}) | {c['lin_ratio']:.2f} |")
    L.append("")
    L.append("Members of the two consistent sets: " + cmp["tubes consistent with the matrix"]["members"] + ".\n")
    ab = df[~df.fluorochrome.isin(VIABILITY)].copy()
    ab["div_si"] = geo1(ab.stain_index); ab["div_geom"] = geo1(ab.vendor_rank.map(BRIGHT_MAPS["B_geom"])); ab["div_lin"] = geo1(ab.vendor_rank.map(BRIGHT_MAPS["B_lin"]))
    L.append("Divisors with geometric mean 1, as arm B applies them before arcsinh, for the 11 antibody conjugates. Arm B divides the "
             "channel by the divisor. A value above 1 pulls a bright channel down, and a value below 1 lifts a dim channel up.\n")
    L.append("| Fluorochrome | Marker | Measured SI | Divisor, measured SI | Divisor, 2^rank | Divisor, linear rank | Spill check |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in ab.sort_values("stain_index", ascending=False).iterrows():
        L.append(f"| {r.fluorochrome} | {r.marker} | {r.stain_index:,.0f} | {r.div_si:.2f} | {r.div_geom:.2f} | {r.div_lin:.2f} | {'fails' if r.spill_deviation >= 0.5 else 'passes'} |")
    L.append("")
    L.append(VERDICT)
    path.write_text("\n".join(L))


VERDICT = """## What the numbers say

The measured stain index does not agree with the vendor rank on this run. Spearman rho is
-0.19 on the 11 antibody conjugates and 0.00 on the five antibodies whose tubes pass the
spill check. The vendor table puts PE, APC and PE-Cy7 at rank 5 above everything else. The
beads put BUV737 first at 4,265, then BUV661, BV786 and FITC, with PE sixth at 964, PE-Cy7
tenth at 481 and APC last at 83.

The measured spread is also far wider than either mapping. Brightest over dimmest is 52 on the 11 antibodies and 15 on the five that pass the spill
check. The 2^rank mapping spans 4 and the linear mapping spans 1.67.

The 2^rank mapping in pilot.py is therefore not close to the measured spread, and the
direction of the divisors is also different. On this panel the rank mapping divides PE, APC
and PE-Cy7 by 2.13 and lifts FITC, Krome Orange and the APC tandem by 0.53. The measured
divisors divide BUV737 by 5.14 and the three Brilliant dyes by 2.2 to 5.1, leave PE at 1.16,
and lift PE-Cy7 to 0.58 and APC to 0.10. A paper that runs arm B with vendor ranks and a paper
that runs it with these single stains make opposite moves on three of the eight DURAClone
channels.

Part of the gap is the instrument and not the dye. The stain index of a bead is the dye
times the conjugate times the detector setting. The BUV737 detector ran at 687 volts and the
PE detector at 410, and on the LSRFortessa the UV detectors need the high voltage, so the
BUV737 index carries the gain. That is the correct quantity for arm B, because the stained
samples carry the same gain, and arm B is meant to remove the scale the acquisition put on
each channel. The vendor rank cannot see the voltage at all.

Part of the gap is this run. Six tubes show a spill pattern that the stored SPILL matrix does
not describe (spill check above 0.5): the ECD tube shows 2.3 times more signal in the PE
detector than in its own, which is the pattern of a PE conjugate or of a broken PE tandem,
and the APC, Alexa Fluor 700, APC-Alexa Fluor 750, BV786 and BUV661 tubes show 0.9 to 3.6
times their own signal in the BUV737 detector, where the stored matrix says 0.05 to 0.27.
Those five all emit in the far red, and a UV detector that reads a far-red emitter at that
level points to a laser delay or window setting, which the files do not record. The stored
matrix came from some other acquisition of the controls, and this export cannot say which one
is right. The ECD, APC and BUV737 values
in the table are therefore the least trustworthy, and the five-antibody set is the one to
quote. On that set the rank correlation is still zero and the spread is still 15.

## Divisor rule for arm B when single stains exist

Use the measured stain index as the divisor on the linear scale before arcsinh, with the
geometric mean of the divisors set to 1 so that the cofactor keeps its meaning:

    divisor_j = SI_j / geomean(SI)        x'_j = x_j / divisor_j        arcsinh(x'_j / cofactor)

Compute SI_j on the uncompensated primary detector of each single-stain control, with the
unstained tube as the negative and the same positive threshold as this script. Run the spill
check first, and drop a control whose pattern fails it; for that channel fall back to the
vendor rank and say so in the table. Report the rank correlation between the measured index
and the vendor rank as a quality figure, because a value near zero, as here, says the rank
mapping was not a usable stand-in for that run. Fall back to the 2^rank mapping only when no
controls exist, and say in the Methods that the rank mapping spans a factor of 4 where
measured values span a factor of 15 to 50.

## Limits

A stain index is a property of one conjugate on one instrument at one setting, and it is not
a property of the dye alone. Each control here is a different antibody, so the index mixes
the dye brightness with the fluorochrome-to-protein ratio of that lot and with the amount of
antibody the capture bead binds. On cells the index also mixes in the antigen density, and the
stained samples of this run are not single stains, so this run cannot separate the two. A
divisor from bead controls corrects for the conjugate and the detector setting, and a divisor
from cell controls corrects for the antigen density as well. The paper has to say which one it
uses and why, and on a cell control the Krome Orange CD45 channel, the dimmest conjugate
here, resolves well because CD45 is among the densest antigens on a leukocyte.

Two smaller limits. The viability dye is not an antibody, and its bead control (BV480
detector, index 790) is a stain of the bead surface, so it is in the table and out of the
antibody comparisons. And the unstained tube reads 568 on the APC detector with a robust SD of
62, against 0 to 40 on nine other detectors, so the APC index of 83 rests on a negative that
is itself unusual.
"""


if __name__ == "__main__":
    main()
