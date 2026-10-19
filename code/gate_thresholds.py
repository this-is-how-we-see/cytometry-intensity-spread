"""Supervised per-marker thresholds from the manual labels of Levine_32dim.

This is the upper bound the plan names in section 5, method 1. It is circular by construction:
the labels define which populations are positive and negative for a marker, and the cut on
arcsinh(x / 5) is the one that best separates those populations. The paper reports it as the
upper bound and never as the recommended method.

For each marker the truth table below names the populations the label set defines as positive
and as negative. A population with intermediate expression (CD4 on monocytes, CD7 on T cells,
CD20 on pre B cells) is left out of both lists on purpose, so the cut is not pulled toward it.
A marker that no labeled population defines (erythroid CD235ab, platelet CD41 and CD61, the
continuous markers CXCR4, CD44, CD47, CD49d, CD321, and CD15, whose granulocytes are not in
the label set) falls back to the two-component Gaussian mixture threshold from pilot.py.
CD45 falls back as well, because every labeled population is positive.

The cut maximizes balanced accuracy, 0.5 * (true positive rate + true negative rate), over a
grid of 2,000 candidate values between the 0.1st and 99.9th percentile of the labeled cells.
Balanced accuracy is used so that a 300-cell population and a 26,000-cell population count
the same.

Usage: python3 gate_thresholds.py <project_dir>
Writes results/gate_thresholds.json.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pilot import load_levine32, gmm_threshold

T = ["CD4_T_cells", "CD8_T_cells"]
NK = ["CD16-_NK_cells", "CD16+_NK_cells"]
HSPC = ["CD34+CD38+CD123-_HSPCs", "CD34+CD38+CD123+_HSPCs"]
HSC = "CD34+CD38lo_HSCs"
B = ["Mature_B_cells", "Pre_B_cells", "Pro_B_cells"]

# marker: (positive populations, negative populations, note). Populations use the names in
# population_names_Levine_32dim.txt. An empty positive or negative list means fallback.
TRUTH = {
    "CD45RA": (["Mature_B_cells"] + NK, [HSC, "CD34+CD38+CD123-_HSPCs", "Basophils", "Pro_B_cells"],
               "T cells left out: naive and memory mix"),
    "CD133": ([HSC, "CD34+CD38+CD123-_HSPCs"], T + ["Mature_B_cells"] + NK + ["Monocytes"],
              "stain is dim on HSCs in this file, median arcsinh 0.2"),
    "CD19": (B + ["Plasma_B_cells"], T + NK + ["Monocytes", HSC, "Basophils", "pDCs"], ""),
    "CD22": (["Mature_B_cells"], T + NK + ["Monocytes", HSC, "pDCs"], "pre B and pro B left out: dim"),
    "CD11b": (["Monocytes"], T + B + [HSC] + HSPC + ["pDCs"], "NK cells and basophils left out: dim"),
    "CD4": (["CD4_T_cells"], ["CD8_T_cells", "Mature_B_cells"] + NK + [HSC, "Pre_B_cells"],
            "monocytes and pDCs left out: intermediate, the hypothesis 3 case"),
    "CD8": (["CD8_T_cells"], ["CD4_T_cells", "Mature_B_cells", "Monocytes", HSC, "Pre_B_cells", "pDCs"],
            "NK cells left out: CD8 dim"),
    "CD34": ([HSC] + HSPC + ["Pro_B_cells"], T + ["Mature_B_cells"] + NK + ["Monocytes", "Basophils", "pDCs", "Plasma_B_cells", "Pre_B_cells"], ""),
    "Flt3": (["CD34+CD38+CD123+_HSPCs", "pDCs"], T + NK + ["Mature_B_cells"], "stain is dim, median arcsinh 0.3 on the positive populations"),
    "CD20": (["Mature_B_cells"], T + NK + ["Monocytes", HSC, "Pro_B_cells", "Plasma_B_cells"], "pre B left out: dim"),
    "CXCR4": ([], [], "continuous marker, no population in the label set defines it"),
    "CD235ab": ([], [], "erythroid marker, no erythroid population in the label set"),
    "CD45": ([], [], "every labeled population is positive"),
    "CD123": (["pDCs", "Basophils", "CD34+CD38+CD123+_HSPCs"], T + ["Mature_B_cells"] + NK + ["CD34+CD38+CD123-_HSPCs", HSC],
              "monocytes left out: dim"),
    "CD321": ([], [], "continuous marker, no population in the label set defines it"),
    "CD14": (["Monocytes"], T + ["Mature_B_cells"] + NK + [HSC, "pDCs", "Basophils"], "stain is dim on monocytes in this file, median arcsinh 0.1"),
    "CD33": (["Monocytes"], T + ["Mature_B_cells"] + NK + ["Plasma_B_cells"], "basophils left out: median arcsinh 0.1 in this file"),
    "CD47": ([], [], "continuous marker, no population in the label set defines it"),
    "CD11c": (["Monocytes"], T + B + [HSC, "pDCs"], ""),
    "CD7": (NK, ["Mature_B_cells", "Monocytes", "Basophils", "pDCs", HSC, "Plasma_B_cells"], "T cells left out: CD7 dim in this file"),
    "CD15": ([], [], "granulocyte marker, no granulocyte population in the label set"),
    "CD16": (["CD16+_NK_cells"], ["CD16-_NK_cells"] + T + ["Mature_B_cells", HSC, "pDCs", "Basophils"], "monocytes left out: classical and non-classical mix"),
    "CD44": ([], [], "continuous marker, no population in the label set defines it"),
    "CD38": (HSPC + ["Plasma_B_cells", "Pre_B_cells", "Pro_B_cells"], [HSC], "the label set defines CD38 only through the CD38+ against CD38lo split of the CD34+ cells"),
    "CD13": (["Monocytes"], T + ["Mature_B_cells"] + NK, "stain is dim on monocytes, median arcsinh 0.6"),
    "CD3": (T, NK + B + ["Plasma_B_cells", "Monocytes", "Basophils", "pDCs", HSC] + HSPC, ""),
    "CD61": ([], [], "platelet marker, no megakaryocyte population in the label set"),
    "CD117": ([HSC, "CD34+CD38+CD123-_HSPCs"], T + ["Mature_B_cells"] + NK + ["Monocytes", "Pre_B_cells", "Plasma_B_cells"], ""),
    "CD49d": ([], [], "continuous marker, no population in the label set defines it"),
    "HLA-DR": (B + ["Monocytes", "pDCs"] + HSPC, T + NK + ["Basophils", "Plasma_B_cells"], "HSCs left out: intermediate"),
    "CD64": (["Monocytes"], T + ["Mature_B_cells"] + NK + ["pDCs", "Basophils", HSC], ""),
    "CD41": ([], [], "platelet marker, no megakaryocyte population in the label set"),
}


def balanced_accuracy_cut(pos, neg):
    """Cut on a grid that maximizes 0.5 * (TPR + TNR). Ties go to the middle of the tied run."""
    both = np.concatenate([pos, neg])
    lo, hi = np.percentile(both, [0.1, 99.9])
    grid = np.linspace(lo, hi, 2000)
    pos_s, neg_s = np.sort(pos), np.sort(neg)
    tpr = 1 - np.searchsorted(pos_s, grid, side="right") / len(pos_s)
    tnr = np.searchsorted(neg_s, grid, side="right") / len(neg_s)
    ba = 0.5 * (tpr + tnr)
    best = np.flatnonzero(ba >= ba.max() - 1e-12)
    i = best[len(best) // 2]
    return float(grid[i]), float(ba[i]), float(tpr[i]), float(tnr[i])


def balanced_accuracy_at(cut, pos, neg):
    tpr = float(np.mean(pos > cut)); tnr = float(np.mean(neg <= cut))
    return 0.5 * (tpr + tnr), tpr, tnr


def main(root):
    root = Path(root)
    ds = load_levine32(root)
    A1 = np.arcsinh(ds["X"] / ds["cofactor"])
    lab = ds["labels"]
    name_to_id = {v: k for k, v in ds["pop_names"].items()}
    out = {}
    for j, m in enumerate(ds["markers"]):
        pos_pops, neg_pops, note = TRUTH[m]
        gmm = gmm_threshold(A1[:, j])
        rec = dict(mixture_threshold=round(gmm, 4), note=note, positive_populations=pos_pops, negative_populations=neg_pops)
        if pos_pops and neg_pops:
            pos = A1[np.isin(lab, [name_to_id[p] for p in pos_pops]), j]
            neg = A1[np.isin(lab, [name_to_id[p] for p in neg_pops]), j]
            cut, ba, tpr, tnr = balanced_accuracy_cut(pos, neg)
            ba_g, tpr_g, tnr_g = balanced_accuracy_at(gmm, pos, neg)
            rec.update(method="gate", threshold=round(cut, 4), balanced_accuracy=round(ba, 4), tpr=round(tpr, 4), tnr=round(tnr, 4),
                       n_positive_cells=int(len(pos)), n_negative_cells=int(len(neg)),
                       mixture_balanced_accuracy=round(ba_g, 4), mixture_tpr=round(tpr_g, 4), mixture_tnr=round(tnr_g, 4))
        else:
            rec.update(method="mixture_fallback", threshold=round(gmm, 4))
        rec["fraction_positive_all_cells"] = round(float(np.mean(A1[:, j] > rec["threshold"])), 4)
        rec["fraction_positive_mixture"] = round(float(np.mean(A1[:, j] > gmm)), 4)
        out[m] = rec
        print(f"{m:8s} {rec['method']:17s} thr={rec['threshold']:6.3f} gmm={gmm:6.3f}"
              + (f"  BA={rec['balanced_accuracy']:.3f} (mixture {rec['mixture_balanced_accuracy']:.3f})" if "balanced_accuracy" in rec else ""), flush=True)
    meta = dict(dataset="Levine_32dim", transform="arcsinh(x / 5)", objective="balanced accuracy on the labeled cells of the named populations",
                grid="2000 points between the 0.1st and 99.9th percentile", labeled_cells=int((lab > 0).sum()),
                n_gate=sum(r["method"] == "gate" for r in out.values()), n_fallback=sum(r["method"] != "gate" for r in out.values()))
    json.dump(dict(meta=meta, markers=out), open(root / "results/gate_thresholds.json", "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
