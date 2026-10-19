"""Print the markdown tables for results/gate_thresholds.md from the gate run outputs.

Usage: python3 gate_report.py <project_dir>
Reads results/gate_thresholds.json, results/gate_metrics.csv, results/gate_perpop.csv,
results/leiden_metrics.csv. Prints markdown to stdout. The prose in gate_thresholds.md is
written by hand; this script only keeps the numbers honest.
"""
import json, sys
from pathlib import Path
import pandas as pd

ARM_NAMES = {"A1": "A1 arcsinh", "A2": "A2 z-score", "A3": "A3 percentile", "C": "C lineage mixture + percentile",
             "C2": "C2 lineage mixture + z-score", "D": "D all mixture",
             "C_gate": "C_gate lineage gate + z-score", "D_gate": "D_gate all gate", "D_gate_clean": "D_gate_clean 15 clean markers"}
ORDER = ["A1", "A2", "A3", "C2", "C_gate", "D", "D_gate", "D_gate_clean", "C"]


def f(v, d=2):
    return "" if pd.isna(v) else f"{v:.{d}f}"


def main(root):
    root = Path(root); res = root / "results"
    g = json.load(open(res / "gate_thresholds.json"))["markers"]
    print("## Thresholds per marker\n")
    print("| Marker | Method | Gate cut | Mixture cut | Balanced accuracy, gate | Balanced accuracy, mixture | Positive fraction, gate | Positive fraction, mixture | Positive populations | Negative populations | Note |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for m, r in g.items():
        pos = ", ".join(p.replace("_", " ") for p in r["positive_populations"]) or "none"
        neg = ", ".join(p.replace("_", " ") for p in r["negative_populations"]) or "none"
        ba = f(r.get("balanced_accuracy")); bam = f(r.get("mixture_balanced_accuracy"))
        cut = f(r["threshold"]) if r["method"] == "gate" else ""
        print(f"| {m} | {r['method'].replace('_', ' ')} | {cut} | {f(r['mixture_threshold'])} | {ba} | {bam} | {f(r['fraction_positive_all_cells'])} | {f(r['fraction_positive_mixture'])} | {pos} | {neg} | {r['note']} |")
    mt = pd.read_csv(res / "gate_metrics.csv")
    for name, k in [("Levine_32dim", 14), ("Levine_32dim_k40", 40)]:
        d = mt[(mt.dataset == name) & mt.arm.isin(ORDER)].set_index("arm").loc[[a for a in ORDER if a in set(mt[mt.dataset == name].arm)]]
        print(f"\n## FlowSOM, {k} metaclusters, full data, seed 1\n")
        print("| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | Occupied phenotypes | Phenotypes for 95 percent of cells |")
        print("|---|---|---|---|---|---|---|")
        for a, r in d.iterrows():
            occ = "" if pd.isna(r.occupied_phenotypes) else f"{int(r.occupied_phenotypes):,}"
            c95 = "" if pd.isna(r.phenotypes_covering_95pct) else f"{int(r.phenotypes_covering_95pct):,}"
            print(f"| {ARM_NAMES[a]} | {f(r.mean_f1)} | {f(r.ari)} | {f(r.nmi)} | {f(r.stability_ari)} | {occ} | {c95} |")
    pp = pd.read_csv(res / "gate_perpop.csv")
    for name, k in [("Levine_32dim_k40", 40), ("Levine_32dim", 14)]:
        d = pp[pp.dataset == name].pivot(index="population", columns="arm", values="f1")
        n = pp[pp.dataset == name].groupby("population").n.first()
        d = d.join(n).sort_values("n", ascending=False)
        cols = [a for a in ["A1", "C2", "C_gate", "D", "D_gate", "D_gate_clean"] if a in d.columns]
        print(f"\n## Per-population F1, FlowSOM, {k} metaclusters\n")
        print("| Population | Cells | " + " | ".join(ARM_NAMES[a].split(" ")[0] for a in cols) + " |")
        print("|---|---|" + "---|" * len(cols))
        for p, r in d.iterrows():
            print(f"| {p.replace('_', ' ')} | {int(r.n):,} | " + " | ".join(f(r[a]) for a in cols) + " |")
    lf = res / "leiden_metrics.csv"
    if lf.exists():
        lm = pd.read_csv(lf)
        print("\n## Leiden on a kNN graph, 50,000-cell subsample, 15 neighbors\n")
        print("| Arm | Distance | Features | Resolution | Clusters | Mean F1 | ARI | NMI | Stability ARI, 10 seeds |")
        print("|---|---|---|---|---|---|---|---|---|")
        for _, r in lm.iterrows():
            print(f"| {r.arm} | {r.metric} | {r.n_features} | {r.resolution} | {r.n_clusters} | {f(r.mean_f1)} | {f(r.ari)} | {f(r.nmi)} | {f(r.stability_ari)} |")


if __name__ == "__main__":
    main(sys.argv[1])
