"""Print the markdown tables for results/phase2.md (or any later run, with --tag) from the run outputs.

Usage: python3 phase2_report.py <project_dir> [--tag phase2]
Reads results/<tag>_metrics.csv and results/<tag>_perpop.csv. Prints markdown to stdout.
The prose in phase2.md is written by hand; this script only keeps the numbers honest.
"""
import sys
from pathlib import Path
import pandas as pd

ARM_NAMES = {"A1": "A1 arcsinh", "A2": "A2 z-score", "A3": "A3 percentile",
             "C": "C lineage mixture + percentile", "C2": "C2 lineage mixture + z-score", "D": "D all mixture",
             "W": "W resolution-weighted arcsinh", "W_z": "W_z resolution-weighted z-score", "W2": "W2 sharpness-weighted arcsinh", "W2s": "W2s sharpness squared", "P": "P posterior + z-score",
             "B_geom": "B brightness, 2^rank", "B_lin": "B brightness, linear rank", "B_meas": "B measured divisor"}
ORDER = ["A1", "A2", "A3", "W", "W_z", "W2", "W2s", "B_geom", "B_lin", "B_meas", "P", "C2", "D", "C"]


def jobs_in(mt):
    """(dataset, k) pairs in the order the metrics file lists them."""
    seen = []
    for _, r in mt.iterrows():
        if (r.dataset, int(r.k)) not in seen: seen.append((r.dataset, int(r.k)))
    return seen


def f(v, d=2):
    return "" if pd.isna(v) else f"{v:.{d}f}"


def main(root, tag="phase2"):
    res = Path(root) / "results"
    mt = pd.read_csv(res / f"{tag}_metrics.csv")
    pp = pd.read_csv(res / f"{tag}_perpop.csv")
    JOBS = jobs_in(mt)
    for name, k in JOBS:
        sub = mt[mt.dataset == name]
        if sub.empty: continue
        d = sub[sub.arm.isin(ORDER)].set_index("arm").loc[[a for a in ORDER if a in set(sub.arm)]]
        agree = {r.arm.split("~")[1]: r.arm_agreement_ari for _, r in sub[sub.arm.str.startswith("A1~")].iterrows()}
        print(f"\n## {name.replace('_k40', '')}, FlowSOM, {k} metaclusters, full data, seed 1\n")
        print("| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |")
        print("|---|---|---|---|---|---|---|---|---|")
        for a, r in d.iterrows():
            occ = "" if pd.isna(r.get("occupied_phenotypes")) else f"{int(r.occupied_phenotypes):,}"
            pos = "" if pd.isna(r.get("possible_phenotypes")) else f"{int(r.possible_phenotypes):,}"
            c95 = "" if pd.isna(r.get("phenotypes_covering_95pct")) else f"{int(r.phenotypes_covering_95pct):,}"
            ag = "" if a == "A1" else f(agree.get(a))
            print(f"| {ARM_NAMES[a]} | {f(r.mean_f1)} | {f(r.ari)} | {f(r.nmi)} | {f(r.stability_ari)} | {ag} | {occ} | {pos} | {c95} |")
    for name, k in JOBS:
        sub = pp[pp.dataset == name]
        if sub.empty: continue
        d = sub.pivot(index="population", columns="arm", values="f1")
        n = sub.groupby("population").n.first()
        d = d.join(n).sort_values("n", ascending=False)
        cols = [a for a in ORDER if a in d.columns]
        print(f"\n## Per-population F1, {name.replace('_k40', '')}, {k} metaclusters\n")
        print("| Population | Cells | " + " | ".join(a for a in cols) + " |")
        print("|---|---|" + "---|" * len(cols))
        for p, r in d.iterrows():
            print(f"| {p.replace('_', ' ')} | {int(r.n):,} | " + " | ".join(f(r[a]) for a in cols) + " |")


if __name__ == "__main__":
    tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else "phase2"
    main(sys.argv[1], tag)
