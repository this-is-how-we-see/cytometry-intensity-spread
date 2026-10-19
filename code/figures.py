"""Figures, one per hypothesis. Reads results/{pilot,phase2,mosmannB}_*.csv, writes results/fig*.png.

Palette: the dataviz reference categorical palette, fixed slot per arm (never cycled), light mode.
Usage: python3 figures.py <project_dir>
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
ARM_COLOR = {"A1": "#2a78d6", "A2": "#eb6834", "A3": "#1baf7a", "B_geom": "#eda100",
             "B_lin": "#e87ba4", "C": "#008300", "C2": "#e34948", "D": "#6250d6"}
ARM_LABEL = {"A1": "A1 arcsinh", "A2": "A2 z-score", "A3": "A3 percentile", "B_geom": "B brightness (2^rank)",
             "B_lin": "B brightness (linear rank)", "C": "C lineage binary + percentile", "C2": "C2 lineage binary + z-score", "D": "D all binary"}
ARMS = list(ARM_COLOR)
TAGS = ["pilot", "phase2", "mosmannB", "samusik_all", "hipc"]
DATASETS = ["Nilsson_rare", "Mosmann_rare", "Levine_32dim", "Levine_32dim_k40", "Levine_13dim", "Levine_13dim_k40",
            "Samusik_01", "Samusik_01_k40", "Samusik_all", "Samusik_all_k40", "HIPC_tcell", "HIPC_tcell_k30"]
SEQ = LinearSegmentedColormap.from_list("blue", ["#f0efec", "#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})

def style(ax, ylabel=None, xlabel=None):
    ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
    if ylabel: ax.set_ylabel(ylabel)
    if xlabel: ax.set_xlabel(xlabel)

def bars(ax, cats, arms, values, width=None, labels=True):
    """Grouped bars: one group per category, one bar per arm, 2px gap, fixed arm colors."""
    n = len(arms); width = width or 0.8 / n
    for i, arm in enumerate(arms):
        x = np.arange(len(cats)) + (i - (n - 1) / 2) * width
        v = [values.get((c, arm), np.nan) for c in cats]
        ax.bar(x, v, width * 0.92, color=ARM_COLOR[arm], label=ARM_LABEL[arm], linewidth=0)
        for xi, vi in zip(x, v):
            if labels and not np.isnan(vi): ax.text(xi, vi + 0.01, f"{vi:.2f}", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_xticks(np.arange(len(cats))); ax.set_xticklabels(cats, rotation=20 if len(cats) > 4 else 0, ha="right" if len(cats) > 4 else "center")

def main(root):
    root = Path(root); res = root / "results"
    # Merge the pilot run with the later runs that share its arms. A later file wins on a
    # (dataset, arm) pair, so the mosmannB run replaces the pilot's Mosmann rows (A1 to A3 are
    # identical, seed 1) and adds its B arms.
    ms, ps, meta = [], [], {}
    for tag in TAGS:
        if (res / f"{tag}_metrics.csv").exists():
            ms.append(pd.read_csv(res / f"{tag}_metrics.csv"))
            try:
                ps.append(pd.read_csv(res / f"{tag}_perpop.csv"))
            except pd.errors.EmptyDataError:   # a rare-cell run writes an empty per-population file
                pass
            meta.update(json.load(open(res / f"{tag}_meta.json")))
    m = pd.concat(ms, ignore_index=True).drop_duplicates(["dataset", "arm"], keep="last")
    per = pd.concat(ps, ignore_index=True).drop_duplicates(["dataset", "arm", "population"], keep="last")
    datasets = [d for d in DATASETS if d in set(m.dataset)]
    single = m[~m.arm.str.contains("~")]; pairs = m[m.arm.str.contains("~")]

    # ---- Figure 1: hypothesis 1, the brightness term
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0), gridspec_kw=dict(width_ratios=[1, 1.4]))
    ax = axes[0]
    div = meta.get("Nilsson_rare", {}).get("divisors", {}).get("B_geom", {})
    fl = meta.get("Nilsson_rare", {}).get("fluor", {})
    if div:
        names = sorted(div, key=lambda k: div[k])
        ax.barh(names, [div[k] for k in names], color="#2a78d6", height=0.7, linewidth=0)
        for i, k in enumerate(names): ax.text(div[k] + 0.03, i, fl.get(k, ""), va="center", fontsize=8, color=INK2)
        ax.axvline(1, color=INK2, linewidth=1, linestyle=":")
        ax.set_title("Nilsson_rare: per-marker brightness divisor (arm B, 2^rank)", fontsize=10, loc="left", color=INK)
        ax.grid(axis="x", color=GRID); ax.set_axisbelow(True); ax.set_xlabel("divisor, geometric mean 1")
    ax = axes[1]
    want = ["A1~A2", "A1~A3", "A1~B_geom", "A2~B_geom", "A1~C", "A1~D", "C~D"]
    rows = []
    for d in datasets:
        sub = pairs[pairs.dataset == d].set_index("arm")["arm_agreement_ari"]
        rows.append([sub.get(w, np.nan) for w in want])
    M = np.array(rows, dtype=float)
    im = ax.imshow(M, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(want))); ax.set_xticklabels([w.replace("~", " vs ") for w in want], rotation=30, ha="right")
    ax.set_yticks(range(len(datasets))); ax.set_yticklabels(datasets)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if not np.isnan(M[i, j]): ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8, color=INK if M[i, j] < 0.6 else SURFACE)
    ax.set_title("Agreement between arms on the full clustering (ARI)", fontsize=10, loc="left", color=INK)
    ax.grid(False)
    fig.tight_layout(); fig.savefig(res / "fig1_brightness.png", dpi=160); plt.close(fig)

    # ---- Figure 2: hypothesis 2, recovery and stability
    fig, axes = plt.subplots(1, 2, figsize=(18, 4.6))
    vals = {}
    for _, r in single.iterrows():
        vals[(r.dataset, r.arm)] = r.mean_f1 if not np.isnan(r.get("mean_f1", np.nan)) else r.rare_f1
    arms_present = [a for a in ARMS if a in set(single.arm)]
    bars(axes[0], datasets, arms_present, vals); style(axes[0], "F1")
    axes[0].set_ylim(0, 1.05); axes[0].set_title("Recovery of gated populations\n(mean F1, Hungarian, for Levine; best-cluster F1 for the rare sets)", fontsize=9, loc="left", color=INK)
    stab = {(r.dataset, r.arm): r.stability_ari for _, r in single.iterrows()}
    bars(axes[1], datasets, arms_present, stab); style(axes[1], "ARI")
    axes[1].set_ylim(0, 1.05); axes[1].set_title("Run-to-run stability\n(mean pairwise ARI over 10 seeds on 50,000 cells)", fontsize=9, loc="left", color=INK)
    axes[1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(); fig.savefig(res / "fig2_f1.png", dpi=160); plt.close(fig)

    # ---- Figure 3: hypothesis 3, per-population F1 on Levine_32dim
    if len(per):
        p = per[per.dataset == ("Levine_32dim_k40" if "Levine_32dim_k40" in set(per.dataset) else "Levine_32dim")]
        arms_l = [a for a in ARMS if a in set(p.arm)]
        pops = p.groupby("population")["n"].first().sort_values(ascending=False)
        M = np.array([[p[(p.population == pop) & (p.arm == a)].f1.values[0] if len(p[(p.population == pop) & (p.arm == a)]) else np.nan for a in arms_l] for pop in pops.index])
        fig, ax = plt.subplots(figsize=(9, 0.4 * len(pops) + 1.8))
        ax.imshow(M, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(arms_l))); ax.set_xticklabels([ARM_LABEL[a] for a in arms_l], rotation=20, ha="right")
        ax.set_yticks(range(len(pops))); ax.set_yticklabels([f"{n}  (n={c:,})" for n, c in pops.items()])
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8, color=INK if M[i, j] < 0.6 else SURFACE)
        ax.set_title(f"{p.dataset.iloc[0]}: per-population F1 by arm\n(populations sorted by size)", fontsize=10, loc="left", color=INK)
        ax.grid(False); fig.tight_layout(); fig.savefig(res / "fig3_failure.png", dpi=160); plt.close(fig)

    # ---- Figure 4: dimensionality
    fig, axes = plt.subplots(1, 2, figsize=(18, 4.6))
    ax = axes[0]
    dsets4 = [d for d in datasets if not d.endswith("_k40")]   # same cells and same bits as the k=14 run
    for i, d in enumerate(dsets4):
        for j, arm in enumerate(["C", "D"]):  # C2 has the same binary columns as C
            r = single[(single.dataset == d) & (single.arm == arm)]
            if not len(r) or np.isnan(r.occupied_phenotypes.values[0]): continue
            occ, pos, c95 = r.occupied_phenotypes.values[0], r.possible_phenotypes.values[0], r.phenotypes_covering_95pct.values[0]
            x = i + (j - 0.5) * 0.36
            ax.bar(x, pos, 0.3, color=GRID, linewidth=0)
            ax.bar(x, occ, 0.3, color=ARM_COLOR[arm], linewidth=0)
            ax.bar(x, c95, 0.3, color=INK, linewidth=0, alpha=0.85)
            ax.text(x, pos * 1.3, f"{int(occ):,}\nof {int(pos):,}\n95% in {int(c95):,}", ha="center", fontsize=7, color=INK2)
    ax.set_yscale("log"); ax.set_ylim(top=10 ** 12); ax.set_xticks(range(len(dsets4))); ax.set_xticklabels(dsets4)
    ax.set_xlim(-0.6, len(dsets4) - 0.4)
    ax.set_title("Boolean phenotypes\n(grey: possible 2^k, color: occupied, black: enough to cover 95% of cells)", fontsize=9, loc="left", color=INK)
    style(ax, "count, log scale")
    ax = axes[1]
    idv = {(r.dataset, r.arm): r.two_nn_id for _, r in single.iterrows()}
    bars(ax, datasets, arms_present, idv, labels=False); style(ax, "two-NN intrinsic dimension")
    ax.set_title("Intrinsic dimension per arm", fontsize=10, loc="left", color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(); fig.savefig(res / "fig4_dimensionality.png", dpi=160); plt.close(fig)
    print("figures written")

if __name__ == "__main__":
    main(sys.argv[1])
