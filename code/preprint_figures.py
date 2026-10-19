"""Figures for the preprint. Reads results/*.csv and the two Stanford HIPC files, writes
manuscript/figures/fig*.png and .pdf. Usage: python3 preprint_figures.py <project_dir>"""
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#0b0b0b", "#52514e", "#dddcd7"
COL = {"A1": "#2a78d6", "A2": "#eb6834", "A3": "#1baf7a", "B_geom": "#eda100", "B_lin": "#e87ba4", "C": "#008300", "C2": "#e34948", "D": "#6250d6", "manual": "#777777"}
LBL = {"A1": "arcsinh (A1)", "A2": "z-score (A2)", "A3": "percentile (A3)", "B_geom": "brightness 2^rank (B)", "B_lin": "brightness linear (B)",
       "C": "lineage bits + percentile (C)", "C2": "lineage bits + z-score (C2)", "D": "all bits (D)", "manual": "central manual gates"}
DS = [("Levine_32dim_k40", "Levine 32"), ("Levine_13dim_k40", "Levine 13"), ("Samusik_01_k40", "Samusik 01"), ("Samusik_all_k40", "Samusik all"), ("Nilsson_rare", "Nilsson rare"), ("Mosmann_rare", "Mosmann rare")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 200, "savefig.dpi": 300})


def load_metrics(res):
    ms = []
    for tag in ["pilot", "phase2", "mosmannB", "samusik_all"]:
        m = pd.read_csv(res / f"{tag}_metrics.csv"); m["tag"] = tag; ms.append(m)
    m = pd.concat(ms, ignore_index=True); single = m[~m.arm.str.contains("~")].drop_duplicates(["dataset", "arm"], keep="last")
    single["score"] = single.mean_f1.where(single.mean_f1.notna(), single.rare_f1)
    return single, m[m.arm.str.contains("~")].drop_duplicates(["dataset", "arm"], keep="last")


def bars(ax, cats, arms, val, ylabel, labels=True):
    n = len(arms); w = 0.8 / n
    for i, a in enumerate(arms):
        x = np.arange(len(cats)) + (i - (n - 1) / 2) * w
        v = [val.get((c, a), np.nan) for c in cats]
        ax.bar(x, v, w * 0.9, color=COL[a], label=LBL[a], linewidth=0)
        if labels:
            for xi, vi in zip(x, v):
                if not np.isnan(vi): ax.text(xi, vi + 0.012, f"{vi:.2f}", ha="center", va="bottom", fontsize=5.5, color=INK2)
    ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats, rotation=20 if len(cats) > 4 else 0, ha="right" if len(cats) > 4 else "center"); ax.set_ylim(0, 1.05); ax.set_ylabel(ylabel)
    ax.grid(axis="y", color=GRID, linewidth=0.6); ax.set_axisbelow(True)


def save(fig, res, name):
    fig.savefig(res / f"{name}.png", bbox_inches="tight"); fig.savefig(res / f"{name}.pdf", bbox_inches="tight"); plt.close(fig)


def fig1(single, out):
    """Scaling cost: F1 and ARI for A1, A2, A3 per dataset."""
    cats = [l for _, l in DS]; names = dict(DS)
    f1 = {(names[r.dataset], r.arm): r.score for _, r in single.iterrows() if r.dataset in names}
    ari = {(names[r.dataset], r.arm): r.ari for _, r in single.iterrows() if r.dataset in names}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6))
    bars(axes[0], cats, ["A1", "A2", "A3"], f1, "F1"); axes[0].set_title("a  Recovery of the gated populations", loc="left", fontsize=8, color=INK)
    bars(axes[1], cats[:4], ["A1", "A2", "A3"], ari, "ARI"); axes[1].set_title("b  Adjusted Rand index, labeled sets", loc="left", fontsize=8, color=INK)
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower left")
    fig.tight_layout(); save(fig, out, "fig1_scaling_cost")


def fig2(res, out):
    """Where the loss falls: per-population F1 heatmap across arms."""
    pp = pd.concat([pd.read_csv(res / f"{t}_perpop.csv").assign(tag=t) for t in ["pilot", "phase2", "samusik_all"]]).drop_duplicates(["dataset", "arm", "population"], keep="last")
    want = [("Levine_32dim_k40", "CD4_T_cells"), ("Levine_32dim_k40", "CD34+CD38lo_HSCs"), ("Levine_32dim_k40", "Pro_B_cells"), ("Levine_32dim_k40", "Plasma_B_cells"), ("Levine_32dim_k40", "CD34+CD38+CD123+_HSPCs"),
            ("Levine_13dim_k40", "CD11b-_Monocyte_cells"), ("Levine_13dim_k40", "CD11bmid_Monocyte_cells"), ("Levine_13dim_k40", "CD11bhi_Monocyte_cells"), ("Levine_13dim_k40", "Mature_CD38mid_B_cells"),
            ("Samusik_all_k40", "Classical Monocytes"), ("Samusik_all_k40", "Intermediate Monocytes"), ("Samusik_all_k40", "Non-Classical Monocytes")]
    arms = ["A1", "A2", "A3", "C2", "D"]; names = dict(DS)
    rows = []; labels = []
    for ds, pop in want:
        sub = pp[(pp.dataset == ds) & (pp.population == pop)]
        if not len(sub): continue
        rows.append([float(sub[sub.arm == a].f1.iloc[0]) if (sub.arm == a).any() else np.nan for a in arms]); labels.append(f"{pop.replace('_', ' ')}  ·  {names[ds]}  ·  n {int(sub.n.iloc[0]):,}")
    # the rare sets as two more rows: best-cluster F1
    single, _ = load_metrics(res)
    for ds, lab in [("Mosmann_rare", "activated memory CD4 T cells  ·  Mosmann rare  ·  n 109"), ("Nilsson_rare", "HSC  ·  Nilsson rare  ·  n 358")]:
        s = single[single.dataset == ds].set_index("arm"); rows.append([float(s.loc[a].score) if a in s.index else np.nan for a in arms]); labels.append(lab)
    M = np.array(rows)
    fig, ax = plt.subplots(figsize=(5.2, 0.28 * len(rows) + 0.9))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#f4f3ef", "#cde2fb", "#6da7ec", "#2a78d6", "#0d366b"])
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(arms))); ax.set_xticklabels([LBL[a] for a in arms], rotation=25, ha="right", fontsize=7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if not np.isnan(M[i, j]): ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6, color="white" if M[i, j] > 0.6 else INK)
    ax.set_title("Per-population F1 at the free cluster number", loc="left", fontsize=8, color=INK); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    fig.tight_layout(); save(fig, out, "fig2_where_the_loss_falls")


def fig3(single, pairs, out):
    """The boundary: every arm per dataset, and arm B agreement."""
    cats = [l for _, l in DS]; names = dict(DS)
    f1 = {(names[r.dataset], r.arm): r.score for _, r in single.iterrows() if r.dataset in names}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(width_ratios=[2.2, 1]))
    bars(axes[0], cats, ["A1", "A2", "A3", "B_geom", "B_lin", "C2", "D"], f1, "F1", labels=False)
    axes[0].set_title("a  Recovery under every representation", loc="left", fontsize=8, color=INK); axes[0].legend(frameon=False, fontsize=6, ncol=2, loc="upper right")
    ax = axes[1]; want = ["A1~A2", "A1~B_geom", "A1~B_lin", "A2~B_geom", "A2~B_lin"]
    for i, ds in enumerate(["Nilsson_rare", "Mosmann_rare"]):
        sub = pairs[pairs.dataset == ds].set_index("arm").arm_agreement_ari
        ax.barh(np.arange(len(want)) + (i - 0.5) * 0.38, [sub.get(w, np.nan) for w in want], 0.36, color=["#2a78d6", "#5b8fd9"][i], label=names[ds], linewidth=0)
    ax.set_yticks(range(len(want))); ax.set_yticklabels([w.replace("~", " vs ").replace("B_geom", "B 2^rank").replace("B_lin", "B linear") for w in want], fontsize=6.5)
    ax.set_xlim(0, 1); ax.set_xlabel("ARI between clusterings"); ax.invert_yaxis(); ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    ax.set_title("b  Brightness division against scaling", loc="left", fontsize=8, color=INK); ax.grid(axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    fig.tight_layout(); save(fig, out, "fig3_representations")


def fig4(res, out):
    """HIPC between-site spread per arm at both granularities, and the manual reference."""
    b12 = pd.read_csv(res / "hipc_between_site_k12.csv"); b30 = pd.read_csv(res / "hipc_between_site.csv")
    order = ["manual", "A1", "A2", "A3", "B_geom", "B_lin", "C2", "D"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axes[0]; x = np.arange(len(order))
    v12 = [float(b12[b12.arm == a].median_between_site_cv_t.iloc[0]) for a in order]; v30 = [float(b30[b30.arm == a].median_between_site_cv_t.iloc[0]) if (b30.arm == a).any() else np.nan for a in order]
    ax.bar(x - 0.19, v12, 0.36, color=[COL[a] for a in order], linewidth=0, label="12 metaclusters")
    ax.bar(x + 0.19, v30, 0.36, color=[COL[a] for a in order], alpha=0.45, linewidth=0, label="30 metaclusters")
    within = float(b12[b12.arm == "manual"].median_within_site_cv_t.iloc[0]); ax.axhline(within, color=INK2, linestyle=":", linewidth=1, label=f"replicate noise within a center, {within:.2f}")
    ax.set_xticks(x); ax.set_xticklabels([LBL[a] for a in order], fontsize=6, rotation=25, ha="right"); ax.set_ylabel("median between-site CV, among T cells")
    ax.set_ylim(0, 0.68); ax.legend(frameon=False, fontsize=6.5, loc="upper left"); ax.grid(axis="y", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    ax.set_title("a  Seven centers, three lots, per representation", loc="left", fontsize=8, color=INK)
    man = pd.read_csv(res / "hipc_between_site_manual_k12.csv"); man = man[man.denominator == "T cells"]
    piv = man.groupby("cluster").between_site_cv.median().sort_values(); ax = axes[1]
    ax.barh(range(len(piv)), piv.values, color=COL["manual"], linewidth=0); ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=6)
    ax.set_xlabel("between-site CV, of CD3"); ax.set_title("b  The central manual gates themselves", loc="left", fontsize=8, color=INK); ax.grid(axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    fig.tight_layout(); save(fig, out, "fig4_between_site")


def fig5(root, out):
    """HIPC deposited transform: README logicle fails, FlowJo logicle with per-channel width fits."""
    import flowio
    from flowutils import transforms as FT
    sys.path.insert(0, str(root / "code")); from hipc_transform_check import read_fj_matrix, read, DETECTOR_TO_MARKER
    hdr, Mws, _ = read_fj_matrix(root / "data/HIPC/XML/Stanford/Stanford Matrices/Tcell")
    name = "1228-1_A1_A01.fcs"; pr, R, _ = read(root / "data/HIPC/RAW/Stanford" / name); pdd, D, _ = read(root / "data/HIPC/FCS/Stanford" / name)
    X = R[:, [pr.index(d) for d in hdr]] @ np.linalg.inv(Mws)
    W = json.load(open(root / "results/hipc_meta.json"))["HIPC_tcell"] if (root / "results/hipc_meta.json").exists() else None
    widths = {"Live": 0.40, "CCR7": 0.90, "CD4": 0.80, "CD45RA": 0.80, "CD3": 0.15, "HLADR": 0.75, "CD38": 0.60, "CD8": 1.10}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    j = hdr.index("Pacific Blue-A"); k = pdd.index("CD3"); x = X[:, j]; y = D[:, k]; ok = y > -4096
    sub = np.random.default_rng(0).choice(np.where(ok)[0], 20000, replace=False)
    y_readme = FT.logicle(x[sub].reshape(-1, 1), [0], t=16777215, m=7.22471987004958, w=1.16437481092794, a=0).ravel() * 4096
    y_fit = FT.logicle(x[sub].reshape(-1, 1), [0], t=262144, m=4.5, w=widths["CD3"], a=0).ravel() * 4096
    for ax, yy, title in [(axes[0], y_readme, "a  The README iplogicle"), (axes[1], y_fit, "b  Logicle, T 262144, M 4.5, W 0.15")]:
        ax.scatter(y[sub], yy, s=0.5, color="#2a78d6", alpha=0.25, linewidths=0); ax.plot([-1000, 4100], [-1000, 4100], color=INK2, linewidth=0.8, linestyle=":")
        ax.set_xlabel("deposited CD3, 4096 scale"); ax.set_ylabel("reconstructed from raw"); ax.set_title(title, loc="left", fontsize=7.5, color=INK); ax.set_xlim(-1000, 4200); ax.set_ylim(-1000, 4200)
    ax = axes[2]; ms = list(widths); ax.barh(range(len(ms)), [widths[m] for m in ms], color="#2a78d6", linewidth=0); ax.set_yticks(range(len(ms))); ax.set_yticklabels(ms, fontsize=7)
    ax.set_xlabel("fitted width W, decades"); ax.set_title("c  Width per channel, Stanford", loc="left", fontsize=7.5, color=INK); ax.grid(axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    fig.tight_layout(); save(fig, out, "fig5_hipc_transform")


def fig6(res, out):
    """Cofactor sweep: F1 under A1, A2, A3 per dataset across cofactors."""
    f = res / "cofactor_metrics.csv"
    if not f.exists(): return
    m = pd.read_csv(f); m = m[~m.arm.str.contains("~")]; m["score"] = m.mean_f1.where(m.mean_f1.notna(), m.rare_f1)
    single, _ = load_metrics(res)
    base = {("Levine_32dim_k40", 5), ("Samusik_01_k40", 5), ("Levine_13dim_k40", 5), ("Nilsson_rare", 150), ("Mosmann_rare", 150)}
    rows = [dict(base_name=r.dataset.split("@")[0], cf=float(r.dataset.split("@cf=")[1]), arm=r.arm, score=r.score) for _, r in m.iterrows()]
    for ds, cf in base:
        for a in ["A1", "A2", "A3"]:
            s = single[(single.dataset == ds) & (single.arm == a)]
            if len(s): rows.append(dict(base_name=ds, cf=cf, arm=a, score=float(s.score.iloc[0])))
    d = pd.DataFrame(rows); names = dict(DS)
    sets = [x for x in ["Levine_32dim_k40", "Levine_13dim_k40", "Samusik_01_k40", "Nilsson_rare", "Mosmann_rare"] if x in set(d.base_name)]
    fig, axes = plt.subplots(1, len(sets), figsize=(1.5 * len(sets), 2.3), sharey=True)
    for ax, ds in zip(np.atleast_1d(axes), sets):
        sub = d[d.base_name == ds]
        for a in ["A1", "A2", "A3"]:
            s = sub[sub.arm == a].sort_values("cf"); ax.plot(s.cf, s.score, marker="o", markersize=3, color=COL[a], label=LBL[a], linewidth=1.2)
        ax.set_xscale("log"); ax.set_title(names[ds], fontsize=7.5, color=INK); ax.set_xlabel("cofactor"); ax.grid(color=GRID, linewidth=0.5); ax.set_ylim(0, 1)
    np.atleast_1d(axes)[0].set_ylabel("F1"); np.atleast_1d(axes)[0].legend(frameon=False, fontsize=6, loc="lower left")
    fig.tight_layout(); save(fig, out, "fig6_cofactor")


def main(root):
    root = Path(root); res = root / "results"; out = root / "manuscript/figures"; out.mkdir(parents=True, exist_ok=True)
    single, pairs = load_metrics(res)
    fig1(single, out); fig2(res, out); fig3(single, pairs, out); fig4(res, out); fig5(root, out); fig6(res, out)
    print("figures written to", out)


if __name__ == "__main__":
    main(sys.argv[1])
