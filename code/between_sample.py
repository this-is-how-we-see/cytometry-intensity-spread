"""Between-mouse variance of cluster frequencies per arm, on the Samusik_all clusterings.

The plan reserves between-site variance of cluster frequencies for the HIPC data. Samusik_all
holds the same panel on 10 mice, so this is the same measurement with mice in place of sites.

For each arm, and for the manual labels as the reference:
  1. the frequency table, cluster by mouse, as a fraction of that mouse's cells
  2. the coefficient of variation (SD / mean, 10 mice) per cluster, and its median over the
     clusters whose mean frequency is at least 1 percent, so a 20-cell cluster does not set it
  3. per mouse: mean F1 (Hungarian), ARI and NMI against the labels, then the mean over mice
  4. per population: the Pearson correlation across the 10 mice between the frequency of the
     population and the frequency of its Hungarian-matched cluster. A high value says the
     cluster rises and falls with the population from mouse to mouse.

Frequencies for clusters and for labels are both taken among the labeled cells, so the two
are comparable; the all-cell cluster frequency is reported beside it.

Usage: python3 between_sample.py <project_dir> [--labels results/labels_samusik_all_Samusik_all_k40.npz]
Writes results/samusik_all_between.csv (one row per arm), results/samusik_all_between_clusters.csv
(one row per arm and cluster) and results/samusik_all_between_pops.csv (one row per arm and
population). Prints the markdown tables.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ARM_NAMES = {"A1": "A1 arcsinh", "A2": "A2 z-score", "A3": "A3 percentile",
             "C2": "C2 lineage mixture + z-score", "D": "D all mixture", "truth": "Manual labels"}
ORDER = ["truth", "A1", "A2", "A3", "C2", "D"]


def f1_matrix(true, pred, pops, clusters):
    M = np.zeros((len(pops), len(clusters)))
    for i, p in enumerate(pops):
        t = true == p
        for j, q in enumerate(clusters):
            c = pred == q; tp = np.sum(t & c)
            if tp == 0: continue
            prec, rec = tp / c.sum(), tp / t.sum(); M[i, j] = 2 * prec * rec / (prec + rec)
    return M


def hungarian(true, pred):
    pops = sorted(set(true)); clusters = sorted(set(pred))
    M = f1_matrix(true, pred, pops, clusters)
    r, c = linear_sum_assignment(-M)
    match = {pops[i]: clusters[j] for i, j in zip(r, c)}
    per = {pops[i]: M[i, j] for i, j in zip(r, c)}
    for p in pops: per.setdefault(p, 0.0)
    return float(np.mean(list(per.values()))), per, match


def freq_table(assign, sample, mice):
    """Rows: cluster ids. Columns: mice. Value: fraction of that mouse's cells in the cluster."""
    ids = np.unique(assign)
    T = np.zeros((len(ids), len(mice)))
    for j, m in enumerate(mice):
        sel = sample == m
        counts = pd.Series(assign[sel]).value_counts()
        T[:, j] = [counts.get(i, 0) / sel.sum() for i in ids]
    return pd.DataFrame(T, index=ids, columns=mice)


def cv_rows(T, min_mean=0.01):
    mean = T.mean(axis=1); sd = T.std(axis=1, ddof=1)
    cv = sd / mean.replace(0, np.nan)
    return pd.DataFrame({"mean_freq": mean, "sd_freq": sd, "cv": cv, "used": mean >= min_mean})


def main(root, labels_file):
    root = Path(root); res = root / "results"
    z = np.load(root / labels_file, allow_pickle=True)
    truth, sample = z["truth"].astype(int), z["sample"].astype(int)
    mice = sorted(set(sample.tolist()))
    labeled = truth > 0
    pop_names = {int(k): v for k, v in json.load(open(root / "data/HDCytoData/population_names_Samusik.json")).items()}
    arms = [a for a in ORDER if a == "truth" or a in z.files]
    summary, clusters_out, pops_out = [], [], []
    for arm in arms:
        pred = truth if arm == "truth" else z[arm].astype(int)
        # frequencies among labeled cells, and among all cells
        Tl = freq_table(pred[labeled], sample[labeled], mice)
        cv_l = cv_rows(Tl)
        Ta = freq_table(pred, sample, mice) if arm != "truth" else None
        cv_a = cv_rows(Ta) if Ta is not None else None
        row = dict(arm=arm, n_clusters=int(len(Tl)),
                   median_cv_labeled=float(cv_l.cv[cv_l.used].median()),
                   n_clusters_used=int(cv_l.used.sum()),
                   median_cv_all_cells=float(cv_a.cv[cv_a.used].median()) if cv_a is not None else np.nan,
                   frac_clusters_cv_over_0p5=float((cv_l.cv[cv_l.used] > 0.5).mean()))
        for cid, r in cv_l.iterrows():
            clusters_out.append(dict(arm=arm, cluster=int(cid), mean_freq=r.mean_freq, sd_freq=r.sd_freq, cv=r.cv, used=bool(r.used),
                                     name=pop_names.get(int(cid), "") if arm == "truth" else ""))
        if arm != "truth":
            # per-mouse accuracy
            f1s, aris, nmis = [], [], []
            for m in mice:
                sel = labeled & (sample == m)
                f1, _, _ = hungarian(truth[sel], pred[sel])
                f1s.append(f1); aris.append(adjusted_rand_score(truth[sel], pred[sel])); nmis.append(normalized_mutual_info_score(truth[sel], pred[sel]))
            row.update(mean_f1_all=hungarian(truth[labeled], pred[labeled])[0], mean_f1_per_mouse=float(np.mean(f1s)),
                       sd_f1_per_mouse=float(np.std(f1s, ddof=1)), ari_per_mouse=float(np.mean(aris)), nmi_per_mouse=float(np.mean(nmis)))
            # does the matched cluster track the population across mice?
            _, per, match = hungarian(truth[labeled], pred[labeled])
            Tt = freq_table(truth[labeled], sample[labeled], mice)
            rs = []
            for p, cid in match.items():
                a, b = Tt.loc[p].values, Tl.loc[cid].values
                r_ = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan
                rs.append(r_)
                pops_out.append(dict(arm=arm, population=pop_names[p], n=int(np.sum(truth == p)), f1=per[p],
                                     cluster=int(cid), r_across_mice=r_, cv_population=float(Tt.loc[p].std(ddof=1) / Tt.loc[p].mean()),
                                     cv_cluster=float(Tl.loc[cid].std(ddof=1) / Tl.loc[cid].mean())))
            row.update(median_r_across_mice=float(np.nanmedian(rs)), frac_r_over_0p8=float(np.nanmean(np.array(rs) > 0.8)))
        summary.append(row)
    S = pd.DataFrame(summary); C = pd.DataFrame(clusters_out); P = pd.DataFrame(pops_out)
    S.to_csv(res / "samusik_all_between.csv", index=False); C.to_csv(res / "samusik_all_between_clusters.csv", index=False)
    P.to_csv(res / "samusik_all_between_pops.csv", index=False)

    def f(v, d=2): return "" if pd.isna(v) else f"{v:.{d}f}"
    print(f"## Between-mouse variance, {labels_file}, 10 mice\n")
    print("| Arm | Clusters | Clusters at 1 percent or more | Median CV, labeled cells | Median CV, all cells | Share of clusters with CV over 0.5 | Mean F1, all mice pooled | Mean F1 per mouse (SD) | ARI per mouse | Median r across mice | Share of populations with r over 0.8 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in S.iterrows():
        print(f"| {ARM_NAMES[r.arm]} | {r.n_clusters} | {r.n_clusters_used} | {f(r.median_cv_labeled)} | {f(r.median_cv_all_cells)} | {f(r.frac_clusters_cv_over_0p5)} | {f(r.get('mean_f1_all'))} | {f(r.get('mean_f1_per_mouse'))} ({f(r.get('sd_f1_per_mouse'))}) | {f(r.get('ari_per_mouse'))} | {f(r.get('median_r_across_mice'))} | {f(r.get('frac_r_over_0p8'))} |")
    print("\n## Per population: does the matched cluster track the population across mice?\n")
    piv_r = P.pivot(index="population", columns="arm", values="r_across_mice")
    piv_f = P.pivot(index="population", columns="arm", values="f1")
    n = P.groupby("population").n.first(); cvp = P.groupby("population").cv_population.first()
    arms_p = [a for a in ORDER if a in piv_r.columns]
    print("| Population | Cells | CV of the population across mice | " + " | ".join(f"r {a}" for a in arms_p) + " | " + " | ".join(f"F1 {a}" for a in arms_p) + " |")
    print("|---|---|---|" + "---|" * (2 * len(arms_p)))
    for p in n.sort_values(ascending=False).index:
        print(f"| {p} | {int(n[p]):,} | {f(cvp[p])} | " + " | ".join(f(piv_r.loc[p, a]) for a in arms_p) + " | " + " | ".join(f(piv_f.loc[p, a]) for a in arms_p) + " |")


if __name__ == "__main__":
    lab = sys.argv[sys.argv.index("--labels") + 1] if "--labels" in sys.argv else "results/labels_samusik_all_Samusik_all_k40.npz"
    main(sys.argv[1], lab)
