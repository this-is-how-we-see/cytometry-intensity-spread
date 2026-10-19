"""Between-site variance of cluster frequencies per arm on the HIPC Lyoplate T cell panel.

Design of the data: three SeraCare lots (12828, 1349, 1369), each stained and acquired in
triplicate at seven centers with the same lyophilized BD reagents and centrally gated. For a
given lot the true cell composition is the same at every center, so the spread of a cluster's
frequency across centers is instrument and handling variance, and lower is better. This is
the property the Samusik_all run lacked.

For each arm, from the saved labels (results/labels_hipc_HIPC_tcell_k30.npz, which carries site,
lot, replicate and file per cell):
  1. frequency of each cluster in each file, relative to that file's events in the subsample
  2. per lot and cluster: mean over the three replicates within a center, then the coefficient
     of variation across the seven centers (between-site CV); and the CV across the three
     replicates within a center, averaged over centers (within-site CV, the floor)
  3. medians over the clusters with mean frequency of at least 1 percent, per lot and pooled
The reference is the central manual gating in data/HIPC/popstats_manual/popstats-tcell.csv.
Two denominators, because the denominator decides what the CV measures:
  * all events: population Count over the file's total events. This inherits the lymphocyte
    gate, and the debris fraction differs by center, so the manual CVs here are about 0.25.
  * T cells: population Count over the file's CD3 count. This is the number Finak 2016
    reports (CV 0.04 to 0.37). For the clusters the T cell denominator is the events in the
    clusters whose median CD3 lies above the CD3 mixture cut of arm A1, so each arm defines
    its own T cell pool from its own clusters.

Matching clusters to manual populations without per-cell labels: for each manual population,
the cluster whose frequency across the 63 files correlates best with it (Pearson), reported
with the correlation. Three lots differ in composition, so the correlation carries signal, but
it is weaker evidence than a Hungarian match on cell labels and the report says so.

Usage: python3 between_site.py <project_dir> [--labels results/labels_hipc_HIPC_tcell_k30.npz] [--prefix hipc_between_site]
Writes results/<prefix>.csv, results/<prefix>_clusters.csv, results/<prefix>_manual.csv,
results/<prefix>_match.csv and prints markdown.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ARM_NAMES = {"A1": "A1 logicle as deposited", "A2": "A2 z-score", "A3": "A3 percentile", "B_geom": "B brightness, 2^rank",
             "B_lin": "B brightness, linear rank", "B_site": "B brightness, stain index per center",
             "B_pool": "B brightness, stain index pooled over centers", "C2": "C2 lineage mixture + z-score", "D": "D all mixture",
             "C": "C lineage mixture + percentile", "manual": "Central manual gates"}
ORDER = ["manual", "A1", "A2", "A3", "B_geom", "B_lin", "B_pool", "B_site", "C2", "D", "C"]
MIN_MEAN = 0.01


def cv(x):
    x = np.asarray(x, float); m = x.mean()
    return float(x.std(ddof=1) / m) if m > 0 and len(x) > 1 else np.nan


def site_stats(F, meta):
    """F: DataFrame index=file, columns=cluster or population, values=frequency.
    meta: DataFrame index=file with site, lot, replicate. Returns per (lot, column) rows."""
    rows = []
    for lot, g in meta.groupby("lot"):
        Fl = F.loc[g.index]
        site_mean = Fl.groupby(g.site).mean()                      # mean of replicates per site
        within = Fl.groupby(g.site).agg(lambda s: cv(s)).mean(0)    # replicate CV, averaged over sites
        for col in F.columns:
            rows.append(dict(lot=lot, cluster=col, mean_freq=float(Fl[col].mean()),
                             between_site_cv=cv(site_mean[col]), within_site_cv=float(within[col]),
                             n_sites=int(site_mean[col].notna().sum())))
    return pd.DataFrame(rows)


def tcell_mask(root, z, arms):
    """Per arm, the events in clusters whose median CD3 (deposited scale) is above the A1
    mixture cut. Reloads the dataset with the same seed, so the subsample matches the labels."""
    sys.path.insert(0, str(Path(root) / "code")); import pilot
    ds = pilot.load_hipc_tcell(Path(root))
    cd3 = ds["X"][:, ds["markers"].index("CD3")]
    assert len(cd3) == len(z["site"]), "subsample does not match the labels file"
    cut = pilot.gmm_threshold(cd3)
    out = {}
    for arm in arms:
        pred = z[arm].astype(int)
        med = pd.Series(cd3).groupby(pred).median()
        out[arm] = np.isin(pred, med.index[med > cut].values)
    return out, float(cut)


def main(root, labels_file, prefix="hipc_between_site"):
    root = Path(root); res = root / "results"
    z = np.load(root / labels_file, allow_pickle=True)
    site, lot, rep, fname = (z[k].astype(str) for k in ("site", "lot", "replicate", "fcs_file"))
    key = pd.Series([f"{s}|{l}|{r}|{f}" for s, l, r, f in zip(site, lot, rep, fname)])
    meta = pd.DataFrame({"key": key, "site": site, "lot": lot, "replicate": rep}).drop_duplicates("key").set_index("key")
    arms = [a for a in ORDER if a in z.files]
    # manual reference: population count / total events per file
    man = pd.read_csv(root / "data/HIPC/popstats_manual/popstats-tcell.csv")
    tot = man[man.Population == "Lymphocytes"].set_index("Filename").Count_Parent          # root count per file
    man["freq"] = man.Count / man.Filename.map(tot)
    cd3n = man[man.Population == "CD3"].set_index("Filename").Count
    man["freq_t"] = man.Count / man.Filename.map(cd3n)          # relative to CD3, the Finak denominator
    man["key"] = man.Center.str.replace("Baylor", "Baylor") + "|" + man.Sample.astype(str) + "|" + man.Replicate.astype(str) + "|" + man.Filename
    # the labels carry the center name as in pData (Baylor); popstats uses the same names
    Fm = man.pivot_table(index="key", columns="Population", values="freq")
    Fm_t = man.pivot_table(index="key", columns="Population", values="freq_t").drop(columns=["Lymphocytes", "CD3"], errors="ignore")
    common = [k for k in meta.index if k in Fm.index]
    missing = [k for k in meta.index if k not in Fm.index]
    Fm = Fm.loc[common]; Fm_t = Fm_t.loc[common]; meta_m = meta.loc[common]
    tmask, cd3_cut = tcell_mask(root, z, arms)
    summary, clusters_out, match_out = [], [], []
    def summarize(arm, F, Ft):
        s = site_stats(F, meta_m if arm == "manual" else meta); s["arm"] = arm; s["denominator"] = "all events"
        st = site_stats(Ft, meta_m if arm == "manual" else meta); st["arm"] = arm; st["denominator"] = "T cells"
        clusters_out.extend([s, st])
        u, ut = s[s.mean_freq >= MIN_MEAN], st[st.mean_freq >= MIN_MEAN]
        return s, dict(arm=arm, n_clusters=F.shape[1], n_used=int(u.cluster.nunique()),
                       median_between_site_cv=float(u.between_site_cv.median()), median_within_site_cv=float(u.within_site_cv.median()),
                       **{f"between_cv_lot_{l}": float(u[u.lot == l].between_site_cv.median()) for l in sorted(meta.lot.unique())},
                       n_used_t=int(ut.cluster.nunique()), median_between_site_cv_t=float(ut.between_site_cv.median()),
                       median_within_site_cv_t=float(ut.within_site_cv.median()))
    sm, row = summarize("manual", Fm, Fm_t); summary.append(row)
    for arm in arms:
        pred = z[arm].astype(int)
        F = pd.crosstab(key, pred, normalize="index").loc[meta.index]          # all events of the subsample
        tm = tmask[arm]
        Ft = pd.crosstab(key[tm], pred[tm], normalize="index").reindex(meta.index).fillna(0)   # within the arm's T cell clusters
        s, row = summarize(arm, F, Ft)
        row["t_cell_fraction"] = float(tm.mean())
        # match manual populations to clusters by correlation across the files both have
        Fc = F.loc[common]; rs = []
        for p in Fm.columns:
            a = Fm[p].values
            cors = [float(np.corrcoef(a, Fc[c].values)[0, 1]) if Fc[c].std() > 0 and a.std() > 0 else np.nan for c in Fc.columns]
            j = int(np.nanargmax(cors)); c = Fc.columns[j]
            rs.append(cors[j])
            mcv = sm[sm.cluster == p].between_site_cv.median(); ccv = s[s.cluster == c].between_site_cv.median()
            match_out.append(dict(arm=arm, population=p, cluster=int(c), r=cors[j], manual_mean_freq=float(Fm[p].mean()),
                                  cluster_mean_freq=float(Fc[c].mean()), manual_between_cv=float(mcv), cluster_between_cv=float(ccv)))
        row.update(median_r_match=float(np.nanmedian(rs)))
        summary.append(row)
    S = pd.DataFrame(summary); C = pd.concat(clusters_out, ignore_index=True); M = pd.DataFrame(match_out)
    S.to_csv(res / f"{prefix}.csv", index=False); C.to_csv(res / f"{prefix}_clusters.csv", index=False)
    M.to_csv(res / f"{prefix}_match.csv", index=False)
    C[C.arm == "manual"].to_csv(res / f"{prefix}_manual.csv", index=False)

    def f(v, d=2): return "" if pd.isna(v) else f"{v:.{d}f}"
    lots = sorted(meta.lot.unique())
    print(f"## Between-site variance, {labels_file}, {meta.site.nunique()} centers, {len(lots)} lots, {len(meta)} files ({len(missing)} files without manual statistics)\n")
    print(f"CD3 cut on the deposited scale for the T cell denominator: {cd3_cut:.0f}\n")
    print("| Arm | Clusters | At 1 percent or more | Median between-site CV, all events | " + " | ".join(f"Lot {l}" for l in lots) + " | Median within-site CV, replicates | T cell fraction | Clusters at 1 percent of T cells | Median between-site CV, of T cells | Median within-site CV, of T cells | Median r to best cluster |")
    print("|---|---|---|---|" + "---|" * len(lots) + "---|---|---|---|---|---|")
    for _, r in S.iterrows():
        print(f"| {ARM_NAMES.get(r.arm, r.arm)} | {r.n_clusters} | {r.n_used} | {f(r.median_between_site_cv)} | " + " | ".join(f(r[f'between_cv_lot_{l}']) for l in lots) + f" | {f(r.median_within_site_cv)} | {f(r.get('t_cell_fraction'))} | {r.n_used_t} | {f(r.median_between_site_cv_t)} | {f(r.median_within_site_cv_t)} | {f(r.get('median_r_match'))} |")
    print("\n## Manual populations: between-site CV per lot (the reference)\n")
    piv = sm.pivot(index="cluster", columns="lot", values="between_site_cv"); mf = sm.groupby("cluster").mean_freq.mean()
    smt = C[(C.arm == "manual") & (C.denominator == "T cells")]; pivt = smt.pivot(index="cluster", columns="lot", values="between_site_cv")
    print("| Population | Mean frequency of all events | " + " | ".join(f"Lot {l}, all events" for l in lots) + " | Within-site CV | " + " | ".join(f"Lot {l}, of CD3" for l in lots) + " |")
    print("|---|---|" + "---|" * (2 * len(lots) + 1))
    for p in mf.sort_values(ascending=False).index:
        tcols = " | ".join(f(pivt.loc[p, l]) if p in pivt.index else "" for l in lots)
        print(f"| {p} | {f(mf[p], 3)} | " + " | ".join(f(piv.loc[p, l]) for l in lots) + f" | {f(sm[sm.cluster == p].within_site_cv.mean())} | {tcols} |")
    print("\n## Best-matching cluster per manual population, per arm (r across files; between-site CV of the cluster against the population)\n")
    arms_p = [a for a in ORDER if a in set(M.arm)]
    print("| Population | Manual CV | " + " | ".join(f"{a} r / CV" for a in arms_p) + " |")
    print("|---|---|" + "---|" * len(arms_p))
    for p in mf.sort_values(ascending=False).index:
        cells = []
        for a in arms_p:
            m = M[(M.arm == a) & (M.population == p)]
            cells.append(f"{f(m.r.values[0])} / {f(m.cluster_between_cv.values[0])}" if len(m) else "")
        print(f"| {p} | {f(sm[sm.cluster == p].between_site_cv.median())} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    lab = sys.argv[sys.argv.index("--labels") + 1] if "--labels" in sys.argv else "results/labels_hipc_HIPC_tcell_k30.npz"
    pre = sys.argv[sys.argv.index("--prefix") + 1] if "--prefix" in sys.argv else "hipc_between_site"
    main(sys.argv[1], lab, pre)
