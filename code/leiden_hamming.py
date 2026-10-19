"""Leiden on a kNN graph, with a distance that suits each representation, on a 50,000-cell
subsample of Levine_32dim.

FlowSOM uses Euclidean distance on every arm. A binary arm has a natural distance, Hamming,
so this script asks whether the binary arms lose under FlowSOM because of the representation or
because of the clusterer.

Graphs
  A1            Euclidean on arcsinh(x / 5), the scanpy default
  C_gate        Euclidean on the mixed matrix (0/1 lineage columns plus z-scored functional
                columns), and Hamming on the 25 lineage columns alone. Hamming on a mixed
                matrix is not defined, because every continuous column counts as a mismatch.
  D_gate        Hamming on all 32 binary columns. Euclidean on 0/1 data ranks neighbors the
                same way as Hamming, so one run covers both.
  D_gate_clean  Hamming on the binary columns whose gate cut reaches balanced accuracy 0.9.

Leiden has no cluster count. The script runs a resolution grid and reports the number of
clusters each resolution finds, so the reader can pick the rows nearest 14 and 40. Stability is
the mean pairwise ARI over 10 Leiden seeds on the same graph.

Usage: python3 leiden_hamming.py <project_dir>
Writes results/leiden_metrics.csv and results/leiden_perpop.csv.
"""
import sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pilot import load_levine32, load_gate_thresholds, build_arms, hungarian_f1, SUB_N

warnings.filterwarnings("ignore")
RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5]
SEEDS = list(range(1, 11))
N_NEIGHBORS = 15


def graph(X, metric, seed=0):
    a = ad.AnnData(X=np.ascontiguousarray(X, dtype=np.float32))
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, use_rep="X", metric=metric, random_state=seed)
    return a


def main(root):
    root = Path(root)
    ds = load_levine32(root)
    gate = load_gate_thresholds(root, ds)
    assert gate is not None, "run code/gate_thresholds.py first"
    arms, thr, lin, _ = build_arms(ds, gate)
    rng = np.random.default_rng(20261018)
    sub = rng.choice(len(ds["X"]), SUB_N, replace=False)
    truth = ds["labels"][sub]; labeled = truth > 0
    jobs = [("A1", "euclidean", arms["A1"][sub]),
            ("C_gate", "euclidean", arms["C_gate"][sub]),
            ("C_gate", "hamming", arms["C_gate"][sub][:, lin]),
            ("D_gate", "hamming", arms["D_gate"][sub]),
            ("D_gate_clean", "hamming", arms["D_gate_clean"][sub])]
    rows, perpop = [], []
    for arm, metric, X in jobs:
        t0 = time.time()
        a = graph(X, metric)
        for r in RESOLUTIONS:
            labs = []
            for s in SEEDS:
                sc.tl.leiden(a, resolution=r, random_state=s, key_added="l", flavor="igraph", n_iterations=2, directed=False)
                labs.append(a.obs["l"].astype(int).to_numpy())
            pred = labs[0]
            f1, per = hungarian_f1(truth[labeled], pred[labeled])
            pair = [adjusted_rand_score(labs[i], labs[j]) for i in range(len(labs)) for j in range(i + 1, len(labs))]
            row = dict(arm=arm, metric=metric, n_features=X.shape[1], resolution=r, n_clusters=int(len(np.unique(pred))),
                       mean_f1=f1, ari=adjusted_rand_score(truth[labeled], pred[labeled]),
                       nmi=normalized_mutual_info_score(truth[labeled], pred[labeled]), stability_ari=float(np.mean(pair)))
            rows.append(row)
            for p, v in per.items():
                perpop.append(dict(arm=arm, metric=metric, resolution=r, population=ds["pop_names"][p], f1=v))
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        print(f"{arm} {metric} done in {time.time() - t0:.0f} s", flush=True)
    pd.DataFrame(rows).to_csv(root / "results/leiden_metrics.csv", index=False)
    pd.DataFrame(perpop).to_csv(root / "results/leiden_perpop.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
