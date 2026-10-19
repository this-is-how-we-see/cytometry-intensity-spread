"""How much does arm B on HIPC depend on the per-channel logicle widths?

Only Stanford deposited raw files, so the widths fitted there (HIPC_LOGICLE["w"]) are applied to
every center. This script measures the cost of that assumption on Stanford itself: arm B built
with the fitted widths against arm B built with one flat width (0.5 everywhere) and with the
widths all doubled, clustered the same way (FlowSOM, 12 metaclusters, seed 1), compared by ARI.
It also reports the ARI of arm B against A1, so the width effect can be read against the effect
arm B itself has.

Usage: python3 hipc_width_check.py <project_dir>
"""
import sys
from pathlib import Path
import numpy as np
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).parent)); import pilot

root = Path(sys.argv[1])
ds = pilot.load_hipc_tcell(root, per_file=20_000, centers=("Stanford",))
variants = {"fitted": dict(pilot.HIPC_LOGICLE), "flat 0.5": dict(pilot.HIPC_LOGICLE, w=0.5),
            "doubled": dict(pilot.HIPC_LOGICLE, w={k: min(2 * v, 2.0) for k, v in pilot.HIPC_LOGICLE["w"].items()})}
labels = {}
for name, L in variants.items():
    arms, thr, lin, div = pilot.build_arms(dict(ds, logicle=L))
    if "A1" not in labels: labels["A1"] = pilot.cluster(arms["A1"], 12, seed=1)
    labels[name] = pilot.cluster(arms["B_geom"], 12, seed=1)
    print(f"{name:9} built; B_geom CD3 median {np.median(arms['B_geom'][:, ds['markers'].index('CD3')]):.0f}", flush=True)
print("\n| Comparison | ARI |\n|---|---|")
for a, b in [("fitted", "flat 0.5"), ("fitted", "doubled"), ("A1", "fitted"), ("A1", "flat 0.5")]:
    print(f"| B_geom {a} vs {b} | {adjusted_rand_score(labels[a], labels[b]):.3f} |")
