"""Reproduce the deposited (compensated, transformed) HIPC Stanford T cell files from the raw
Diva exports, and fit the transform the depositor applied.

What the depositor did (established 2026-10-19 with this script): flowWorkspace parsed the
FlowJo 9.6.3 central-analysis workspace, compensated with the matrix the analyst edited in
FlowJo (the first block of `XML/<Center>/<Center> Matrices/Tcell`, not the acquisition matrix
and not `FCS/compensation/<Center>.csv`), and applied FlowJo's biexponential, which is a
logicle with T 262144, M 4.5 decades, A 0, scaled to 4096, and a width W chosen per parameter.
The README's iplogicle (T 16777215, M 7.22, W 1.16) does not reproduce the files.

For each fluorescence channel the script searches W on a grid and reports the best width and
the median absolute error against the deposited values, event by event (flowCore keeps event
order). Run it on two files to see that the widths are a property of the workspace, not of
the file.

Usage: python3 hipc_transform_check.py <project_dir> [raw_file_name ...]
"""
import re, sys
from pathlib import Path
import numpy as np, flowio
from flowutils import transforms as FT

DETECTOR_TO_MARKER = {"FITC-A": "Live", "PE-A": "CCR7", "PerCP-Cy5-5-A": "CD4", "PE-Cy7-A": "CD45RA",
                      "Pacific Blue-A": "CD3", "Am Cyan-A": "HLADR", "APC-A": "CD38", "APC-Cy7-A": "CD8"}
T, M, A, SCALE = 262144.0, 4.5, 0.0, 4096.0


def read_fj_matrix(path):
    """FlowJo 9 'Matrices/<panel>' export, CR line endings: name, '<\\t>', the detector header, the
    k rows of the matrix the workspace used, then 'Original' and the acquisition matrix."""
    txt = open(path, encoding="latin-1").read().replace("\r\n", "\n").replace("\r", "\n")
    lines = [l for l in txt.split("\n") if l.strip()]
    hdr = lines[2].split("\t"); k = len(hdr); num = re.compile(r"^-?[\d.]+(e-?\d+)?$")
    def block(start):
        rows = []; i = start
        while len(rows) < k:
            cells = lines[i].split("\t"); i += 1
            if all(num.match(c.strip()) for c in cells[:k]): rows.append([float(c) for c in cells[:k]])
        return np.array(rows)
    orig_i = next(i for i, l in enumerate(lines) if l.startswith("Original"))
    return hdr, block(3), block(orig_i + 1)


def read(p):
    f = flowio.FlowData(str(p)); t = f.text; n = f.channel_count
    return [t.get(f"p{i}n", t.get(f"$P{i}N")) for i in range(1, n + 1)], np.reshape(f.events, (-1, n)).astype(np.float64), t


def main(root, names):
    root = Path(root)
    hdr, Mws, Morig = read_fj_matrix(root / "data/HIPC/XML/Stanford/Stanford Matrices/Tcell")
    print("Workspace matrix differs from the acquisition matrix in these cells (row spills into column):")
    for i, r in enumerate(hdr):
        for j, c in enumerate(hdr):
            if abs(Mws[i, j] - Morig[i, j]) > 0.005: print(f"  {r} -> {c}: workspace {Mws[i, j]:.3f}, acquisition {Morig[i, j]:.3f}")
    Ws = np.round(np.arange(0.05, 2.05, 0.05), 2)
    print("\n| File | Marker | Best W | Median abs error, 4096 scale | r |\n|---|---|---|---|---|")
    for name in names:
        pr, R, _ = read(root / "data/HIPC/RAW/Stanford" / name); pdd, D, _ = read(root / "data/HIPC/FCS/Stanford" / name)
        assert len(R) == len(D), "event counts differ"
        X = R[:, [pr.index(d) for d in hdr]] @ np.linalg.inv(Mws)
        sub = np.random.default_rng(0).choice(len(X), min(30000, len(X)), replace=False)
        for j, d in enumerate(hdr):
            b = D[sub, pdd.index(DETECTOR_TO_MARKER[d])]; ok = b > -4096; x = X[sub, j]
            best = min(((np.median(np.abs(FT.logicle(x.reshape(-1, 1), [0], t=T, m=M, w=W, a=A).ravel() * SCALE - b)[ok]), W) for W in Ws))
            a = FT.logicle(x.reshape(-1, 1), [0], t=T, m=M, w=best[1], a=A).ravel() * SCALE
            print(f"| {name} | {DETECTOR_TO_MARKER[d]} | {best[1]:.2f} | {best[0]:.1f} | {np.corrcoef(a[ok], b[ok])[0, 1]:.5f} |")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or ["1228-1_A1_A01.fcs", "1369-3_A9_A09.fcs"])
