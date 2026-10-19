"""Extract expression matrices and name vectors from HDCytoData .rda files without R.

HDCytoData stores each dataset as a SummarizedExperiment (cells in rows, markers in columns)
with rowData holding population_id and colData holding channel_name, marker_name, marker_class.
The assay lives inside an R environment, so this walks the parsed object tree generically.
Usage: python3 rda_extract.py <file.rda> <out.npz>
"""
import sys, dataclasses
import numpy as np, rdata
from rdata.parser import RObjectType as T

def attr_dict(o):
    """Return the attributes pairlist of an RObject as {name: RObject}."""
    d = {}
    a = o.attributes
    while a is not None and getattr(a, "info", None) is not None and a.info.type == T.LIST:
        if a.tag is not None:
            key = a.tag
            while hasattr(key, "value"): key = key.value
            if isinstance(key, bytes): key = key.decode()
            d[key] = a.value[0]
        a = a.value[1] if isinstance(a.value, (list, tuple)) and len(a.value) == 2 else None
    return d

def strvec(o):
    return [v.value.decode() if isinstance(v.value, bytes) else v.value for v in o.value]

def walk(o, found, depth=0):
    if depth > 60: return
    if hasattr(o, "info"):
        t = o.info.type
        if t == T.REAL and o.value is not None and len(o.value) >= 1000:
            dims = attr_dict(o).get("dim")
            found["real"].append((np.asarray(o.value), list(dims.value) if dims is not None else None))
        elif t == T.INT and o.value is not None and len(o.value) >= 1000:
            found["int"].append(np.asarray(o.value))
        elif t == T.STR and o.value is not None and 3 <= len(o.value) <= 200:
            found["str"].append(strvec(o))
        for name in ("value", "attributes", "tag"):
            walk(getattr(o, name, None), found, depth + 1)
    elif dataclasses.is_dataclass(o):
        for f in dataclasses.fields(o):
            walk(getattr(o, f.name), found, depth + 1)
    elif isinstance(o, (list, tuple)):
        for c in o: walk(c, found, depth + 1)

def extract(path):
    parsed = rdata.parser.parse_file(path)
    found = {"real": [], "int": [], "str": []}
    walk(parsed.object, found)
    mats = [(m, d) for m, d in found["real"] if d is not None and len(d) == 2]
    mats.sort(key=lambda x: -x[0].size)
    if mats:
        mat, dims = mats[0]
    else:  # fall back: largest REAL vector, rows from the population_id vector
        mat = max((m for m, d in found["real"]), key=lambda m: m.size)
        nrow = max(len(a) for a in found["int"])
        dims = [nrow, mat.size // nrow]
    expr = mat.reshape(dims[1], dims[0]).T  # R is column-major
    strs = found["str"]
    # the colData vectors: channel_name, marker_name (same length as ncol), marker_class levels
    ncol = dims[1]
    cands = [s for s in strs if len(s) == ncol]
    import re
    def channel_score(s):  # fluorescence channels end in -A/-H/-W; mass channels look like (Pd102)Di
        return sum(x.endswith(("-A", "-H", "-W")) or bool(re.match(r"\(\w+\d+\)D[id]", x)) for x in s)
    cands.sort(key=lambda s: -channel_score(s))
    channel = cands[0]
    marker = next(s for s in cands if s != channel)
    # row data: every integer vector with one entry per cell. The population id is the one with
    # the fewest distinct values above 1 that is not a plain sample index; keep all of them anyway.
    rowvecs = [a for a in found["int"] if len(a) == dims[0]]
    pop = None
    if rowvecs:
        scored = sorted(rowvecs, key=lambda a: len(np.unique(a)))
        # prefer a vector whose values include 0 (unassigned) or whose distinct count exceeds 10
        pop = next((a for a in scored if 0 in np.unique(a) or len(np.unique(a)) > 10), scored[0])
    return expr, channel, marker, pop, found, rowvecs

if __name__ == "__main__":
    expr, channel, marker, pop, found, rowvecs = extract(sys.argv[1])
    print("expr", expr.shape, "channels", channel, "markers", marker)
    print("row vectors (distinct counts):", [len(np.unique(a)) for a in rowvecs])
    if pop is not None:
        v, c = np.unique(pop, return_counts=True); print("population_id", dict(zip(v.tolist(), c.tolist())))
    print("other str vectors:", [s[:6] for s in found["str"] if len(s) != len(channel)][:8])
    extra = {f"rowvec{i}": a for i, a in enumerate(rowvecs)}
    np.savez_compressed(sys.argv[2], expr=expr, channel=np.array(channel), marker=np.array(marker), pop=pop if pop is not None else np.array([]), **extra)
