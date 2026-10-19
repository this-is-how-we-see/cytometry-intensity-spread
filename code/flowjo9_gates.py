"""Apply the HIPC Lyoplate central manual gates (FlowJo 9.6.3 Mac workspaces, "Save as XML",
Workspace version 2.0) to the deposited T cell panel FCS files, so that every event gets a
manual label and the dataset can be scored with Hungarian F1 like the others.

What the workspace holds, and how this module reads it (established 2026-10-19):

* `SampleList/Sample`: one per FCS file, with `Keywords` ($FIL is the file name), `Parameter`
  elements (range 4096, gain 0.015625 for scatter, log 1 and calibrationIndex 1 for the
  fluorescence channels) and a `SampleNode` that holds the gating tree as nested `Population`
  elements. Each Population has one `PolygonGate` with a `Polygon` or `PolyRect` child that
  names the axes (`xAxisName`, `yAxisName`; angle brackets mean the compensated parameter),
  carries `negated="1"` for a NOT gate, and lists `Vertex` elements.
* Vertex coordinates are on the linear data scale, not in display units: a scatter vertex is
  0 to 262144 and a fluorescence vertex is the linear value under the display transform
  (so a rectangle edge dragged past the plot border has a value like -2.8e5 or 8.6e7).
* `CalibrationTables/Table`: one table per group and compensated channel, named like
  "Tcell <PE-A>", with 4097 pairs (display channel 0..4096, linear value). This is the
  biexponential FlowJo used for that channel in that group (the attributes give its width
  basis), and it is what flowWorkspace applied when it wrote the deposited files. A gate
  vertex maps into display space by inverse interpolation of the table; outside the table
  the ends are extended linearly, which only matters for edges far past the plot border.
* `Groups/GroupNode[name="Tcell"]/Group/SampleRefs`: which samples belong to the T cell panel.

The deposited files are already compensated and in display space on the 4096 scale
(`data/README.md`, HIPC section), with FSC-A and SSC-A linear on 262144 and the fluorescence
columns renamed to markers. So gating is: map each vertex to display space, then test the
deposited values. Scatter gates test the deposited scatter directly. The singlets gate needs
FSC-H or FSC-W, which the deposited files do not carry; when a raw file of the same name sits
in `data/HIPC/RAW/<Center>/` with the same event count (Stanford only), the missing scatter
columns come from it, otherwise the singlets gate is skipped and reported.

Labels: 1..8 for the eight memory subsets in `LEAF_NAMES` order (CD4 then CD8: naive,
central memory, effector memory, effector), 0 for every other event (debris, doublets, dead,
non-T, double positive and double negative T cells). `activated` is a separate flag (CD38+
HLA-DR+ within CD4 or CD8). The check is the central manual population statistics in
`data/HIPC/popstats_manual/popstats-tcell.csv`, which flowWorkspace computed from these same
gates.

Usage:
  python -I code/flowjo9_gates.py <project_dir> [--centers Stanford,Baylor] [--files name,...]
      [--write] [--report results/hipc_gate_agreement.csv]
"""
import argparse, re, sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np, pandas as pd

CENTERS = {"Baylor": "BIIR", "CIMR": "CIMR", "Miami": "Miami", "NHLBI": "NHLBI", "Stanford": "Stanford", "UCLA": "UCLA", "Yale": "Yale"}
GROUP = "Tcell"
# manual population path -> common name of popstats-tcell.csv (the HIPC mapping-tcell.r is not
# on disk; this mapping is read off the gate names and confirmed by the counts)
QUAD = {"CCR7+ 45RA+": "Naive", "CCR7+ 45RA-": "Central Memory", "CCR7- 45RA-": "Effector Memory", "CCR7- 45RA+": "Effector"}
LEAF_NAMES = ["CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", "CD4 Effector",
              "CD8 Naive", "CD8 Central Memory", "CD8 Effector Memory", "CD8 Effector"]
LINEAGE = {"4+ 8-": "CD4", "4- 8+": "CD8"}
# detector name (without brackets, FlowJo spelling per center) -> marker of the fixed Lyoplate panel
DETECTOR_RULES = [
    (r"FITC|Alexa Fluor 488|B515", "Live"),
    (r"PE[- ]?Cy7|G780", "CD45RA"),
    (r"PerCP|B710", "CD4"),
    (r"^PE\b|PE YG|G560", "CCR7"),
    (r"Pacific Blue|V450", "CD3"),
    (r"Am ?Cyan|V500|V545", "HLADR"),
    (r"APC[- ]?(H7|Cy7)|R780", "CD8"),
    (r"^APC|R660", "CD38"),
]


def detector_marker(axis):
    name = axis.strip("<>")
    if re.match(r"^(FSC|SSC)-[AHW]$", name): return name
    for pat, m in DETECTOR_RULES:
        if re.search(pat, name): return m
    raise KeyError(f"no marker for detector {axis!r}")


# ----------------------------------------------------------------------------- parsing
def parse_workspace(path):
    root = ET.parse(path).getroot()
    tables = {}
    for t in root.find("CalibrationTables"):
        v = np.array([float(x) for x in t.text.replace("\n", "").split(",") if x.strip()])
        tables[t.get("name")] = dict(channel=v[0::2], value=v[1::2], width=float(t.get("biexponentialWidth", "nan")),
                                     decades=float(t.get("biexponentialDecades", "nan")))
    groups = {}
    for g in root.find("Groups"):
        refs = g.find("Group").find("SampleRefs")
        groups[g.get("name")] = {s.get("sampleID") for s in refs} if refs is not None else set()
    samples = []
    for s in root.find("SampleList"):
        kw = {k.get("name"): k.get("value") for k in s.find("Keywords")}
        node = s.find("SampleNode")
        if node is None: continue
        pops = []
        _walk(node, (), pops)
        samples.append(dict(sample_id=s.get("sampleID"), file=kw.get("$FIL"), event_count=int(s.get("eventCount", 0)),
                            groups={g for g, ids in groups.items() if s.get("sampleID") in ids},
                            params={p.get("name"): dict(p.attrib) for p in s.findall("Parameter")}, populations=pops))
    return dict(version=root.get("flowJoVersion"), tables=tables, groups=groups, samples=samples)


def _walk(elem, path, out):
    for p in elem.findall("Population"):
        name = p.get("name"); here = path + (name,)
        g = p.find("PolygonGate")
        shape = None
        for tag in ("Polygon", "PolyRect"):
            e = None if g is None else g.find(tag)
            if e is not None:
                shape = dict(type=tag, x_axis=e.get("xAxisName"), y_axis=e.get("yAxisName"), negated=e.get("negated") == "1",
                             vertices=np.array([[float(v.get("x")), float(v.get("y"))] for v in e.iter("Vertex")]))
                break
        if shape is None:   # BooleanGate and the like, used in the other panels; refused at gating time
            shape = dict(type=[c.tag for c in p if c.tag.endswith("Gate")], x_axis=None, y_axis=None, negated=False, vertices=None)
        out.append(dict(path=here, name=name, parent=path, count=int(p.get("count", -1)), **shape))
        _walk(p, here, out)


# ----------------------------------------------------------------------------- transform
def group_table(tables, axis, group=GROUP):
    """The group's table for a compensated axis: named "Tcell <PE-A>" at six centers and
    "Tcell wizard1 <G560-A>" at NHLBI, so match on the group prefix and the axis suffix."""
    hits = [k for k in tables if k.startswith(group) and k.endswith(f" {axis}")]
    if len(hits) != 1: raise KeyError(f"{len(hits)} calibration tables for {group} {axis}: {hits}")
    return tables[hits[0]]


def to_display(values, axis, tables, group=GROUP):
    """Linear vertex coordinate -> display channel (0..4096 scale) for a compensated axis,
    by inverse interpolation of the group's calibration table; scatter passes through."""
    if not axis.startswith("<"): return np.asarray(values, float)
    t = group_table(tables, axis, group); ch, val = t["channel"], t["value"]
    v = np.asarray(values, float); out = np.interp(v, val, ch)
    lo, hi = v < val[0], v > val[-1]
    out[lo] = ch[0] + (v[lo] - val[0]) * (ch[1] - ch[0]) / (val[1] - val[0])
    out[hi] = ch[-1] + (v[hi] - val[-1]) * (ch[-1] - ch[-2]) / (val[-1] - val[-2])
    return out


def from_display(display, axis, tables, group=GROUP):
    """Display channel -> linear value through the table (the forward direction FlowJo used)."""
    t = group_table(tables, axis, group); ch, val = t["channel"], t["value"]
    d = np.asarray(display, float); out = np.interp(d, ch, val)
    lo, hi = d < ch[0], d > ch[-1]
    out[lo] = val[0] + (d[lo] - ch[0]) * (val[1] - val[0]) / (ch[1] - ch[0])
    out[hi] = val[-1] + (d[hi] - ch[-1]) * (val[-1] - val[-2]) / (ch[-1] - ch[-2])
    return out


# ----------------------------------------------------------------------------- gating
def in_polygon(x, y, verts):
    """Even-odd ray casting, vectorised over events."""
    inside = np.zeros(len(x), bool)
    n = len(verts)
    for i in range(n):
        x0, y0 = verts[i]; x1, y1 = verts[(i + 1) % n]
        cross = (y0 > y) != (y1 > y)          # a horizontal edge never crosses, so its xi is unused
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
        inside ^= cross & (x < xi)
    return inside


def _edges(v, axis):
    top = 4096.0 if axis.startswith("<") else 262144.0
    return (-np.inf if v.min() < 0 else v.min()), (np.inf if v.max() > top else v.max())


def gate_mask(pop, columns, tables):
    """Boolean mask of events inside this population's own gate (parent not applied).
    `columns` maps a marker or scatter name to its display-space array; None means the
    axis is unavailable and the gate cannot be applied."""
    if pop["vertices"] is None: raise ValueError(f"{'/'.join(pop['path'])}: gate type {pop['type']} is not supported")
    xm, ym = detector_marker(pop["x_axis"]), detector_marker(pop["y_axis"])
    if columns.get(xm) is None or columns.get(ym) is None: return None
    vx = to_display(pop["vertices"][:, 0], pop["x_axis"], tables)
    vy = to_display(pop["vertices"][:, 1], pop["y_axis"], tables)
    x, y = columns[xm], columns[ym]
    if pop["type"] == "PolyRect":
        # an edge dragged past the plot border means no limit: FlowJo piles events beyond the
        # display range onto the axis, so they fall inside such a gate (checked against popstats
        # on Baylor, where the depositor's export put CD45RA values below the plot)
        (xlo, xhi), (ylo, yhi) = (_edges(vx, pop["x_axis"]), _edges(vy, pop["y_axis"]))
        m = (x >= xlo) & (x <= xhi) & (y >= ylo) & (y <= yhi)
    else:
        m = in_polygon(x, y, np.column_stack([vx, vy]))
    return ~m if pop["negated"] else m


def apply_tree(sample, columns, tables):
    """Masks for every population path, parent applied. Returns (masks, skipped) where skipped
    lists populations whose axes the data lacks (their mask is the parent's, pass-through)."""
    masks, skipped = {(): np.ones(len(next(v for v in columns.values() if v is not None)), bool)}, []
    for pop in sample["populations"]:
        own = gate_mask(pop, columns, tables)
        if own is None: skipped.append(pop["path"]); own = np.ones_like(masks[pop["parent"]])
        masks[pop["path"]] = masks[pop["parent"]] & own
    return masks, skipped


def find_path(masks, *parts):
    """The one population path whose last elements are `parts` (names repeat under CD4 and CD8)."""
    hits = [p for p in masks if p[-len(parts):] == parts]
    if len(hits) != 1: raise KeyError(f"{parts}: {len(hits)} matches")
    return hits[0]


def label_events(masks):
    """label 1..8 per LEAF_NAMES, 0 elsewhere; activated flag; the common-name count table."""
    n = len(next(iter(masks.values()))); label = np.zeros(n, np.int8); activated = np.zeros(n, bool)
    counts = {"Lymphocytes": int(masks[find_path(masks, "LYM")].sum()), "Live": int(masks[find_path(masks, "not dead")].sum()),
              "CD3": int(masks[find_path(masks, "CD3")].sum())}
    for gate, lin in LINEAGE.items():
        counts[lin] = int(masks[find_path(masks, gate)].sum())
        for q, sub in QUAD.items():
            m = masks[find_path(masks, gate, q)]; name = f"{lin} {sub}"
            label[m] = LEAF_NAMES.index(name) + 1; counts[name] = int(m.sum())
        m = masks[find_path(masks, gate, "38+ DR+")]; activated |= m; counts[f"{lin} Activated"] = int(m.sum())
    return label, activated, counts


# ----------------------------------------------------------------------------- data
def read_fcs(path):
    import flowio
    f = flowio.FlowData(str(path)); t = f.text; n = f.channel_count
    names = [t.get(f"p{i}n", t.get(f"$P{i}N")) for i in range(1, n + 1)]
    return names, np.reshape(f.events, (-1, n)).astype(np.float64)


def columns_for(center, fname, root):
    """Display-space columns of a deposited file, keyed by marker and scatter name; FSC-H and
    FSC-W from the raw file of the same name when it exists with the same event count."""
    names, D = read_fcs(root / "data/HIPC/FCS" / center / fname)
    cols = {n: D[:, i] for i, n in enumerate(names)}
    raw = root / "data/HIPC/RAW" / center / fname
    source = "deposited"
    if raw.exists():
        rn, R = read_fcs(raw)
        if len(R) == len(D):
            for n in ("FSC-H", "FSC-W", "SSC-H", "SSC-W"):
                if n in rn and n not in cols: cols[n] = R[:, rn.index(n)]
            source = "deposited + raw scatter"
    return cols, source


def tcell_samples(ws, center, root):
    pdata = pd.read_csv(root / f"data/HIPC/pData/{center}.csv"); pdata.columns = ["file", "institution", "panel", "replicate", "sample"]
    tfiles = set(pdata[pdata.panel == "T cell"].file)
    out = [s for s in ws["samples"] if s["file"] in tfiles and GROUP in s["groups"]]
    missing = tfiles - {s["file"] for s in out}
    return out, missing


def check_axes(sample):
    """Every gate must sit on the detectors its name implies; raise otherwise."""
    expect = {"not dead": ("Live", "SSC-A"), "CD3": ("CD3", "SSC-A"), "4+ 8-": ("CD4", "CD8"), "4- 8+": ("CD4", "CD8"),
              "DPT": ("CD4", "CD8"), "DNT": ("CD4", "CD8"), "LYM": ("FSC-A", "SSC-A")}
    for q in QUAD: expect[q] = ("CCR7", "CD45RA")
    for q in ("38- DR+", "38+ DR+", "38+ DR-", "38- DR-"): expect[q] = ("CD38", "HLADR")
    for p in sample["populations"]:
        got = (detector_marker(p["x_axis"]), detector_marker(p["y_axis"]))
        if p["name"] in expect and got != expect[p["name"]]:
            raise ValueError(f"{sample['file']} {'/'.join(p['path'])}: axes {p['x_axis']}, {p['y_axis']} read as {got}, expected {expect[p['name']]}")


# ----------------------------------------------------------------------------- main
def run(root, centers, files=None, write=False, report=None, verbose=True):
    pop = pd.read_csv(root / "data/HIPC/popstats_manual/popstats-tcell.csv")
    rows = []
    for center in centers:
        ws = parse_workspace(root / "data/HIPC/XML" / CENTERS[center] / f"CA c3 v2 {CENTERS[center]}.xml")
        samples, missing = tcell_samples(ws, center, root)
        if missing and verbose: print(f"{center}: T cell files with no sample in the Tcell group: {sorted(missing)}", file=sys.stderr)
        for s in samples:
            if files and s["file"] not in files: continue
            if not (root / "data/HIPC/FCS" / center / s["file"]).exists(): continue
            check_axes(s)
            cols, source = columns_for(center, s["file"], root)
            if len(next(iter(cols.values()))) != s["event_count"]:
                raise ValueError(f"{center} {s['file']}: {len(next(iter(cols.values())))} events on disk, workspace says {s['event_count']}")
            masks, skipped = apply_tree(s, cols, ws["tables"])
            label, activated, counts = label_events(masks)
            ref = pop[(pop.Center == center) & (pop.Filename == s["file"])]
            ref_counts = dict(zip(ref.Population, ref.Count))
            live_parent = ref[ref.Parent == "Live"].Count_Parent
            if len(live_parent): ref_counts["Live"] = int(live_parent.iloc[0])
            xml_counts = {}
            for p in s["populations"]:
                if p["name"] in ("LYM", "not dead", "CD3", "4+ 8-", "4- 8+") or p["name"] in QUAD or p["name"] == "38+ DR+":
                    key = {"LYM": "Lymphocytes", "not dead": "Live", "CD3": "CD3", "4+ 8-": "CD4", "4- 8+": "CD8"}.get(p["name"])
                    if key is None:
                        lin = LINEAGE.get(p["path"][-2]); key = f"{lin} {QUAD.get(p['name'], 'Activated')}"
                    xml_counts[key] = p["count"]
            for name, n in counts.items():
                r = ref_counts.get(name); x = xml_counts.get(name)
                rows.append(dict(center=center, file=s["file"], population=name, ours=n, popstats=r, flowjo=x,
                                 rel_err=(n - r) / r if r else np.nan, source=source,
                                 singlets_skipped=any(p[-1] == "singlets" for p in skipped)))
            if verbose:
                cd3 = next(r for r in rows[-len(counts):] if r["population"] == "CD3")
                print(f"{center:9s} {s['file']:32s} {source:24s} CD3 ours {cd3['ours']:7d} popstats {cd3['popstats']:7d} "
                      f"({100 * cd3['rel_err']:+.2f} %)" + (f"  skipped {['/'.join(p) for p in skipped]}" if skipped else ""), flush=True)
            if write:
                out = root / "data/HIPC/labels" / center; out.mkdir(parents=True, exist_ok=True)
                flags = {k: masks[find_path(masks, k)] for k in ("LYM", "not dead", "CD3", "4+ 8-", "4- 8+", "DPT", "DNT")}
                np.savez_compressed(out / f"{s['file']}.npz", label=label, activated=activated,
                                    pop_names=np.array(["other"] + LEAF_NAMES), source=source,
                                    singlets_applied=not any(p[-1] == "singlets" for p in skipped),
                                    **{k.replace(" ", "_").replace("+", "p").replace("-", "n"): v for k, v in flags.items()})
    df = pd.DataFrame(rows)
    if report: df.to_csv(root / report, index=False)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root"); ap.add_argument("--centers", default=",".join(CENTERS)); ap.add_argument("--files", default=None)
    ap.add_argument("--write", action="store_true"); ap.add_argument("--report", default=None)
    a = ap.parse_args()
    df = run(Path(a.root), a.centers.split(","), set(a.files.split(",")) if a.files else None, a.write, a.report)
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
