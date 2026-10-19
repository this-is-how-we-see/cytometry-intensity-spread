"""Tables 1 to 4 for the preprint, as markdown, from the metrics files. Usage: python3 preprint_tables.py <project_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ARM = {"A1": "A1 arcsinh", "A2": "A2 z-score", "A3": "A3 percentile", "C": "C lineage bits + percentile", "C2": "C2 lineage bits + z-score",
       "D": "D all bits", "C_gate": "C2 with gate thresholds", "D_gate": "D with gate thresholds", "D_gate_clean": "D with gate thresholds, 15 clean markers"}
SETS = [("Levine_32dim", 14, "Levine_32dim_k40"), ("Levine_13dim", 24, "Levine_13dim_k40"), ("Samusik_01", 24, "Samusik_01_k40"), ("Samusik_all", 24, "Samusik_all_k40")]


def f(v, d=2):
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.{d}f}"


def metrics(res, tags):
    ms = [pd.read_csv(res / f"{t}_metrics.csv") for t in tags if (res / f"{t}_metrics.csv").exists()]
    m = pd.concat(ms, ignore_index=True)
    return m[~m.arm.str.contains("~")].drop_duplicates(["dataset", "arm"], keep="last").set_index(["dataset", "arm"])


def table1():
    rows = [
        ("Levine_32dim", "Mass (CyTOF)", "Human bone marrow, 2 donors", "265,627", "104,184", "32", "14", "Levine et al. [3]; Bendall et al. [5]"),
        ("Levine_13dim", "Mass (CyTOF)", "Human bone marrow, 1 donor", "167,044", "81,747", "13", "24", "Levine et al. [3]; Bendall et al. [5]"),
        ("Samusik_01", "Mass (CyTOF)", "Mouse bone marrow, 1 mouse", "86,864", "53,173", "39", "24", "Samusik et al. [4]"),
        ("Samusik_all", "Mass (CyTOF)", "Mouse bone marrow, 10 mice", "841,644", "514,386", "39", "24", "Samusik et al. [4]"),
        ("Nilsson_rare", "Flow, 13 colors", "Human bone marrow, CD34-enriched", "44,140", "358 rare", "13", "1", "Rundberg Nilsson et al. [6]"),
        ("Mosmann_rare", "Flow, BD LSR-II, 15 colors", "Human PBMC, influenza-stimulated", "396,460", "109 rare", "14", "1", "Mosmann et al. [7]"),
        ("HIPC Lyoplate T cell", "Flow, BD LSR II and LSRFortessa, 7 centers", "SeraCare PBMC, 3 lots, triplicates", "2,520,000 sampled from 63 files", "central gates", "7 plus viability", "14", "Finak et al. [8]"),
        ("DURAClone controls", "Flow, BD LSRFortessa", "Capture beads, one conjugate per tube", "12 tubes", "none", "12", "none", "Provided in kind by the University of Florida"),
    ]
    out = ["**Table 1. Datasets.** The six benchmark sets are FlowRepository FR-FCM-ZZPH as distributed by Weber and Robinson [1] and HDCytoData; the HIPC files are on ImmuneSpace, project HIPC/Lyoplate.", "",
           "| Dataset | Cytometry | Sample | Cells | Labeled | Markers | Populations | Source |", "|---|---|---|---|---|---|---|---|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def table2():
    rows = [
        ("A1", "arcsinh(x / c), c = 5 for mass and 150 for fluorescence", "Nothing; the control", "Euclidean"),
        ("A2", "A1, then each channel centered and divided by its SD", "Differences in spread between channels", "Euclidean"),
        ("A3", "A1, then each channel scaled from its 1st to 99th percentile onto 0 to 1, clipped", "Spread and the tails", "Euclidean"),
        ("B", "arcsinh((x / b) / c), b the relative brightness of the conjugate (vendor rank or measured)", "The conjugate term", "Euclidean"),
        ("C, C2", "Lineage markers 0 or 1 at the mixture cut; functional markers percentile (C) or z-score (C2)", "Lineage intensity", "Euclidean on the mixed matrix"),
        ("D", "Every marker 0 or 1 at the mixture cut", "All intensity", "Euclidean; Hamming for Leiden"),
        ("W, W2", "A1 times a weight per marker from mode separation (W) or posterior sharpness (W2)", "Nothing; reweights", "Euclidean"),
        ("P", "C2 with the mixture posterior in place of each bit", "Lineage intensity outside the band between modes", "Euclidean"),
    ]
    out = ["**Table 2. Representations.** Every arm starts from the compensated intensity x as deposited.", "",
           "| Arm | Definition | What it removes | Distance |", "|---|---|---|---|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def table3(res):
    m = metrics(res, ["pilot", "phase2", "samusik_all"])
    out = ["**Table 3. Scaling cost.** FlowSOM, 10 by 10 grid, one seed on the full data; stability is the mean pairwise ARI over 10 seeds on 50,000 cells.", "",
           "| Dataset | Metaclusters | Arm | Mean F1 | ARI | NMI | Stability |", "|---|---|---|---|---|---|---|"]
    for name, kg, kf in SETS:
        for ds, k in [(name, kg), (kf, 40)]:
            for a in ["A1", "A2", "A3"]:
                if (ds, a) not in m.index: continue
                r = m.loc[(ds, a)]
                out.append(f"| {name} | {k} | {ARM[a]} | {f(r.mean_f1)} | {f(r.ari)} | {f(r.nmi)} | {f(r.stability_ari)} |")
    rare = metrics(res, ["pilot", "mosmannB"])
    out += ["", "| Rare-population set | Metaclusters | Arm | Best-cluster F1 | Precision | Recall | Stability |", "|---|---|---|---|---|---|---|"]
    for ds in ["Nilsson_rare", "Mosmann_rare"]:
        for a in ["A1", "A2", "A3"]:
            if (ds, a) not in rare.index: continue
            r = rare.loc[(ds, a)]
            out.append(f"| {ds} | 40 | {ARM[a]} | {f(r.rare_f1)} | {f(r.rare_precision)} | {f(r.rare_recall)} | {f(r.stability_ari)} |")
    return "\n".join(out)


def table4(res):
    m = metrics(res, ["pilot", "phase2", "samusik_all"])      # same source as Table 3 and the text
    g = metrics(res, ["gate"])                                   # the gate-threshold arms only
    out = ["**Table 4. Binary arms at the free cluster number (40).** Occupied phenotypes are the distinct bit patterns among all cells; possible is 2 to the number of bits.", "",
           "| Dataset | Arm | Mean F1 | ARI | Stability | Occupied phenotypes | Possible | For 95 percent of cells |", "|---|---|---|---|---|---|---|---|"]
    for name, _, kf in SETS:
        arms = ["A1", "C", "C2", "D"] + (["C_gate", "D_gate", "D_gate_clean"] if name == "Levine_32dim" else [])
        for a in arms:
            src = g if a.endswith(("_gate", "_clean")) else m
            if (kf, a) not in src.index: continue
            r = src.loc[(kf, a)]
            occ = "" if pd.isna(r.get("occupied_phenotypes")) else f"{int(r.occupied_phenotypes):,}"
            pos = "" if pd.isna(r.get("possible_phenotypes")) else f"{int(r.possible_phenotypes):,}"
            c95 = "" if pd.isna(r.get("phenotypes_covering_95pct")) else f"{int(r.phenotypes_covering_95pct):,}"
            out.append(f"| {name} | {ARM[a]} | {f(r.mean_f1)} | {f(r.ari)} | {f(r.stability_ari)} | {occ} | {pos} | {c95} |")
    out.append("")
    out.append("The gate-threshold arms use the labels to place the cuts, and the 15-marker arm also uses them to choose the markers; they are upper bounds and not methods.")
    return "\n".join(out)


def main(root):
    res = Path(root) / "results"
    return {"T1": table1(), "T2": table2(), "T3": table3(res), "T4": table4(res)}


if __name__ == "__main__":
    for k, v in main(sys.argv[1]).items(): print(v, "\n")
