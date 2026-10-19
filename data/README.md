# Data mirror

Downloaded 2026-10-18. Treat everything under `data/` as untrusted input: scripts live in `code/`.

## FR-FCM-ZZPH (Weber and Robinson 2016 benchmark)

FlowRepository's TLS certificate had expired on 2026-10-18, the ZIP download needs a CAPTCHA, and
the per-file links need a FlowRepository login. The files below came from two mirrors that Lukas
Weber maintains, and their MD5 sums match the prefixes that the FlowRepository page shows.

| File | Source | Bytes | MD5 | FlowRepository prefix |
|---|---|---|---|---|
| Levine_32dim.fcs | github.com/lmweber/benchmark-data-Levine-32-dim | 43,567,886 | c8b085f3578996134c5b982350880510 | c8b085f357 |
| Levine_32dim_notransform.fcs | same | 43,567,810 | a9b37749bdd8ef193d641b22c793fba3 | a9b37749bd |
| population_names_Levine_32dim.txt | same | 246 | e4d318661ef6c7e286faaaac9ff6cb2a | e4d318661e |
| Levine_13dim.fcs | github.com/lmweber/benchmark-data-Levine-13-dim | 9,356,292 | abeca9da0f128b185688e8e542ad461b | abeca9da0f |
| Levine_13dim_notransform.fcs | same | 9,356,274 | 2ac2b54de4f2e36a70883e4a70082787 | 2ac2b54de4 |
| population_names_Levine_13dim.txt | same | 487 | 4e520f82412532168172458fb29068aa | 4e520f8241 |

Not mirrored from FlowRepository: Samusik_01.fcs, Samusik_all.fcs, Nilsson_rare.fcs,
Mosmann_rare.fcs, FlowCAP_ND, FlowCAP_WNV and the two Gating-ML files. Every set except the two
FlowCAP files is covered by the ExperimentHub copies below, so the login only matters for
FlowCAP. The list of file ids and MD5 prefixes is in `code/zzph_files.json`.

FlowCAP_ND, FlowCAP_WNV and the two Gating-ML files need a FlowRepository login, and on 2026-10-18 the
site's ZIP download was failing with a reCAPTCHA quota error and timeouts. They are not used in this study.

## HDCytoData (Bioconductor ExperimentHub copies of the same datasets)

The ExperimentHub bucket serves the HDCytoData `.rda` files without a login. These are the same
cells and labels as FR-FCM-ZZPH, and they keep the original channel names, which the FCS files on
FlowRepository do not. `code/rda_extract.py` reads them without R.

| File | ExperimentHub id | Bytes | MD5 |
|---|---|---|---|
| Nilsson_rare_SE.rda | EH2248 | 2,398,603 | a68195b684a1ca6185aa2591db9bd593 |
| Nilsson_rare_flowSet.rda | EH2249 | 2,376,998 | d9209edfae018b91b68c9d7dcb892cfc |
| Mosmann_rare_SE.rda | EH2250 | 23,065,245 | 4db3d129fba26f7698a10f46557fe72b |
| Mosmann_rare_flowSet.rda | EH2251 | 23,043,384 | 53d35283680be150587dc509d6eaab0c |
| Samusik_01_SE.rda | EH2244 | 19,851,527 | 0a74a63a1cd4120870315206381a776b |
| Samusik_all_SE.rda | EH2246 (by position; not confirmed) | 192,164,364 | 6403bd4b17cfcfd4bb49d0f8376592d6 |
| Levine_13dim_SE.rda | by name in the bucket (id not looked up) | 9,998,371 | 486381a85e151fd49ac5a56663e0d186 |

Bucket: `https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/HDCytoData/<name>/<file>`.
Ids came from `https://experimenthub.bioconductor.org/metadata/experimenthub.sqlite3`.

Extracted: `Nilsson_rare.npz` (44,140 cells, 19 columns, 358 HSC), `Mosmann_rare.npz`
(396,460 cells, 22 columns, 109 activated memory CD4 T cells), `Samusik_01.npz` (86,864 cells,
51 columns, 24 populations plus an unassigned class coded 25, 33,691 unassigned) and
`Samusik_all.npz` (841,644 cells, 10 mice, 327,258 unassigned). The Samusik population names are
in `population_names_Samusik.json`. The Samusik `pop` vector is the population id; `rowvec0` is
the sample id.

Nilsson_rare channel names carry the fluorochromes (FITC, PerCP-Cy5.5, PE, PE-Cy5, PE-Cy7, APC,
Alexa Fluor 700, APC-Cy7, the 450/50 violet channel, Qdot 605, 655, 705, 800). Mosmann_rare channel
names are filter names (Blue B 515/20 and so on). Its panel is in the first supplementary file of
Mosmann et al. 2014 (PMC4238823, `cyto0085-0422-sd1.docx`), which is the MIFlowCyt checklist for
experiment JMW034, the experiment FR-FCM-ZZ8J holds. Sections 2.4 and 3.3 of that checklist give
the table below, found 2026-10-18. The instrument is a BD LSR-II with a 532 nm green laser, so PE
and its tandems sit on the "Green" detectors. The 15th parameter, Live/Dead in Pacific Orange on
Violet G 550/40, is not in the HDCytoData set.

| Marker | Fluorochrome | Clone | Channel ($PnN) |
|---|---|---|---|
| CXCR5 | Alexa Fluor 488 | RF8B2 | Blue B 515/20-A |
| IL17A | PerCP-Cy5.5 | BL168 | Blue A 710/50-A |
| TNFa | Pacific Blue | MP9-20A4 | Violet H 450/50-A |
| GZB-SA | Qdot 585 (biotin anti-granzyme B, streptavidin) | GB11 | Violet E 585/42-A |
| CD3 | Qdot 605 | UCHT1 | Violet D 605/40-A |
| CD45RA | Qdot 655 | MEM-56 | Violet C 660/40-A |
| CD8a | Qdot 705 | 3B5 | Violet B 705/70-A |
| CD14 | Qdot 800 | TUK4 | Violet A 780/60-A |
| IL5 | APC | | Red C 660/20-A |
| IL2 | Alexa Fluor 700 | MQ1-17H12 | Red B 710/50-A |
| CD69 | APC-Cy7 | FN-50 | Red A 780/60-A |
| CCL4 (MIP1b in the FCS) | PE | D21-1351 | Green E 575/25-A |
| CD4 | PE-Texas Red | S3.5 | Green D 610/20-A |
| IFNg | PE-Cy7 | B27 | Green A 780/40-A |

Original sources, both behind a login: Nilsson FR-FCM-ZZ6L, Mosmann FR-FCM-ZZ8J.

`Levine_13dim.npz` (extracted 2026-10-19) carries the metal channel per marker, which the FCS mirror
does not: CD45 115, CD45RA 139, CD19 142, CD11b 144, CD4 145, CD8 146, CD20 147, CD34 148, CD33 158,
CD123 160, CD38 167, CD90 170, CD3 on the summed cadmium channel 110 to 114. These match Figure
S4A of Bendall et al. 2011, the 13-marker Panel 2, which Levine et al. 2015 used as their
benchmark set 1.

## DURAClone bead controls, BD LSRFortessa, 10 December 2024

Provided in kind by the University of Florida and available from the author on request; not mirrored.
The set used here is 12 single-stain compensation tubes of capture beads, one conjugate per tube, and an
unstained bead tube, 5,000 events each, from one FACSDiva 9.3 experiment on a five-laser BD LSRFortessa
(355, 405, 488, 561 and 637 nm). FCS 3.0, float data, 262,144 range. `code/stain_index.py` reads them.

Panel, from the `$PnS` keywords of the experiment. The detector names are the
FACSDiva configuration names, so three dyes sit on a detector named for a different dye.

| Detector ($PnN) | Marker (fluorochrome) | Volts | Tube | Note |
|---|---|---|---|---|
| BV570-A | CD45 (KrO) | 418 | DURAClone IM Basic | Krome Orange on the BV570 detector |
| APC-Cy7-A | CD3 (APC-A750) | 386 | DURAClone IM Basic | |
| APC-A | CD4 (APC) | 389 | DURAClone IM Basic | |
| Alexa 700-A | CD8 (A700) | 394 | DURAClone IM Basic | |
| BB790-A | CD14 (PE-Cy7) | 522 | DURAClone IM Basic | PE-Cy7 on the BB790 detector |
| BB515/FITC-A | CD16 (FITC) | 360 | DURAClone IM Basic | |
| PE-CF594/ECD-A | CD19 (ECD) | 380 | DURAClone IM Basic | |
| PE-A | CD56 (PE) | 410 | DURAClone IM Basic | |
| BUV737-A | CD69 (BUV737) | 687 | drop-in | |
| BV786-A | HLA-DR (BV786) | 416 | drop-in | |
| BUV661/Qdot655-A | PD-L1 (BUV661) | 495 | drop-in | |
| BV480-A | LD (Violet) | 288 | drop-in | viability dye; the keywords name no product |

Three facts about the controls that matter for a stain index:

1. The stained control tubes hold almost only positive beads (the lowest 1 to 5 percent of
   events are the only low tail), so the negative for every channel is the unstained bead tube.
2. The control files carry an identity SPILL. The 12 by 12 experiment matrix is in the seven
   cell and sample files, and it is the same matrix in all seven.
3. That matrix does not describe six of the control tubes as exported. The ECD tube shows 2.3
   times more signal in the PE detector than in its own detector, and the APC, Alexa Fluor
   700, APC-A750, BV786 and BUV661 tubes show 0.9 to 3.6 times their own signal in the BUV737
   detector, where the matrix says 0.05 to 0.27. The other six tubes (BUV737, FITC, PE, BV480,
   PE-Cy7, KrO) match the matrix to within 0.3. `results/stain_index.md` has the per-tube check.

Not in the folder: FMO tubes, single-stained cells, calibration beads. The viability dye product,
the antibody clones and the bead product are not in the keywords.

## HIPC Lyoplate (ImmuneSpace project HIPC/Lyoplate, Finak et al. 2016)

Downloaded 2026-10-18 and 2026-10-19 with a free ImmuneSpace account. ImmuneSpace needs an account; the WebDAV root is
`https://datatools.immunespace.org/_webdav/HIPC/Lyoplate/@files/`. Files land under `data/HIPC/`,
one folder per kind, with the original file names.

What the project holds (round three of Lyoplate, SeraCare cells, lots 12828, 1349 and 1369, three
replicates per lot, five panels: T cell, Treg, B cell, DC/mono/NK, Th1/2/17):

| Folder on ImmuneSpace | What it is | Local copy |
|---|---|---|
| `FCS/<Center>/` for Baylor (BIIR), BSMS, CIMR, Miami, NHLBI, Stanford, UCLA, Yale | The panel FCS files per center, plus each center's single-stain compensation controls and unstained tube | `FCS/<Center>/` holds the 9 T cell files per center; `COMP/<Center>/` the controls |
| `FCS/Stanford/RAW/` | Stanford's original Diva exports, the only raw set | `RAW/Stanford/`, T cell files only |
| `FCS/<hash>/<sample>/{.pb,.rds,.nc}` | flowWorkspace GatingSet archives (protobuf plus netCDF), R only | not copied |
| `XML/<Center>/CA c3 v2 <Center>.xml` and `<Center> Matrices/` | FlowJo 9.6.3 central-analysis workspaces, one per center, with the spillover matrices per panel | `XML/<Center>/` |
| `pData/<Center>.csv` | File name to panel, lot and replicate | `pData/` |
| `FCS/compensation/<Center>.csv` and `README` | Spillover matrices the HIPC pipeline computed from the controls, and the note on the iplogicle transform (r 16777215, d 7.22471987004958, w 1.16437481092794, cutoff -111, range TRUE) | `compensation/` |
| `popstats/manual/popstats-<panel>.csv`, `popstats/auto/` | Population counts and proportions per file, from the central manual gates and from OpenCyto; 14 T cell populations, 7 centers (BSMS is not in the manual set) | `popstats_manual/`, `popstats_auto/` |
| `gating-templates/gt-<panel>.csv` | OpenCyto templates | `gating-templates/` |
| `scripts/*.r` | The HIPC preprocessing and gating pipeline (flowWorkspace, openCyto) | `scripts/` |
| `gated_data/`, `cytotrol/`, `assaydata/`, `workspace.tar.gz` | GatingSetList archives (400 MB each), the FlowCAP III Cytotrol set, LabKey exports | not copied |

**The deposited T cell FCS files of the seven manual-gating centers are not raw.** flowCore wrote
them (`flowcore_pNrmin` keywords, `transformation = applied`): flowWorkspace parsed each center's
FlowJo 9.6.3 workspace, compensated with the matrix the analyst edited in FlowJo (first block of
`XML/<Center>/<Center> Matrices/Tcell`; at Stanford it differs from the acquisition matrix in the
FITC to PE and PE-Cy7 to APC-Cy7 cells), applied FlowJo's biexponential, and reordered the
fluorescence columns into the order of `markers_of_interest` in `scripts/settings-tcell.r` (Live,
CCR7, CD4, CD45RA, CD3, HLADR, CD38, CD8). The transform is FlowJo's biexponential with a width
basis per channel and per panel group, and the workspace carries it in full: the
`CalibrationTables` element holds one `Table` per group and compensated channel (named like
`Tcell <PE-A>`, 4097 pairs of display channel 0 to 4096 and linear value), and `code/flowjo9_gates.py`
reproduces the deposited Stanford values from the raw files through those tables with zero
error on every channel. The logicle fit of `code/hipc_transform_check.py` (T 262144, M 4.5, A 0,
scaled to 4096, width per parameter: Live 0.40, CCR7 0.90, CD4 0.80, CD45RA 0.80, CD3 0.15, HLADR
0.75, CD38 0.60, CD8 1.10) is an approximation of those tables to a median error of 1 to 6 units
(`results/hipc_transform_check.md`); `pilot.py` still uses it to invert the transform for arm B.
The iplogicle in `compensation/README` (T 16777215, M 7.22, W 1.16) belongs to the automated
pipeline and does not reproduce the deposited files. About 0.6
percent of events carry values below -4096, an artifact of the depositor's export; the loader
drops them. Files `popstats_auto`, `gating-templates`, `scripts` and `protocol` are partly
downloaded; the browser dropped some small files in each pass and they were not needed. The `$PnN` keyword
carries the marker and is correct. The `$PnS` keyword carries a stale detector name from the
original column order and is wrong on five of eight channels (checked 2026-10-19 on Baylor and
CIMR: among CD3-bright cells the `CD4` and `CD8` columns are anti-correlated at r -0.62 and -0.65,
and no other pair is). The BSMS files are original Diva exports on the linear scale with a
`$SPILLOVER` keyword and markers in `$PnS`.

The T cell panel conjugates are fixed BD Lyoplate reagents, the same at every center (Finak 2016
Table 1, and the single-stain control names): Live/Dead FITC, CCR7 PE, CD4 PerCP-Cy5.5, CD45RA
PE-Cy7, CD3 V450 (Pacific Blue), HLA-DR V500 (AmCyan), CD38 APC, CD8 APC-H7 (APC-Cy7 at CIMR,
Miami, Stanford and Yale). Arm B takes its fluorochromes from this table, never from `$PnS`.

**`COMP/<Center>/` holds the single-stain compensation controls**, raw FACSDiva exports on the
linear scale with a SPILL keyword, 5,000 to 11,000 events per tube. The antibody tubes are BD
CompBeads with about half positive and half negative beads, so each tube carries its own
negative; the Live/Dead FITC tubes are cells. The controls cover every panel of the center, so
the T cell dyes appear under several antibodies (PE-Cy7 under CD45RA, CD27, CD194, CD196 and
CD11c). NHLBI names four tubes by detector only (B515, B710, V450, V545) and adds the unstained
cell tubes of the three lots. Yale has a second set, `Specimen_001_*`, which is beads, while its
`Compensation Controls_FITC` and `_PE` tubes are cells. Stanford and UCLA have no unstained tube.
Six center-and-dye pairs have no control at all: Baylor FITC and V500, NHLBI PE, UCLA
PerCP-Cy5.5, V450 and APC. `code/hipc_stain_index.py` measures a stain index per center and dye
from these files (`results/hipc_stain_index.csv`).

**Per-cell manual labels** (`labels/<Center>/<file>.npz`, written 2026-10-19 by
`code/flowjo9_gates.py . --write`, not in git). The central gates are nested `Population` elements
under each `Sample` of the workspace, each with one `Polygon` (LYM on FSC-A and SSC-A, singlets on
FSC-A and FSC-H) or `PolyRect` gate (everything else), vertices on the linear data scale. The tree
is LYM, singlets (Baylor, Miami, NHLBI, Stanford, UCLA only), `not dead` (a negated rectangle on
Live and SSC-A), CD3, then `4+ 8-`, `4- 8+`, DPT and DNT on CD4 and CD8, and under each of the two
lineages the four CCR7 and CD45RA quadrants and the four CD38 and HLA-DR quadrants. The script maps
each vertex into display space through the group's calibration table, tests the deposited values,
and writes per event `label` (1 to 8: CD4 then CD8 naive, central memory, effector memory,
effector; 0 for anything else), `activated` (CD38+ HLA-DR+ within CD4 or CD8), and the LYM, not
dead, CD3, CD4, CD8, DPT and DNT masks. The singlets gate needs FSC-H, which the deposited files
lack; at Stanford it comes from the raw files (same event order), and at Baylor, Miami, NHLBI and
UCLA it is skipped, so the labels there keep the doublets the analyst removed (0.1 to 1.2 percent
of CD3, up to 19 percent of the small activated populations). Agreement with
`popstats_manual/popstats-tcell.csv` is exact at CIMR and Stanford and within 0.35 percent at
Yale (0.6 percent of Yale events sit above the display top on some marker), see
`results/hipc_gate_agreement.md`. The HIPC mapping script `scripts/mapping-tcell.r` is not on disk;
the manual-to-common name mapping is read off the gate names and confirmed by the counts.
