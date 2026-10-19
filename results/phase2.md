# Phase 2 — the remaining mass cytometry sets and arm B on Mosmann_rare

**Run:** 2026-10-18 · **Code:** `code/pilot.py --datasets Levine_13dim,Levine_13dim_k40,Samusik_01,Samusik_01_k40 --tag phase2` and `--datasets Mosmann_rare --arms A1,A2,A3,B_geom,B_lin --tag mosmannB`, tables by `code/phase2_report.py` · **Metrics:** `results/phase2_metrics.csv`, `results/phase2_perpop.csv`, `results/mosmannB_metrics.csv`

## The question

The pilot settled hypothesis 2 on Levine_32dim and left a caveat: Levine_13dim, Samusik and
the conventional flow sets had not run. This run asks whether the three results of the pilot
hold on the two mass cytometry sets that remain, and it runs arm B on Mosmann_rare now that
the panel is known.

The three results hold. Binary input loses to arcsinh intensity at the free cluster number on
both sets. Per-channel scaling costs mean F1 at the free cluster number on both sets. The
populations the binary arms lose most are the ones a marker defines at three levels.

## What ran

| Dataset | Cells | Labeled | Markers | Lineage bits | Metaclusters | Ground truth |
|---|---|---|---|---|---|---|
| Levine_13dim | 167,044 | 81,747 | 13 | 12 (CD38 continuous) | 24 (the gated number) and 40 (free) | 24 populations, human bone marrow, 1 donor, mass cytometry |
| Samusik_01 | 86,864 | 53,173 | 39 | 28 (11 activation and maturation markers continuous) | 24 and 40 | 24 populations, mouse bone marrow, 1 of 10 mice, mass cytometry. Population 25, unassigned, scores as unlabeled. |

Same arms, same FlowSOM settings, same seeds and same scoring as the pilot. Arm B does not run on
mass cytometry. The lineage and continuous split for Samusik_01 is a judgment made in this
session: CD44, CD25, CD27, Ly6C, Sca1, CD150, CD43, CD103, CCR7, CD16_32 and Foxp3 stay
continuous because the label set uses them at more than one level.

## Hypothesis 1: arm B on Mosmann_rare, the second set with known fluorochromes

The panel came from the JMW034 MIFlowCyt checklist in the supplement of Mosmann et al. 2014
(`data/README.md` has the table). Three ranks come from the fluor-ai table (Alexa Fluor 488
at 3, PE-Texas Red at 4) and Qdot 585 is assumed at 3 with the other dim Qdots. The 2^rank
divisors, geometric mean 1, are 0.64 for the nine rank-3 dyes, 1.28 for the three Qdot 605,
655 and PE-Texas Red channels (CD3, CD45RA, CD4) and 2.56 for APC, PE and PE-Cy7 (IL5, CCL4,
IFNg).

| Arm | Rare F1 | Precision | Recall | Stability ARI | ARI against A1 | ARI against A2 |
|---|---|---|---|---|---|---|
| A1 arcsinh | 0.57 | 0.40 | 0.98 | 0.52 |  | 0.58 |
| A2 z-score | 0.61 | 0.44 | 0.99 | 0.57 | 0.58 |  |
| A3 percentile | 0.09 | 0.05 | 0.72 | 0.59 | 0.52 | 0.55 |
| B brightness, 2^rank | 0.13 | 0.07 | 0.94 | 0.49 | 0.56 | 0.47 |
| B brightness, linear rank | 0.13 | 0.07 | 0.95 | 0.57 | 0.78 | 0.60 |

Membership changes again (ARI 0.56 and 0.78 against A1), and arm B agrees with A1 better than
with the z-score control A2 (0.47 and 0.60), so on the second dataset as on the first, dividing
by brightness is not per-channel scaling. The hypothesis 1 statement now rests on two
conventional flow sets.

Arm B also costs the rare population on both sets. Nilsson_rare went from 0.50 to 0.22 in the
pilot. Mosmann_rare goes from 0.57 to 0.13 here, with recall intact (0.94) and precision gone
(0.07): the activated memory CD4 T cells stay together and 14 times as many other cells join
them. The three cytokines that define the rare population, IL5, CCL4 and IFNg, are on APC, PE
and PE-Cy7, the three brightest dyes in the panel, and arm B divides them by 2.56 while it
multiplies the nine dim channels by 1.56. The panel designer put the bright dyes on the
markers the question depends on, and dividing by brightness removes that weight. The same
thing happened on Nilsson_rare, where CD34 and CD38, the two markers that define the HSC, are on
PerCP-Cy5.5 and FITC and got multiplied up with their noise. The paper has to say this
plainly. Brightness normalization makes the representation independent of the panel, and the
panel was designed to make the representation good. Which of those a reader wants depends on
whether the question is the one the panel was built for.

Two things this does not settle. The divisors are vendor ranks on another instrument, with
five of 14 assumed, and a measured stain index from single stains can move them a lot. And
both datasets score one rare population, so a cell that moves cluster anywhere else in the
file costs nothing here and the F1 column understates what changed. Levine-style full labels
on a conventional flow set with known fluorochromes are the data need for hypothesis 1, and no
public set in tier 1 has them. The HIPC Lyoplate set (tier 2) does.

## Hypothesis 2: binary input against arcsinh intensity

Levine_13dim, 40 metaclusters. Mean F1 0.53 for A1 against 0.25 to 0.32 for the binary arms, and
the binary arms agree with A1 at ARI 0.26 to 0.29. Samusik_01, 40 metaclusters. Mean F1 0.70
for A1 against 0.38 to 0.48, agreement 0.26 to 0.34. On both sets the binary arms are also less
stable across seeds. The pattern is the Levine_32dim pattern: the mixture threshold on a mass
cytometry channel cuts inside the zero spike, the bits are noisy, and the phenotype count says
so. Samusik_01 has 28 lineage bits and the cells occupy 65,032 phenotypes among 86,864 cells.
Levine_13dim has 12 lineage bits, 4,096 possible phenotypes, and the cells occupy 3,226 of
them with 1,480 needed for 95 percent of cells. Nilsson_rare, 12 conventional flow bits, needed
200. Mass cytometry channels binarize worse than fluorescence channels under the same mixture
rule, and that is a result in its own right for the Methods.

At the gated cluster number the order is less clean. On Levine_13dim at 24, A1 merges the naive
and mature CD4 T cells (F1 0.02 and 0.80) and the mature CD8 T cells into the naive ones (0.00),
and A3 keeps them apart (0.96, 0.97, but 0.00 for mature CD8). Mean F1 at 24 is 0.39 for A3,
0.36 for A1, 0.30 for C2. On Samusik_01 at 24, A1 and A3 tie at 0.62 and A1 merges the
intermediate monocytes into the classical ones (0.00). These are merge-order effects at a
forced cluster number, the same thing the pilot saw with CD8 T cells on Levine_32dim at 14, and
they reverse at 40.

### Levine_13dim, FlowSOM, 24 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.36 | 0.67 | 0.80 | 0.83 |  |  |  |  |
| A2 z-score | 0.33 | 0.65 | 0.78 | 0.82 | 0.85 |  |  |  |
| A3 percentile | 0.39 | 0.74 | 0.82 | 0.81 | 0.81 |  |  |  |
| C2 lineage mixture + z-score | 0.30 | 0.47 | 0.64 | 0.53 | 0.39 | 3,176 | 4,096 | 1,431 |
| D all mixture | 0.29 | 0.49 | 0.63 | 0.44 | 0.39 | 5,461 | 8,192 | 2,535 |
| C lineage mixture + percentile | 0.20 | 0.26 | 0.49 | 0.40 | 0.28 | 3,176 | 4,096 | 1,431 |

### Levine_13dim, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.53 | 0.79 | 0.83 | 0.70 |  |  |  |  |
| A2 z-score | 0.41 | 0.59 | 0.75 | 0.65 | 0.53 |  |  |  |
| A3 percentile | 0.48 | 0.61 | 0.77 | 0.69 | 0.54 |  |  |  |
| C2 lineage mixture + z-score | 0.32 | 0.37 | 0.61 | 0.45 | 0.29 | 3,226 | 4,096 | 1,480 |
| D all mixture | 0.29 | 0.32 | 0.58 | 0.40 | 0.26 | 5,562 | 8,192 | 2,640 |
| C lineage mixture + percentile | 0.25 | 0.30 | 0.55 | 0.39 | 0.26 | 3,226 | 4,096 | 1,480 |

### Samusik_01, FlowSOM, 24 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.62 | 0.60 | 0.78 | 0.73 |  |  |  |  |
| A2 z-score | 0.54 | 0.56 | 0.73 | 0.78 | 0.81 |  |  |  |
| A3 percentile | 0.62 | 0.61 | 0.78 | 0.75 | 0.75 |  |  |  |
| C2 lineage mixture + z-score | 0.29 | 0.36 | 0.54 | 0.66 | 0.49 | 65,172 | 268,435,456 | 60,829 |
| D all mixture | 0.47 | 0.47 | 0.63 | 0.58 | 0.32 | 82,338 | 549,755,813,888 | 77,995 |
| C lineage mixture + percentile | 0.42 | 0.41 | 0.60 | 0.58 | 0.26 | 65,172 | 268,435,456 | 60,829 |

### Samusik_01, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.70 | 0.85 | 0.83 | 0.75 |  |  |  |  |
| A2 z-score | 0.60 | 0.48 | 0.66 | 0.71 | 0.46 |  |  |  |
| A3 percentile | 0.66 | 0.60 | 0.76 | 0.67 | 0.59 |  |  |  |
| C2 lineage mixture + z-score | 0.38 | 0.38 | 0.54 | 0.55 | 0.34 | 65,032 | 268,435,456 | 60,689 |
| D all mixture | 0.48 | 0.29 | 0.58 | 0.44 | 0.26 | 82,394 | 549,755,813,888 | 78,051 |
| C lineage mixture + percentile | 0.48 | 0.31 | 0.56 | 0.45 | 0.28 | 65,032 | 268,435,456 | 60,689 |

### Per-population F1, Levine_13dim, 24 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D | C |
|---|---|---|---|---|---|---|---|
| Mature CD4+ T cells | 13,964 | 0.80 | 0.80 | 0.97 | 0.58 | 0.71 | 0.30 |
| Erythroblast cells | 12,030 | 0.80 | 0.72 | 0.80 | 0.61 | 0.62 | 0.36 |
| Naive CD8+ T cells | 9,564 | 0.71 | 0.70 | 0.70 | 0.44 | 0.57 | 0.35 |
| Mature CD8+ T cells | 7,821 | 0.00 | 0.01 | 0.00 | 0.48 | 0.67 | 0.31 |
| Mature CD38lo B cells | 7,796 | 0.97 | 0.96 | 0.96 | 0.93 | 0.73 | 0.57 |
| Naive CD4+ T cells | 6,987 | 0.02 | 0.00 | 0.96 | 0.87 | 0.64 | 0.49 |
| CD11bhi Monocyte cells | 6,779 | 0.94 | 0.95 | 0.95 | 0.68 | 0.84 | 0.80 |
| NK cells | 3,864 | 0.91 | 0.92 | 0.77 | 0.69 | 0.74 | 0.53 |
| Megakaryocyte cells | 3,684 | 0.00 | 0.10 | 0.00 | 0.31 | 0.11 | 0.31 |
| Myelocyte cells | 3,025 | 0.37 | 0.26 | 0.38 | 0.21 | 0.18 | 0.14 |
| CD11bmid Monocyte cells | 1,278 | 0.00 | 0.00 | 0.03 | 0.18 | 0.10 | 0.13 |
| Pre-B II cells | 994 | 0.66 | 0.67 | 0.67 | 0.51 | 0.49 | 0.06 |
| CD11b- Monocyte cells | 912 | 0.67 | 0.68 | 0.71 | 0.30 | 0.26 | 0.18 |
| Mature CD38mid B cells | 608 | 0.57 | 0.01 | 0.01 | 0.07 | 0.03 | 0.06 |
| Immature B cells | 502 | 0.00 | 0.01 | 0.02 | 0.03 | 0.02 | 0.00 |
| Plasma cell cells | 468 | 0.00 | 0.06 | 0.00 | 0.00 | 0.16 | 0.00 |
| Plasmacytoid DC cells | 293 | 0.65 | 0.62 | 0.62 | 0.18 | 0.03 | 0.02 |
| HSC cells | 261 | 0.00 | 0.00 | 0.24 | 0.06 | 0.08 | 0.00 |
| CMP cells | 253 | 0.00 | 0.00 | 0.00 | 0.01 | 0.00 | 0.00 |
| Pre-B I cells | 240 | 0.50 | 0.43 | 0.43 | 0.00 | 0.00 | 0.00 |
| MEP cells | 194 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.10 |
| MPP cells | 152 | 0.04 | 0.04 | 0.12 | 0.02 | 0.00 | 0.00 |
| GMP cells | 73 | 0.00 | 0.02 | 0.02 | 0.05 | 0.00 | 0.01 |
| Platelet cells | 5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

### Per-population F1, Levine_13dim, 40 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D | C |
|---|---|---|---|---|---|---|---|
| Mature CD4+ T cells | 13,964 | 0.94 | 0.81 | 0.77 | 0.52 | 0.42 | 0.50 |
| Erythroblast cells | 12,030 | 0.78 | 0.59 | 0.64 | 0.44 | 0.30 | 0.32 |
| Naive CD8+ T cells | 9,564 | 0.98 | 0.83 | 0.78 | 0.57 | 0.40 | 0.37 |
| Mature CD8+ T cells | 7,821 | 0.96 | 0.83 | 0.82 | 0.81 | 0.60 | 0.40 |
| Mature CD38lo B cells | 7,796 | 0.97 | 0.71 | 0.77 | 0.81 | 0.82 | 0.63 |
| Naive CD4+ T cells | 6,987 | 0.90 | 0.93 | 0.96 | 0.45 | 0.51 | 0.63 |
| CD11bhi Monocyte cells | 6,779 | 0.70 | 0.72 | 0.69 | 0.48 | 0.62 | 0.78 |
| NK cells | 3,864 | 0.81 | 0.80 | 0.79 | 0.59 | 0.74 | 0.66 |
| Megakaryocyte cells | 3,684 | 0.65 | 0.23 | 0.23 | 0.27 | 0.29 | 0.27 |
| Myelocyte cells | 3,025 | 0.27 | 0.23 | 0.29 | 0.15 | 0.16 | 0.20 |
| CD11bmid Monocyte cells | 1,278 | 0.32 | 0.19 | 0.14 | 0.16 | 0.13 | 0.12 |
| Pre-B II cells | 994 | 0.72 | 0.67 | 0.68 | 0.56 | 0.48 | 0.22 |
| CD11b- Monocyte cells | 912 | 0.57 | 0.64 | 0.71 | 0.26 | 0.37 | 0.20 |
| Mature CD38mid B cells | 608 | 0.58 | 0.21 | 0.46 | 0.14 | 0.26 | 0.08 |
| Immature B cells | 502 | 0.04 | 0.06 | 0.16 | 0.11 | 0.08 | 0.08 |
| Plasma cell cells | 468 | 0.63 | 0.08 | 0.15 | 0.77 | 0.23 | 0.01 |
| Plasmacytoid DC cells | 293 | 0.65 | 0.62 | 0.62 | 0.06 | 0.07 | 0.05 |
| HSC cells | 261 | 0.12 | 0.18 | 0.14 | 0.09 | 0.07 | 0.10 |
| CMP cells | 253 | 0.00 | 0.02 | 0.12 | 0.21 | 0.06 | 0.04 |
| Pre-B I cells | 240 | 0.50 | 0.43 | 0.58 | 0.00 | 0.09 | 0.01 |
| MEP cells | 194 | 0.38 | 0.00 | 0.51 | 0.10 | 0.20 | 0.13 |
| MPP cells | 152 | 0.04 | 0.04 | 0.18 | 0.02 | 0.02 | 0.01 |
| GMP cells | 73 | 0.23 | 0.03 | 0.20 | 0.06 | 0.14 | 0.09 |
| Platelet cells | 5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

### Per-population F1, Samusik_01, 24 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D | C |
|---|---|---|---|---|---|---|---|
| Classical Monocytes | 13,607 | 0.68 | 0.67 | 0.76 | 0.57 | 0.66 | 0.47 |
| Intermediate Monocytes | 12,045 | 0.00 | 0.04 | 0.50 | 0.14 | 0.42 | 0.40 |
| IgDpos IgMpos B cells | 8,768 | 0.91 | 0.92 | 0.92 | 0.70 | 0.93 | 0.94 |
| Eosinophils | 5,173 | 0.89 | 0.85 | 0.89 | 0.16 | 0.60 | 0.52 |
| IgD- IgMpos B cells | 3,508 | 0.28 | 0.60 | 0.63 | 0.09 | 0.72 | 0.72 |
| pDCs | 2,974 | 0.96 | 0.69 | 0.98 | 0.94 | 0.95 | 0.98 |
| IgM- IgD- B-cells | 1,123 | 0.56 | 0.00 | 0.01 | 0.19 | 0.01 | 0.00 |
| Non-Classical Monocytes | 1,016 | 1.00 | 0.85 | 0.98 | 0.92 | 0.86 | 0.81 |
| NKT cells | 961 | 0.90 | 0.91 | 0.93 | 0.62 | 0.64 | 0.59 |
| CD8 T cells | 798 | 0.88 | 0.96 | 0.98 | 0.03 | 0.00 | 0.00 |
| GMP | 703 | 0.41 | 0.63 | 0.74 | 0.35 | 0.46 | 0.58 |
| MEP | 420 | 0.78 | 0.45 | 0.82 | 0.74 | 0.80 | 0.74 |
| CD4 T cells | 406 | 0.93 | 0.93 | 0.93 | 0.32 | 0.82 | 0.00 |
| Basophils | 292 | 0.87 | 0.82 | 0.89 | 0.01 | 0.72 | 0.71 |
| B-cell Frac A-C (pro-B cells) | 266 | 0.05 | 0.01 | 0.01 | 0.01 | 0.03 | 0.00 |
| Macrophages | 199 | 0.88 | 0.69 | 0.90 | 0.02 | 0.46 | 0.40 |
| CMP | 197 | 0.46 | 0.00 | 0.00 | 0.05 | 0.01 | 0.02 |
| MPP | 197 | 0.85 | 0.00 | 0.48 | 0.42 | 0.47 | 0.43 |
| NK cells | 189 | 0.66 | 0.56 | 0.69 | 0.25 | 0.63 | 0.59 |
| Plasma Cells | 118 | 0.99 | 0.97 | 0.96 | 0.52 | 0.49 | 0.59 |
| mDCs | 116 | 0.78 | 0.93 | 0.85 | 0.00 | 0.49 | 0.36 |
| CLP | 59 | 0.19 | 0.52 | 0.14 | 0.01 | 0.19 | 0.20 |
| gd T cells | 35 | 0.08 | 0.02 | 0.00 | 0.02 | 0.00 | 0.00 |
| HSC | 3 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

### Per-population F1, Samusik_01, 40 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D | C |
|---|---|---|---|---|---|---|---|
| Classical Monocytes | 13,607 | 0.96 | 0.73 | 0.73 | 0.58 | 0.41 | 0.37 |
| Intermediate Monocytes | 12,045 | 0.91 | 0.45 | 0.61 | 0.21 | 0.28 | 0.29 |
| IgDpos IgMpos B cells | 8,768 | 0.98 | 0.75 | 0.98 | 0.83 | 0.72 | 0.85 |
| Eosinophils | 5,173 | 0.89 | 0.85 | 0.89 | 0.51 | 0.40 | 0.50 |
| IgD- IgMpos B cells | 3,508 | 0.56 | 0.59 | 0.62 | 0.62 | 0.72 | 0.59 |
| pDCs | 2,974 | 0.53 | 0.64 | 0.65 | 0.74 | 0.86 | 0.95 |
| IgM- IgD- B-cells | 1,123 | 0.43 | 0.28 | 0.31 | 0.21 | 0.01 | 0.43 |
| Non-Classical Monocytes | 1,016 | 1.00 | 0.85 | 0.98 | 0.87 | 0.86 | 0.75 |
| NKT cells | 961 | 0.90 | 0.79 | 0.93 | 0.81 | 0.64 | 0.64 |
| CD8 T cells | 798 | 0.88 | 0.96 | 0.98 | 0.74 | 0.00 | 0.00 |
| GMP | 703 | 0.78 | 0.61 | 0.80 | 0.32 | 0.58 | 0.57 |
| MEP | 420 | 0.78 | 0.81 | 0.78 | 0.67 | 0.81 | 0.75 |
| CD4 T cells | 406 | 0.93 | 0.93 | 0.93 | 0.33 | 0.87 | 0.86 |
| Basophils | 292 | 0.87 | 0.82 | 0.89 | 0.02 | 0.78 | 0.69 |
| B-cell Frac A-C (pro-B cells) | 266 | 0.12 | 0.01 | 0.04 | 0.04 | 0.17 | 0.01 |
| Macrophages | 199 | 0.88 | 0.69 | 0.79 | 0.02 | 0.62 | 0.56 |
| CMP | 197 | 0.44 | 0.25 | 0.25 | 0.24 | 0.30 | 0.19 |
| MPP | 197 | 0.85 | 0.00 | 0.70 | 0.43 | 0.55 | 0.49 |
| NK cells | 189 | 0.66 | 0.65 | 0.58 | 0.21 | 0.62 | 0.61 |
| Plasma Cells | 118 | 0.99 | 0.97 | 0.96 | 0.32 | 0.38 | 0.68 |
| mDCs | 116 | 0.78 | 0.93 | 0.85 | 0.45 | 0.56 | 0.47 |
| CLP | 59 | 0.54 | 0.52 | 0.48 | 0.02 | 0.32 | 0.18 |
| gd T cells | 35 | 0.08 | 0.19 | 0.04 | 0.00 | 0.00 | 0.01 |
| HSC | 3 | 0.00 | 0.00 | 0.00 | 0.00 | 0.06 | 0.00 |

## Hypothesis 3: the three-level markers

Levine_13dim labels monocytes at three CD11b levels and mature B cells at two CD38 levels.
Samusik_01 labels monocytes at three Ly6C levels (classical, intermediate, non-classical). Per-
population F1 at 40 metaclusters:

| Population | Cells | A1 | A2 | A3 | C2 | D |
|---|---|---|---|---|---|---|
| CD11bhi monocytes (Levine_13dim) | 6,779 | 0.70 | 0.72 | 0.69 | 0.48 | 0.62 |
| CD11bmid monocytes | 1,278 | 0.32 | 0.19 | 0.14 | 0.16 | 0.13 |
| CD11b- monocytes | 912 | 0.57 | 0.64 | 0.71 | 0.26 | 0.37 |
| Mature CD38lo B cells | 7,796 | 0.97 | 0.71 | 0.77 | 0.81 | 0.82 |
| Mature CD38mid B cells | 608 | 0.58 | 0.21 | 0.46 | 0.14 | 0.26 |
| Classical monocytes (Samusik_01) | 13,607 | 0.96 | 0.73 | 0.73 | 0.58 | 0.41 |
| Intermediate monocytes | 12,045 | 0.91 | 0.45 | 0.61 | 0.21 | 0.28 |
| Non-classical monocytes | 1,016 | 1.00 | 0.85 | 0.98 | 0.87 | 0.86 |

The middle level is the one that goes. CD11bmid monocytes fall from 0.32 to 0.13 to 0.16, and
the intermediate monocytes on Samusik_01 fall from 0.91 to 0.21 to 0.28, while the high and
the negative levels hold up better. A single cut on a three-level marker puts the middle level
with one of its neighbors, and the cell count of that level then decides whether the merged
cluster matches it or the neighbor. The 12,045 intermediate monocytes on Samusik_01 are the
clearest case, because they are the second-largest population in the set and no binary arm
reaches 0.3 on them. CD38mid B cells (0.58 to 0.14) are the same mechanism on the one marker
that the Levine_13dim arms keep continuous, which says the loss there comes from the lineage
bits around CD38 and not from CD38 itself.

## The scaling result

Z-scoring costs mean F1 at the free cluster number on every mass cytometry set so far:

| Dataset, 40 metaclusters | A1 | A2 | Cost | A3 | ARI A1 against A2 |
|---|---|---|---|---|---|
| Levine_32dim | 0.79 | 0.63 | 0.16 | 0.64 | 0.67 |
| Levine_13dim | 0.53 | 0.41 | 0.12 | 0.48 | 0.53 |
| Samusik_01 | 0.70 | 0.60 | 0.10 | 0.66 | 0.46 |

The agreement between A1 and A2 is 0.46 to 0.67, so the z-scored clustering is a different
partition of the same cells and not a small perturbation. On Samusik_01 the ARI against the
labels falls from 0.85 to 0.48. At the gated cluster number the cost is smaller or absent (0.03
and 0.08), and A3 beats A1 on Levine_13dim there, so the statement for the paper is about the
free cluster number, which is how FlowSOM is used in practice.

## Caveats

- One seed for every full-data result, as in the pilot.
- The stability numbers for a dataset move between runs, because the 50,000-cell subsample
  comes from the module's random state and that state depends on which jobs ran before. A1 on
  Mosmann_rare scored 0.56 in the pilot and 0.52 in the mosmannB run with the same clustering
  (rare F1 0.571 in both). Read stability to one decimal.
- No gate thresholds for Levine_13dim or Samusik_01. The Levine_32dim run showed gate
  thresholds do not change the verdict, and the phenotype counts here say the same noise is
  present, so this run did not repeat that test.
- The Samusik lineage split is a judgment. Moving Ly6C into the lineage set makes the
  intermediate monocyte result worse for the binary arms, not better.
- Samusik_all, the other 9 mice, has not run. Samusik_01 is one mouse.
