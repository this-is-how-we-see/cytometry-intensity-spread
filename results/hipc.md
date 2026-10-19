# HIPC Lyoplate T cell panel — between-site variance per arm (phase 3)

**Run:** 2026-10-19 · **Code:** `code/pilot.py --datasets HIPC_tcell,HIPC_tcell_k30 --arms A1,A2,A3,B_geom,B_lin,C2,D --tag hipc`, then `code/between_site.py`, tables by `code/phase2_report.py --tag hipc`; transform work in `code/hipc_transform_check.py` and `code/hipc_width_check.py` · **Metrics:** `results/hipc_metrics.csv`, `results/hipc_between_site.csv`, `results/hipc_between_site_clusters.csv`, `results/hipc_between_site_match.csv`, `results/hipc_transform_check.md`, `results/hipc_width_check.md`; the measured-divisor pass adds `results/hipc_stain_index.csv`, `results/hipc_site_metrics.csv` and `results/hipc_site_k12.csv`, `results/hipc_site_k30.csv` · **Data:** `data/README.md`, HIPC section

## The question

The plan reserves one measurement for this dataset: the same samples ran at seven centers
with the same lyophilized BD reagents and one central analyst gated them, so for a given lot
the true composition is the same everywhere and the spread of a cluster's frequency across
centers is instrument and handling variance. Does any representation lower that spread, and
does any representation raise it?

No arm lowers the between-site spread below the central manual gates, and plain intensity
matches them. At 12 metaclusters, the granularity of the 14 gated populations, the median
between-site CV of cluster frequency among T cells is 0.15 for A1 against 0.14 for the manual
gates, 0.17 for the z-score and percentile arms, 0.17 to 0.19 for the binary arms and 0.23 to
0.27 for arm B. The replicate CV within a center is 0.03 to 0.06 for every arm and 0.04 for the
gates, so the between-site component is three to five times the handling noise under every
representation, and none of them removes it. At 30 metaclusters every arm spreads more (0.26 to
0.50), because finer clusters split along the intensity axes that move between instruments, and
arm B with 2^rank divisors is the worst at 0.50. Dividing by a vendor brightness rank does not
take the instrument out of the data. The one candidate that can is a divisor measured on each
instrument from the single-stain controls on disk. The last section runs that pass, and it does
not lower the spread either: 0.30 of T cells at 12 metaclusters and 0.48 at 30.

## What ran

Seven centers (Baylor, CIMR, Miami, NHLBI, Stanford, UCLA, Yale; BSMS has no central gates),
three SeraCare lots (12828, 1349, 1369), three replicates per lot and center, 63 files. The
deposited files hold 90,000 to 480,000 events each; the loader takes 40,000 events per file
with a fixed seed, so every center and lot weighs the same in the pooled clustering, 2.5
million cells in all. Seven markers: CCR7, CD4, CD45RA, CD3 and CD8 as lineage bits for the
binary arms, CD38 and HLA-DR continuous, Live kept aside. FlowSOM 10 by 10, 12 metaclusters
(the eight memory subsets plus activated, non-T, dead and debris) and 30, seed 1 on the full
pool, 10 seeds on 50,000 cells for stability.

The deposited files are not raw. flowWorkspace wrote them from the FlowJo 9 workspaces:
compensated with the matrix the analyst edited, transformed with FlowJo's biexponential and
reordered into the HIPC marker order, with the fluorochrome keyword left stale. Arm A1 is
therefore the deposited transform as it stands, and arms A2, A3, C2 and D build on it as on
every other dataset. Arm B needs the linear scale, so it inverts the transform.

## The transform, and what it means for arm B

Stanford deposited its raw Diva exports beside the processed files, and the event order
matches, so the transform can be fitted event by event. It is a logicle with T 262144, M 4.5
decades, no extra negative decades, scaled to 4096, and a width that FlowJo chose per
parameter: Live 0.40, CCR7 0.90, CD4 0.80, CD45RA 0.80, CD3 0.15, HLA-DR 0.75, CD38 0.60, CD8
1.10. With the workspace matrix and those widths the median error is 1 to 6 units on the 4096
scale and the same on two files (`results/hipc_transform_check.md`). The iplogicle in the HIPC
README (T 16777215, M 7.22, W 1.16) is the automated pipeline's transform and does not
reproduce the files.

Only Stanford deposited raw files, so the other six centers take Stanford's widths. The cost
of that assumption, measured on Stanford (`results/hipc_width_check.md`, FlowSOM 12
metaclusters, seed 1):

| Comparison | ARI |
|---|---|
| Arm B with the fitted widths against arm B with one flat width 0.5 | 0.80 |
| Arm B with the fitted widths against arm B with all widths doubled | 0.60 |
| Arm A1 against arm B with the fitted widths | 0.82 |
| Arm A1 run-to-run stability, 10 seeds | 0.80 |

The widths move arm B about as much as arm B moves the clustering away from A1, and both
sit at the seed-to-seed noise of FlowSOM on this panel. Two consequences. On seven markers
with this brightness spread, arm B at 12 metaclusters is within seed noise of A1. And for
the six centers without raw files, the arm B result cannot be separated from the inverse-
transform assumption, so this report gives arm B's between-site numbers with that caveat
and does not build a claim on them.

## The reference: central manual gates across centers

The manual population statistics (`data/HIPC/popstats_manual/popstats-tcell.csv`) give each
population's count per file. Two denominators matter. Relative to all events, the median
between-site CV over populations at 1 percent or more is 0.25, and most of it is the
lymphocyte gate: the debris fraction differs by center. Relative to the CD3 parent the CVs are
0.04 (CD4), 0.06 (CD8), 0.09 to 0.18 for the memory subsets and 0.3 to 0.7 for the
populations under 1 percent, which is what Finak 2016 reported. The replicate CV within a
center is 0.05. So the clustering comparison uses both denominators: all events, which every
arm inherits, and the arm's own T cell clusters, which is the Finak number.

## Results

### Stability and agreement between arms

| Arm | Stability ARI, 12 metaclusters | Stability ARI, 30 | ARI against A1, 12 | ARI against A1, 30 | Occupied phenotypes | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|
| A1 logicle as deposited | 0.80 | 0.74 |  |  |  |  |
| A2 z-score | 0.74 | 0.77 | 0.89 | 0.73 |  |  |
| A3 percentile | 0.81 | 0.73 | 0.63 | 0.72 |  |  |
| B brightness, 2^rank | 0.71 | 0.64 | 0.78 | 0.63 |  |  |
| B brightness, linear rank | 0.75 | 0.71 | 0.76 | 0.73 |  |  |
| C2 lineage mixture + z-score | 0.64 | 0.69 | 0.68 | 0.56 | 32 of 32 | 18 |
| D all mixture | 0.55 | 0.76 | 0.45 | 0.43 | 128 of 128 | 45 |

Seven markers give the binary arms 32 and 128 possible phenotypes, and the cells occupy all of
them; 18 and 45 phenotypes hold 95 percent of the cells. This is the clean-panel picture from
the plan, and it is the one conventional flow panel in the study where a bit per marker is not
noise. The stability of the binary arms is still the lowest at 12 metaclusters.

### Between-site variance at 12 metaclusters, the gated granularity

CD3 cut on the deposited scale for the T cell denominator: 1859

| Arm | Clusters | At 1 percent or more | Median between-site CV, all events | Lot 12828 | Lot 1349 | Lot 1369 | Median within-site CV, replicates | T cell fraction | Clusters at 1 percent of T cells | Median between-site CV, of T cells | Median within-site CV, of T cells | Median r to best cluster |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Central manual gates | 14 | 11 | 0.25 | 0.33 | 0.24 | 0.18 | 0.05 |  | 12 | 0.14 | 0.04 |  |
| A1 logicle as deposited | 12 | 8 | 0.20 | 0.25 | 0.19 | 0.16 | 0.06 | 0.48 | 6 | 0.15 | 0.04 | 0.67 |
| A2 z-score | 12 | 8 | 0.21 | 0.27 | 0.19 | 0.15 | 0.04 | 0.47 | 6 | 0.17 | 0.04 | 0.65 |
| A3 percentile | 12 | 11 | 0.21 | 0.26 | 0.19 | 0.15 | 0.04 | 0.47 | 6 | 0.17 | 0.03 | 0.69 |
| B brightness, 2^rank | 12 | 8 | 0.23 | 0.24 | 0.18 | 0.21 | 0.06 | 0.48 | 5 | 0.23 | 0.03 | 0.62 |
| B brightness, linear rank | 12 | 7 | 0.22 | 0.26 | 0.20 | 0.19 | 0.04 | 0.50 | 6 | 0.27 | 0.03 | 0.63 |
| C2 lineage mixture + z-score | 12 | 9 | 0.26 | 0.44 | 0.26 | 0.21 | 0.05 | 0.40 | 5 | 0.19 | 0.04 | 0.68 |
| D all mixture | 12 | 8 | 0.32 | 0.32 | 0.35 | 0.28 | 0.04 | 0.37 | 6 | 0.17 | 0.05 | 0.69 |

The full 12-metacluster tables, including the manual reference and the matched clusters per
population, are in `results/hipc_between_site_k12.md`.

### Between-site variance at 30 metaclusters

CD3 cut on the deposited scale for the T cell denominator: 1859

| Arm | Clusters | At 1 percent or more | Median between-site CV, all events | Lot 12828 | Lot 1349 | Lot 1369 | Median within-site CV, replicates | T cell fraction | Clusters at 1 percent of T cells | Median between-site CV, of T cells | Median within-site CV, of T cells | Median r to best cluster |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Central manual gates | 14 | 11 | 0.25 | 0.33 | 0.24 | 0.18 | 0.05 |  | 12 | 0.14 | 0.04 |  |
| A1 logicle as deposited | 30 | 20 | 0.47 | 0.50 | 0.49 | 0.25 | 0.08 | 0.46 | 13 | 0.30 | 0.06 | 0.77 |
| A2 z-score | 30 | 17 | 0.30 | 0.37 | 0.28 | 0.21 | 0.08 | 0.47 | 13 | 0.30 | 0.07 | 0.78 |
| A3 percentile | 30 | 23 | 0.41 | 0.48 | 0.41 | 0.27 | 0.08 | 0.47 | 16 | 0.35 | 0.06 | 0.81 |
| B brightness, 2^rank | 30 | 19 | 0.55 | 0.55 | 0.54 | 0.55 | 0.07 | 0.47 | 15 | 0.50 | 0.06 | 0.66 |
| B brightness, linear rank | 30 | 21 | 0.32 | 0.41 | 0.29 | 0.27 | 0.08 | 0.46 | 15 | 0.26 | 0.06 | 0.78 |
| C2 lineage mixture + z-score | 30 | 24 | 0.45 | 0.48 | 0.38 | 0.47 | 0.09 | 0.45 | 14 | 0.36 | 0.06 | 0.74 |
| D all mixture | 30 | 16 | 0.38 | 0.41 | 0.37 | 0.28 | 0.06 | 0.46 | 13 | 0.28 | 0.05 | 0.76 |

### Manual populations: between-site CV per lot (the reference)

| Population | Mean frequency of all events | Lot 12828, all events | Lot 1349, all events | Lot 1369, all events | Within-site CV | Lot 12828, of CD3 | Lot 1349, of CD3 | Lot 1369, of CD3 |
|---|---|---|---|---|---|---|---|---|
| Lymphocytes | 0.510 | 0.23 | 0.18 | 0.12 | 0.04 |  |  |  |
| CD3 | 0.340 | 0.30 | 0.23 | 0.15 | 0.05 |  |  |  |
| CD4 | 0.185 | 0.29 | 0.24 | 0.12 | 0.06 | 0.05 | 0.03 | 0.04 |
| CD8 | 0.123 | 0.32 | 0.22 | 0.18 | 0.06 | 0.04 | 0.08 | 0.06 |
| CD4 Central Memory | 0.078 | 0.33 | 0.32 | 0.19 | 0.07 | 0.09 | 0.12 | 0.10 |
| CD4 Naive | 0.070 | 0.30 | 0.22 | 0.11 | 0.07 | 0.09 | 0.13 | 0.11 |
| CD8 Naive | 0.053 | 0.56 | 0.24 | 0.19 | 0.07 | 0.27 | 0.04 | 0.13 |
| CD4 Effector Memory | 0.034 | 0.37 | 0.26 | 0.22 | 0.09 | 0.25 | 0.11 | 0.14 |
| CD8 Effector | 0.025 | 0.47 | 0.47 | 0.27 | 0.13 | 0.41 | 0.40 | 0.31 |
| CD8 Effector Memory | 0.025 | 0.36 | 0.22 | 0.27 | 0.08 | 0.24 | 0.12 | 0.32 |
| CD8 Central Memory | 0.020 | 0.65 | 0.42 | 0.39 | 0.09 | 0.40 | 0.24 | 0.27 |
| CD4 Activated | 0.004 | 0.36 | 0.26 | 0.55 | 0.16 | 0.21 | 0.21 | 0.56 |
| CD8 Activated | 0.003 | 0.34 | 0.23 | 0.40 | 0.12 | 0.35 | 0.21 | 0.36 |
| CD4 Effector | 0.003 | 0.97 | 0.18 | 0.53 | 0.23 | 1.27 | 0.14 | 0.60 |

### Best-matching cluster per manual population, per arm (r across files; between-site CV of the cluster against the population)

| Population | Manual CV | A1 r / CV | A2 r / CV | A3 r / CV | B_geom r / CV | B_lin r / CV | C2 r / CV | D r / CV |
|---|---|---|---|---|---|---|---|---|
| Lymphocytes | 0.18 | 0.74 / 0.18 | 0.72 / 0.16 | 0.73 / 0.24 | 0.65 / 0.55 | 0.74 / 0.16 | 0.71 / 0.56 | 0.77 / 0.25 |
| CD3 | 0.23 | 0.79 / 0.18 | 0.76 / 0.16 | 0.80 / 0.24 | 0.65 / 0.73 | 0.77 / 0.16 | 0.75 / 0.25 | 0.83 / 0.25 |
| CD4 | 0.24 | 0.90 / 0.18 | 0.87 / 0.16 | 0.87 / 0.24 | 0.77 / 0.16 | 0.89 / 0.16 | 0.83 / 0.18 | 0.89 / 0.25 |
| CD8 | 0.22 | 0.85 / 0.19 | 0.85 / 0.19 | 0.85 / 0.20 | 0.84 / 0.32 | 0.85 / 0.22 | 0.90 / 0.22 | 0.87 / 0.24 |
| CD4 Central Memory | 0.32 | 0.87 / 0.18 | 0.85 / 0.16 | 0.88 / 0.24 | 0.66 / 0.19 | 0.89 / 0.16 | 0.82 / 0.25 | 0.88 / 0.25 |
| CD4 Naive | 0.22 | 0.90 / 0.22 | 0.91 / 0.20 | 0.91 / 0.28 | 0.88 / 0.19 | 0.92 / 0.20 | 0.95 / 0.18 | 0.95 / 0.22 |
| CD8 Naive | 0.24 | 0.96 / 0.19 | 0.95 / 0.19 | 0.96 / 0.20 | 0.96 / 0.32 | 0.96 / 0.22 | 0.97 / 0.22 | 0.75 / 0.49 |
| CD4 Effector Memory | 0.26 | 0.75 / 0.23 | 0.81 / 0.59 | 0.83 / 0.56 | 0.41 / 0.73 | 0.79 / 0.18 | 0.68 / 0.34 | 0.52 / 0.16 |
| CD8 Effector | 0.47 | 0.89 / 0.50 | 0.89 / 0.50 | 0.90 / 0.48 | 0.86 / 0.47 | 0.81 / 0.29 | 0.86 / 0.52 | 0.72 / 0.24 |
| CD8 Effector Memory | 0.27 | 0.56 / 0.25 | 0.58 / 0.81 | 0.56 / 0.86 | 0.67 / 0.29 | 0.53 / 0.22 | 0.71 / 0.49 | 0.81 / 0.42 |
| CD8 Central Memory | 0.42 | 0.68 / 2.12 | 0.68 / 2.18 | 0.68 / 2.07 | 0.65 / 2.60 | 0.68 / 2.21 | 0.74 / 0.45 | 0.64 / 1.32 |
| CD4 Activated | 0.36 | 0.41 / 0.47 | 0.36 / 0.25 | 0.41 / 0.18 | 0.52 / 0.47 | 0.52 / 0.47 | 0.56 / 0.48 | 0.55 / 0.51 |
| CD8 Activated | 0.34 | 0.58 / 0.19 | 0.61 / 0.19 | 0.59 / 0.34 | 0.59 / 0.27 | 0.61 / 0.27 | 0.59 / 0.33 | 0.73 / 0.58 |
| CD4 Effector | 0.53 | 0.71 / 0.23 | 0.73 / 0.21 | 0.70 / 0.29 | 0.28 / 0.73 | 0.69 / 0.18 | 0.69 / 0.46 | 0.74 / 0.63 |

### Reading the two granularities together

At 12 metaclusters the clusters are about as coarse as the gated populations, and A1 reproduces
the manual spread (0.15 against 0.14 of T cells; 0.20 against 0.25 of all events, where the
clustering has no lymphocyte gate to inherit). At 30 metaclusters the clusters cut the CD45RA,
CCR7, CD38 and HLA-DR axes finer than the analyst did, and the between-site CV doubles for every
arm. The between-site component sits in the intensity positions along those axes, which shift
by instrument, and a finer partition reads more of it. That is the mechanism the paper has
to state: the spread is a property of where the populations sit on each instrument, and a
representation change moves it only if it moves those positions. Arm B with vendor ranks
moves every channel by the same factor at every site, so it cannot, and its wider spread at
both granularities comes from the dim channels it multiplies up. The matched-cluster
correlations across the 63 files (0.62 to 0.81, median per arm) are driven by the three lots
and separate the arms little; the per-population columns in the k12 and k30 files show the
same order as the CVs.


## Hungarian F1 against the central manual gates (added 2026-10-19)

**Run:** `code/pilot.py . --datasets HIPC_tcell --arms A1,A2,A3,C2,D --tag hipc_f1` · **Labels:** `code/flowjo9_gates.py . --write`, `data/HIPC/labels/<Center>/<file>.npz` · **Metrics:** `results/hipc_f1_metrics.csv`, `results/hipc_f1_perpop.csv`, `results/hipc_gate_agreement.csv`, `results/hipc_gate_agreement.md`

### How the labels were made

The central gates are in the FlowJo 9.6.3 workspaces as nested polygon and rectangle gates with vertices on the linear data scale, and the same workspaces carry the biexponential of each channel as a calibration table of 4097 points. `code/flowjo9_gates.py` maps every vertex into display space through that table and tests the deposited values, which are already in display space. The table reproduces the deposited Stanford values from the raw files with zero error, so the gates land where the analyst drew them. The eight leaf populations (CD4 and CD8 naive, central memory, effector memory and effector) become labels 1 to 8, every other event is 0, and the CD38+ HLA-DR+ activated events carry a separate flag.

The counts agree with the central popstats to the event at CIMR and Stanford, and within 0.35 percent at Yale, where 0.6 percent of events sit above the display top on some marker. Baylor, Miami, NHLBI and UCLA run 0.1 to 1.2 percent high on CD3, because their singlets gate uses FSC-H and the deposited files carry only FSC-A and SSC-A; the labels at those four centers keep the doublets the analyst removed. The small activated populations run up to 19 percent high at Miami for the same reason. A doublet of a T cell and an HLA-DR+ cell looks activated. The full table per center and population is in `results/hipc_gate_agreement.md`.

Of the 2.52 million pooled cells, 0.78 million (31 percent) carry a label; the scorer uses those and ignores the rest, as on every other dataset. The pooled clustering is the same one as in the between-site section. It is FlowSOM 10 by 10 with 12 metaclusters, seed 1, on all seven centers at once.

### Result

| Arm | Mean F1 | ARI | NMI | Stability ARI |
|---|---|---|---|---|
| A1 | 0.401 | 0.604 | 0.687 | 0.803 |
| A2 | 0.411 | 0.636 | 0.727 | 0.736 |
| A3 | 0.421 | 0.658 | 0.738 | 0.813 |
| C2 | 0.395 | 0.507 | 0.576 | 0.644 |
| D | 0.367 | 0.475 | 0.546 | 0.554 |

Per population, F1 of the matched cluster (n is the number of labeled cells in the pool):

| Population | n | A1 | A2 | A3 | C2 | D |
|---|---|---|---|---|---|---|
| CD4 Naive | 174,583 | 0.912 | 0.944 | 0.953 | 0.859 | 0.867 |
| CD4 Central Memory | 198,298 | 0.732 | 0.775 | 0.774 | 0.670 | 0.717 |
| CD4 Effector Memory | 86,504 | 0.054 | 0.051 | 0.009 | 0.041 | 0.205 |
| CD4 Effector | 7,880 | 0.030 | 0.027 | 0.014 | 0.005 | 0.034 |
| CD8 Naive | 136,292 | 0.815 | 0.961 | 0.968 | 0.751 | 0.641 |
| CD8 Central Memory | 50,079 | 0.011 | 0.011 | 0.002 | 0.018 | 0.033 |
| CD8 Effector Memory | 63,410 | 0.651 | 0.520 | 0.521 | 0.109 | 0.372 |
| CD8 Effector | 63,893 | 0.004 | 0.002 | 0.125 | 0.708 | 0.065 |

Every arm scores low. Mean F1 is 0.37 to 0.42, with A3 first and D last, against 0.6 to 0.8 on the single-instrument datasets of phase 1. The order of the arms is the order of the other datasets. The scaled and percentile arms sit at or above A1, and the binary arms sit below it. The spread between arms is small next to the distance to the manual gates.

The per-population column shows where the F1 goes. The two naive subsets and CD4 central memory score 0.7 to 0.97 on every intensity arm, and the five other subsets score 0.0 to 0.2 on most arms, with one exception per arm (CD8 effector memory on A1, CD8 effector on C2 and CD4 effector memory on D). The analyst cut these five subsets with a CCR7 or CD45RA threshold. That threshold sits at a different display position on each instrument, and the between-site section above measured that spread. Two explanations fit this run alone. In a pooled clustering the metaclusters split along site before they split along the threshold, or 12 metaclusters for eight subsets plus non-T, dead, debris and activated are too few to reach the quadrants at all. The subsection below tests both with a per-center clustering and a clustering at 30 metaclusters, and the second explanation holds.

This run does not score arm B. Its inverse transform uses the fitted logicle widths. The workspace calibration tables are now the exact transform, and arm B can switch to them before a later run scores it.

### Per center, and at 30 metaclusters (added 2026-10-19)

**Run:** `code/pilot.py . --datasets HIPC_tcell_<Center> --arms A1,A2,A3,C2,D --tag hipc_f1_site` (seven jobs, one center each, the same 40,000 events per file and k 12) and `--datasets HIPC_tcell_k30 --tag hipc_f1_k30` (the pooled job at 30 metaclusters) · **Metrics:** `results/hipc_f1_site_metrics.csv`, `results/hipc_f1_site_perpop.csv`, `results/hipc_f1_k30_metrics.csv`, `results/hipc_f1_k30_perpop.csv`

Mean F1 per arm. The per-center column is the mean over the seven centers, with the range in brackets. The pooled columns are the seven centers clustered at once, at 12 metaclusters (the section above) and at 30.

| Arm | Per center, mean (min to max) | Pooled, 12 | Pooled, 30 |
|---|---|---|---|
| A1 | 0.371 (0.25 to 0.44) | 0.401 | 0.709 |
| A2 | 0.395 (0.19 to 0.54) | 0.411 | 0.669 |
| A3 | 0.469 (0.40 to 0.54) | 0.421 | 0.740 |
| C2 | 0.377 (0.35 to 0.40) | 0.395 | 0.592 |
| D | 0.312 (0.18 to 0.43) | 0.367 | 0.550 |

Mean F1 and ARI per center and arm:

| Center | A1 F1 / ARI | A2 F1 / ARI | A3 F1 / ARI | C2 F1 / ARI | D F1 / ARI |
|---|---|---|---|---|---|
| Baylor | 0.30 / 0.38 | 0.19 / 0.32 | 0.40 / 0.61 | 0.39 / 0.52 | 0.43 / 0.56 |
| CIMR | 0.35 / 0.49 | 0.43 / 0.62 | 0.43 / 0.64 | 0.39 / 0.56 | 0.18 / 0.32 |
| Miami | 0.25 / 0.33 | 0.34 / 0.43 | 0.41 / 0.60 | 0.37 / 0.42 | 0.23 / 0.31 |
| NHLBI | 0.44 / 0.66 | 0.43 / 0.67 | 0.54 / 0.72 | 0.35 / 0.41 | 0.36 / 0.32 |
| Stanford | 0.42 / 0.65 | 0.43 / 0.69 | 0.43 / 0.68 | 0.38 / 0.54 | 0.37 / 0.41 |
| UCLA | 0.44 / 0.63 | 0.54 / 0.67 | 0.53 / 0.70 | 0.40 / 0.52 | 0.33 / 0.39 |
| Yale | 0.39 / 0.62 | 0.40 / 0.60 | 0.54 / 0.78 | 0.35 / 0.48 | 0.28 / 0.40 |

Per population, F1 of the matched cluster: the mean over the seven per-center runs, against the pooled run at 12 and at 30 metaclusters.

| Population | A1 site / 12 / 30 | A2 site / 12 / 30 | A3 site / 12 / 30 | C2 site / 12 / 30 | D site / 12 / 30 |
|---|---|---|---|---|---|
| CD4 Naive | 0.74 / 0.91 / 0.91 | 0.69 / 0.94 / 0.94 | 0.95 / 0.95 / 0.93 | 0.83 / 0.86 / 0.86 | 0.60 / 0.87 / 0.87 |
| CD4 Central Memory | 0.72 / 0.73 / 0.86 | 0.62 / 0.77 / 0.90 | 0.79 / 0.77 / 0.86 | 0.47 / 0.67 / 0.72 | 0.65 / 0.72 / 0.72 |
| CD4 Effector Memory | 0.03 / 0.05 / 0.86 | 0.19 / 0.05 / 0.84 | 0.15 / 0.01 / 0.86 | 0.24 / 0.04 / 0.58 | 0.10 / 0.21 / 0.21 |
| CD4 Effector | 0.09 / 0.03 / 0.57 | 0.03 / 0.03 / 0.11 | 0.03 / 0.01 / 0.67 | 0.02 / 0.01 / 0.09 | 0.07 / 0.03 / 0.04 |
| CD8 Naive | 0.76 / 0.82 / 0.96 | 0.89 / 0.96 / 0.96 | 0.93 / 0.97 / 0.97 | 0.83 / 0.75 / 0.86 | 0.57 / 0.64 / 0.91 |
| CD8 Central Memory | 0.07 / 0.01 / 0.02 | 0.10 / 0.01 / 0.56 | 0.18 / 0.00 / 0.57 | 0.14 / 0.02 / 0.43 | 0.05 / 0.03 / 0.52 |
| CD8 Effector Memory | 0.40 / 0.65 / 0.66 | 0.23 / 0.52 / 0.18 | 0.28 / 0.52 / 0.19 | 0.24 / 0.11 / 0.49 | 0.22 / 0.37 / 0.42 |
| CD8 Effector | 0.15 / 0.00 / 0.84 | 0.42 / 0.00 / 0.85 | 0.45 / 0.13 / 0.87 | 0.25 / 0.71 / 0.70 | 0.25 / 0.06 / 0.70 |

The per-center run does not support the site-split explanation given above. Clustered one center at a time at 12
metaclusters, the arms score 0.31 to 0.47, the same range as the pooled run. The five low subsets average 0.17 per
center against 0.15 pooled. The same cells clustered at 30 metaclusters score 0.55 to 0.74, and the five low subsets
rise to 0.51.

The limit at 12 metaclusters is the number of clusters, not the mixing of instruments. FlowSOM spends
the 12 metaclusters on the large divisions first (CD4 against CD8, T against non-T, dead, debris), and the CCR7 and
CD45RA quadrants inside each lineage get one cluster each at best. At 30 the quadrants get their own clusters, and
the site spread measured in the between-site section does not pull them apart.

At 30 metaclusters the arms separate more than at 12. A3 scores 0.74, A1 0.71, A2 0.67, C2 0.59 and D 0.55. The
intensity arms recover CD4 effector memory and CD8 effector at 0.84 to 0.87, the binary arm D stays at 0.21 and 0.70
on the same two subsets, and arm C2 sits between them. This is the order the single-instrument datasets gave. The
per-center result is the stronger of the two for the paper, because it holds across seven instruments with one
analyst's gates. Where a per-center number is low for one arm and center (A2 at Baylor 0.19, D at CIMR 0.18), the
cluster of a large subset absorbed its neighbor; the per-population file carries each case.

## Caveats

- One seed on the full pool, as in every run; the stability column is the seed effect.
- Arm B's divisors in the phase 3 tables are vendor ranks (AmCyan rated 2 as V500, the rest
  from the fluor-ai table), and the inverse transform for six centers rests on Stanford's
  widths. The section below runs arm B with divisors measured on each center's single-stain
  controls (`data/HIPC/COMP`); the width assumption applies to those arms too.
- Per-cell manual labels exist since 2026-10-19 (`code/flowjo9_gates.py`, section below); at
  Baylor, Miami, NHLBI and UCLA they keep the doublets the analyst removed, because the deposited
  files carry no FSC-H.
- The 1 percent floor for the CV medians is a choice; `results/hipc_between_site_clusters.csv`
  holds every cluster.
- About 0.6 percent of events carry absurd negative values from the depositor's export and
  are dropped before anything else.

## Arm B with a divisor measured on each instrument

**Run:** 2026-10-19 · **Code:** `code/hipc_stain_index.py`, then `code/pilot.py --datasets HIPC_tcell,HIPC_tcell_k30 --arms A1,B_geom,B_site,B_pool --tag hipc_site`, then `code/between_site.py --prefix hipc_site_k12` and `--prefix hipc_site_k30` · **Metrics:** `results/hipc_stain_index.csv`, `results/hipc_stain_index.md`, `results/hipc_site_metrics.csv`, `results/hipc_site_meta.json` (the divisors), `results/hipc_site_k12.csv`, `results/hipc_site_k30.csv`, `results/hipc_site_k12.md`, `results/hipc_site_k30.md`

The phase 3 caveat names a divisor measured on each instrument as the one version of arm B
that can move the between-site spread. This section measures it. Each center deposited its
single-stain compensation controls (`data/HIPC/COMP/<Center>/`), raw FACSDiva exports on the
linear scale. The antibody controls are BD CompBeads with about half positive beads in the
tube, and the Live/Dead FITC controls are cells. For each center and dye the script takes the
control with the panel's own reagent when the center has one, splits the tube on its primary
detector with a two-component Gaussian mixture on arcsinh(x / 150), and computes the stain
index on the linear scale. The index is the median positive minus the median negative, over
twice the robust SD of the negative. On the cell tubes the positive is the brightest of three components, because
the two-component cut sits inside the negative when the positive fraction is small. Six
center-and-dye pairs have no control: Baylor FITC and V500, NHLBI PE, UCLA PerCP-Cy5.5, V450
and APC. Those cells take the mean over the other centers.

### Stain index per center

The control marker is in parentheses. "none, pooled" marks a center without a control for that
dye. The full table with detector, voltage, event counts and the negative for each tube is in
`results/hipc_stain_index.md`.

| Dye | Baylor | CIMR | Miami | NHLBI | Stanford | UCLA | Yale | Mean over centers | Brightest over dimmest |
|---|---|---|---|---|---|---|---|---|---|
| FITC | none, pooled 144 | 91 (Live, cells) | 98 (Live, cells) | 196 (nan, cells) | 128 (Live, cells) | 138 (Live, cells) | 213 (Live, cells) | 144 | 2.3 |
| PE | 618 (CCR7) | 487 (CCR7) | 600 (CCR7) | none, pooled 579 | 281 (CCR7) | 960 (CCR7) | 530 (nan, cells) | 579 | 3.4 |
| PerCP-Cy5.5 | 30 (CD4) | 25 (CD4) | 19 (CD4) | 26 (nan) | 19 (CD4) | none, pooled 23 | 21 (CD4) | 23 | 1.6 |
| PE-Cy7 | 469 (CD45RA) | 504 (CD45RA) | 158 (CD45RA) | 590 (CD45RA) | 255 (CD194) | 619 (CD194) | 608 (CD45RA) | 457 | 3.9 |
| V450 / Pacific Blue | 112 (CD3) | 59 (CD3) | 98 (CD3) | 110 (nan) | 142 (CD3) | none, pooled 107 | 121 (CD3) | 107 | 2.4 |
| V500 / AmCyan | none, pooled 90 | 44 (HLA-DR) | 130 (HLA-DR) | 97 (nan) | 115 (HLA-DR) | 135 (HLA-DR) | 18 (HLA-DR) | 90 | 7.7 |
| APC | 253 (CD38) | 150 (CD38) | 139 (CD38) | 298 (CD38) | 283 (CD38) | none, pooled 228 | 244 (CD38) | 228 | 2.1 |
| APC-H7 / APC-Cy7 | 98 (CD20) | 67 (CD8) | 56 (CD8) | 83 (CD20) | 51 (CD20) | 138 (CD8) | 113 (CD8) | 87 | 2.7 |

The index carries the conjugate and the detector gain, and not the antigen, because the beads
capture the antibody. The control antibodies differ by center: Baylor, NHLBI and Stanford have
no CD8 control and take CD20 on APC-H7, Stanford and UCLA have no CD45RA control and take
CD194 on PE-Cy7, and the Yale PE control is a cell tube without a marker name. Part of a
between-center difference on one dye is therefore the antibody and not the instrument. The
pooled index agrees with the vendor rank on this panel (Spearman 0.77 on the seven feature
dyes), unlike the DURAClone run in `results/stain_index.md`, but the measured spread from the
brightest to the dimmest dye is 25 against 8 for the 2^rank map. The one value that stands out
is HLA-DR V500 at Yale, 18 against 44 to 135 elsewhere. The negative beads at Yale spread
eight times wider on that detector, and the unstained tube shows the same spread, so the
instrument and not the tube sets that value.

### Two measured arms

B_pool divides every center by the across-center mean index per dye, normalized to geometric
mean 1 over the seven markers, so B_pool against B_geom is the measured-versus-rank effect.
B_site divides each cell by the index of its own center, normalized the same way within the
center, so B_site against B_pool is the per-instrument effect. Both go through the per-channel
inverse logicle of phase 3 (`divide_channels` in `code/pilot.py`), with the same Stanford
widths at every center. The divisors are in `results/hipc_site_meta.json`. Between centers
they differ by a factor of 1.8 (CD4) to 3.4 (CD45RA), and 8.3 for HLA-DR because of Yale.
Against the 2^rank divisors the pooled measured ones divide CD4 by 0.16 instead of 0.61 and
HLA-DR by 0.63 instead of 0.31; the other five are within a factor of 1.7.

| Arm | 12 metaclusters: between-site CV, all events | of T cells | within-site CV, of T cells | 30 metaclusters: between-site CV, all events | of T cells | within-site CV, of T cells | Stability ARI, 12 | Stability ARI, 30 | ARI against A1, 12 | ARI against A1, 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| Central manual gates | 0.25 | 0.14 | 0.04 | 0.25 | 0.14 | 0.04 |  |  |  |  |
| A1 logicle as deposited | 0.20 | 0.15 | 0.04 | 0.47 | 0.30 | 0.06 | 0.80 | 0.73 |  |  |
| B brightness, 2^rank (phase 3) | 0.23 | 0.23 | 0.03 | 0.55 | 0.50 | 0.06 | 0.71 | 0.67 | 0.78 | 0.63 |
| B_pool brightness, stain index pooled over centers | 0.25 | 0.28 | 0.03 | 0.44 | 0.30 | 0.05 | 0.76 | 0.71 | 0.64 | 0.58 |
| B_site brightness, stain index per center | 0.24 | 0.30 | 0.05 | 0.52 | 0.48 | 0.06 | 0.82 | 0.72 | 0.67 | 0.51 |

### Does a per-instrument divisor move the spread?

No. At 12 metaclusters the median between-site CV of T cells is 0.15 for A1, 0.23 for the
2^rank arm, 0.28 for B_pool and 0.30 for B_site. At 30 metaclusters it is 0.30, 0.50, 0.30 and
0.48. The pooled measured divisor matches A1 at 30 metaclusters and is the lowest of the four
arms on all events there (0.44 against 0.47), but it is wider than A1 at 12. A measured divisor
is therefore no better than a rank divisor against the deposited scale, and the measured-
versus-rank effect goes in both directions by granularity. The per-center divisor is wider than
the pooled one at both granularities, by 0.02 at 12 metaclusters and 0.18 at 30, and its
stability is the highest of the four arms at 12 metaclusters (0.82), so the extra spread is not
seed noise. The replicate CV within a center stays at 0.03 to 0.06 for every arm. A per-center
divisor scales the channels of each center by a different factor, and the index of a bead tube
carries the conjugate lot, the antibody and the gain, so that factor moves the populations of
each center by a different amount and the between-site spread goes up. On this panel the
between-site component does not behave as a per-channel gain that a single-stain control can
measure and remove.

