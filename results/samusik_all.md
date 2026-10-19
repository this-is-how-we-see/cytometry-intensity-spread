# Samusik_all — all 10 mice, and between-mouse variance per arm

**Run:** 2026-10-18 · **Code:** `code/pilot.py --datasets Samusik_all,Samusik_all_k40 --arms A1,A2,A3,C2,D --tag samusik_all`, then `code/between_sample.py`, tables by `code/phase2_report.py --tag samusik_all` · **Metrics:** `results/samusik_all_metrics.csv`, `results/samusik_all_perpop.csv`, `results/samusik_all_between.csv`, `results/samusik_all_between_pops.csv`, `results/samusik_all_between_clusters.csv`

## The question

Two questions. First, do the Samusik_01 results hold when the same panel runs on all 10 mice,
841,644 cells, clustered together? Second, the plan reserves between-site variance of cluster
frequencies for the HIPC data. The 10 mice are 10 samples on one instrument, so they give a
first look at that metric with data already on disk. Does any arm lower the between-sample
variance, and does the metric mean what the plan hopes?

The Samusik_01 results hold on 10 mice. Binary input loses at both cluster numbers, z-scoring
costs accuracy at the free cluster number, and the intermediate monocytes are the population
that breaks. The between-mouse coefficient of variation does not separate the arms, and the
reason is in the data: the mice differ from each other, the manual labels show a median CV
of 0.19 across mice, and a clustering cannot go below that without mixing populations. The
one arm that does go below it, C2 at 0.15, is the least accurate arm. The metric that does
separate the arms is the correlation across mice between a population's frequency and its
matched cluster's frequency, and that metric ranks the arms in the same order as F1. The HIPC
design still matters, because there the same sample ran at nine sites, the true frequency is
constant, and the CV is the right measurement.

## What ran

Samusik_all: 841,644 cells from 10 mice, 514,386 labeled with 24 populations, 39 markers, the
same lineage split as Samusik_01 (28 lineage bits, 11 continuous). Cofactor 5. FlowSOM 10 by
10, 24 and 40 metaclusters, one seed on the full data, 10 seeds on 50,000 cells for stability.
All 10 mice were clustered as one file, which is how a cross-sample study pools its data.

## Accuracy against the manual labels

### Samusik_all, FlowSOM, 24 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.58 | 0.64 | 0.78 | 0.72 |  |  |  |  |
| A2 z-score | 0.59 | 0.66 | 0.77 | 0.80 | 0.77 |  |  |  |
| A3 percentile | 0.62 | 0.69 | 0.81 | 0.75 | 0.80 |  |  |  |
| C2 lineage mixture + z-score | 0.28 | 0.33 | 0.50 | 0.67 | 0.36 | 450,229 | 268,435,456 | 408,147 |
| D all mixture | 0.47 | 0.62 | 0.67 | 0.69 | 0.46 | 732,207 | 549,755,813,888 | 690,125 |

### Samusik_all, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.67 | 0.88 | 0.85 | 0.77 |  |  |  |  |
| A2 z-score | 0.62 | 0.63 | 0.75 | 0.75 | 0.57 |  |  |  |
| A3 percentile | 0.68 | 0.78 | 0.81 | 0.69 | 0.71 |  |  |  |
| C2 lineage mixture + z-score | 0.39 | 0.35 | 0.53 | 0.58 | 0.32 | 450,047 | 268,435,456 | 407,965 |
| D all mixture | 0.55 | 0.32 | 0.60 | 0.45 | 0.32 | 732,570 | 549,755,813,888 | 690,488 |

### Per-population F1, Samusik_all, 24 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D |
|---|---|---|---|---|---|---|
| IgDpos IgMpos B cells | 114,412 | 0.94 | 0.97 | 0.96 | 0.32 | 0.95 |
| Classical Monocytes | 99,356 | 0.67 | 0.69 | 0.72 | 0.58 | 0.72 |
| Intermediate Monocytes | 82,384 | 0.00 | 0.23 | 0.28 | 0.15 | 0.54 |
| Eosinophils | 53,495 | 0.93 | 0.91 | 0.93 | 0.46 | 0.61 |
| IgD- IgMpos B cells | 49,717 | 0.64 | 0.75 | 0.73 | 0.37 | 0.69 |
| pDCs | 26,421 | 0.71 | 0.86 | 0.96 | 0.77 | 0.95 |
| IgM- IgD- B-cells | 21,137 | 0.00 | 0.01 | 0.01 | 0.34 | 0.01 |
| CD8 T cells | 12,479 | 0.98 | 0.98 | 0.98 | 0.69 | 0.63 |
| NKT cells | 8,479 | 0.90 | 0.87 | 0.91 | 0.06 | 0.00 |
| Non-Classical Monocytes | 7,855 | 0.00 | 0.00 | 0.97 | 0.91 | 0.81 |
| GMP | 6,068 | 0.36 | 0.50 | 0.58 | 0.02 | 0.50 |
| CD4 T cells | 5,646 | 0.93 | 0.93 | 0.93 | 0.27 | 0.00 |
| MEP | 5,469 | 0.90 | 0.80 | 0.90 | 0.17 | 0.82 |
| B-cell Frac A-C (pro-B cells) | 4,241 | 0.06 | 0.01 | 0.01 | 0.04 | 0.23 |
| Basophils | 4,067 | 0.89 | 0.89 | 0.90 | 0.00 | 0.82 |
| NK cells | 3,324 | 0.80 | 0.81 | 0.82 | 0.38 | 0.73 |
| Macrophages | 2,980 | 0.87 | 0.59 | 0.87 | 0.02 | 0.58 |
| CMP | 2,041 | 0.46 | 0.00 | 0.31 | 0.21 | 0.47 |
| Plasma Cells | 1,606 | 0.72 | 0.98 | 0.97 | 0.52 | 0.71 |
| MPP | 1,229 | 0.74 | 0.51 | 0.45 | 0.39 | 0.00 |
| mDCs | 1,084 | 0.83 | 0.90 | 0.72 | 0.03 | 0.45 |
| CLP | 480 | 0.42 | 0.39 | 0.08 | 0.01 | 0.03 |
| gd T cells | 387 | 0.09 | 0.55 | 0.00 | 0.02 | 0.07 |
| HSC | 29 | 0.00 | 0.00 | 0.00 | 0.01 | 0.00 |

### Per-population F1, Samusik_all, 40 metaclusters

| Population | Cells | A1 | A2 | A3 | C2 | D |
|---|---|---|---|---|---|---|
| IgDpos IgMpos B cells | 114,412 | 0.98 | 0.92 | 0.96 | 0.46 | 0.40 |
| Classical Monocytes | 99,356 | 0.96 | 0.71 | 0.86 | 0.65 | 0.53 |
| Intermediate Monocytes | 82,384 | 0.97 | 0.23 | 0.86 | 0.17 | 0.41 |
| Eosinophils | 53,495 | 0.93 | 0.91 | 0.93 | 0.33 | 0.40 |
| IgD- IgMpos B cells | 49,717 | 0.38 | 0.75 | 0.70 | 0.42 | 0.62 |
| pDCs | 26,421 | 0.53 | 0.57 | 0.68 | 0.69 | 0.95 |
| IgM- IgD- B-cells | 21,137 | 0.43 | 0.02 | 0.24 | 0.30 | 0.50 |
| CD8 T cells | 12,479 | 0.98 | 0.98 | 0.98 | 0.80 | 0.92 |
| NKT cells | 8,479 | 0.90 | 0.87 | 0.87 | 0.73 | 0.83 |
| Non-Classical Monocytes | 7,855 | 0.99 | 0.94 | 0.97 | 0.91 | 0.81 |
| GMP | 6,068 | 0.61 | 0.80 | 0.85 | 0.44 | 0.53 |
| CD4 T cells | 5,646 | 0.78 | 0.74 | 0.93 | 0.39 | 0.84 |
| MEP | 5,469 | 0.82 | 0.81 | 0.82 | 0.51 | 0.82 |
| B-cell Frac A-C (pro-B cells) | 4,241 | 0.13 | 0.05 | 0.04 | 0.06 | 0.24 |
| Basophils | 4,067 | 0.89 | 0.73 | 0.90 | 0.22 | 0.81 |
| NK cells | 3,324 | 0.80 | 0.81 | 0.82 | 0.03 | 0.74 |
| Macrophages | 2,980 | 0.75 | 0.60 | 0.62 | 0.65 | 0.62 |
| CMP | 2,041 | 0.46 | 0.15 | 0.32 | 0.21 | 0.32 |
| Plasma Cells | 1,606 | 0.72 | 0.98 | 0.97 | 0.53 | 0.70 |
| MPP | 1,229 | 0.74 | 0.51 | 0.58 | 0.72 | 0.50 |
| mDCs | 1,084 | 0.74 | 0.90 | 0.72 | 0.03 | 0.43 |
| CLP | 480 | 0.42 | 0.38 | 0.14 | 0.03 | 0.24 |
| gd T cells | 387 | 0.09 | 0.55 | 0.42 | 0.01 | 0.06 |
| HSC | 29 | 0.02 | 0.00 | 0.03 | 0.01 | 0.00 |

Z-scoring at 40 metaclusters costs 0.05 of mean F1 and 0.25 of ARI (0.88 to 0.64). Percentile
scaling costs nothing on mean F1 (0.68 against 0.67) and 0.10 on ARI. The z-score cost is the
smallest of the four mass sets so far (0.16, 0.12, 0.10, 0.05), and the ARI cost is the largest.
The scaling statement for the paper is therefore about z-scoring, and it needs ARI beside F1,
because mean F1 by Hungarian matching gives each population one cluster and does not see the
cells that the other clusters took.

Intermediate monocytes, 82,384 cells, are the hypothesis 3 case again and the clearest one yet.
At 40 metaclusters: A1 0.97, A3 0.86, A2 0.23, C2 0.17, D 0.41. Z-scoring Ly6C, the marker
that splits the three monocyte levels, costs almost as much as binarizing it. At 24
metaclusters A1 merges them into the classical monocytes (0.00), as it did on Samusik_01.

The phenotype count says what the binary arms are doing: 28 lineage bits put the cells in
450,047 phenotypes, and 407,965 of them are needed for 95 percent of cells. On a mass cytometry
panel the mixture threshold makes a bit per noise channel, and a 10-mouse pool does not change
that.

## Between-mouse variance

### Between-mouse variance, results/labels_samusik_all_Samusik_all_k40.npz, 10 mice

| Arm | Clusters | Clusters at 1 percent or more | Median CV, labeled cells | Median CV, all cells | Share of clusters with CV over 0.5 | Mean F1, all mice pooled | Mean F1 per mouse (SD) | ARI per mouse | Median r across mice | Share of populations with r over 0.8 |
|---|---|---|---|---|---|---|---|---|---|---|
| Manual labels | 24 | 13 | 0.19 |  | 0.00 |  |  () |  |  |  |
| A1 arcsinh | 40 | 14 | 0.21 | 0.18 | 0.00 | 0.67 | 0.67 (0.01) | 0.88 | 0.96 | 0.79 |
| A2 z-score | 40 | 14 | 0.20 | 0.18 | 0.00 | 0.62 | 0.62 (0.01) | 0.64 | 0.92 | 0.71 |
| A3 percentile | 40 | 17 | 0.20 | 0.18 | 0.00 | 0.68 | 0.68 (0.01) | 0.78 | 0.95 | 0.71 |
| C2 lineage mixture + z-score | 40 | 17 | 0.15 | 0.14 | 0.00 | 0.39 | 0.39 (0.01) | 0.35 | 0.75 | 0.46 |
| D all mixture | 40 | 23 | 0.20 | 0.17 | 0.00 | 0.55 | 0.55 (0.01) | 0.32 | 0.93 | 0.67 |

### Per population: does the matched cluster track the population across mice?

| Population | Cells | CV of the population across mice | r A1 | r A2 | r A3 | r C2 | r D | F1 A1 | F1 A2 | F1 A3 | F1 C2 | F1 D |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IgDpos IgMpos B cells | 114,412 | 0.16 | 1.00 | 1.00 | 1.00 | 0.67 | 0.82 | 0.98 | 0.92 | 0.96 | 0.46 | 0.40 |
| Classical Monocytes | 99,356 | 0.22 | 1.00 | 0.82 | 0.97 | 0.79 | 0.81 | 0.96 | 0.71 | 0.86 | 0.65 | 0.53 |
| Intermediate Monocytes | 82,384 | 0.22 | 1.00 | 0.90 | 0.99 | 0.55 | 0.93 | 0.97 | 0.23 | 0.86 | 0.17 | 0.41 |
| Eosinophils | 53,495 | 0.19 | 0.99 | 0.99 | 0.99 | 0.59 | 0.51 | 0.93 | 0.91 | 0.93 | 0.33 | 0.40 |
| IgD- IgMpos B cells | 49,717 | 0.18 | 0.94 | 0.96 | 0.97 | 0.86 | 0.95 | 0.38 | 0.75 | 0.70 | 0.42 | 0.62 |
| pDCs | 26,421 | 0.07 | 0.68 | 0.52 | 0.62 | 0.94 | 0.98 | 0.53 | 0.57 | 0.68 | 0.69 | 0.95 |
| IgM- IgD- B-cells | 21,137 | 0.26 | 0.96 | 0.87 | 0.93 | 0.96 | 0.92 | 0.43 | 0.02 | 0.24 | 0.30 | 0.50 |
| CD8 T cells | 12,479 | 0.20 | 1.00 | 1.00 | 1.00 | 0.98 | 1.00 | 0.98 | 0.98 | 0.98 | 0.80 | 0.92 |
| NKT cells | 8,479 | 0.11 | 0.97 | 0.93 | 0.92 | 0.92 | 0.97 | 0.90 | 0.87 | 0.87 | 0.73 | 0.83 |
| Non-Classical Monocytes | 7,855 | 0.24 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.99 | 0.94 | 0.97 | 0.91 | 0.81 |
| GMP | 6,068 | 0.18 | 0.96 | 0.89 | 0.94 | 0.71 | 0.85 | 0.61 | 0.80 | 0.85 | 0.44 | 0.53 |
| CD4 T cells | 5,646 | 0.20 | 0.94 | 0.95 | 0.99 | 0.56 | 0.98 | 0.78 | 0.74 | 0.93 | 0.39 | 0.84 |
| MEP | 5,469 | 0.16 | 0.98 | 0.98 | 0.98 | 0.96 | 0.99 | 0.82 | 0.81 | 0.82 | 0.51 | 0.82 |
| B-cell Frac A-C (pro-B cells) | 4,241 | 0.18 | 0.88 | -0.17 | 0.07 | 0.94 | 0.73 | 0.13 | 0.05 | 0.04 | 0.06 | 0.24 |
| Basophils | 4,067 | 0.20 | 0.99 | 0.97 | 1.00 | 0.61 | 0.95 | 0.89 | 0.73 | 0.90 | 0.22 | 0.81 |
| NK cells | 3,324 | 0.36 | 1.00 | 1.00 | 1.00 | 0.51 | 0.99 | 0.80 | 0.81 | 0.82 | 0.03 | 0.74 |
| Macrophages | 2,980 | 0.34 | 0.92 | 0.91 | 0.87 | 0.97 | 0.95 | 0.75 | 0.60 | 0.62 | 0.65 | 0.62 |
| CMP | 2,041 | 0.18 | 0.43 | 0.48 | -0.25 | 0.53 | 0.22 | 0.46 | 0.15 | 0.32 | 0.21 | 0.32 |
| Plasma Cells | 1,606 | 0.90 | 0.97 | 1.00 | 1.00 | 0.99 | 0.99 | 0.72 | 0.98 | 0.97 | 0.53 | 0.70 |
| MPP | 1,229 | 0.22 | 0.91 | 0.40 | 0.80 | 0.92 | 0.46 | 0.74 | 0.51 | 0.58 | 0.72 | 0.50 |
| mDCs | 1,084 | 0.17 | 0.92 | 0.93 | 0.78 | 0.66 | 0.31 | 0.74 | 0.90 | 0.72 | 0.03 | 0.43 |
| CLP | 480 | 0.23 | 0.59 | 0.59 | 0.09 | -0.61 | -0.33 | 0.42 | 0.38 | 0.14 | 0.03 | 0.24 |
| gd T cells | 387 | 0.22 | -0.03 | 0.51 | 0.41 | 0.47 | 0.52 | 0.09 | 0.55 | 0.42 | 0.01 | 0.06 |
| HSC | 29 | 0.73 | 0.27 | 0.21 | 0.27 | -0.30 | 0.40 | 0.02 | 0.00 | 0.03 | 0.01 | 0.00 |

Read the table this way. "Median CV, labeled cells" is the median over the clusters that hold
at least 1 percent of labeled cells of the coefficient of variation of that cluster's frequency
across the 10 mice. The manual labels give 0.19: that is how much the mice differ from each
other in the gated populations, and it is the floor for an honest clustering. A1, A2, A3 and D
sit at 0.20 to 0.21, so each of them reports the biology. C2 sits at 0.15, below the floor,
and it has the lowest mean F1 (0.39): its clusters mix populations whose frequencies vary in
different directions, and the mixture averages out. A lower between-sample CV is not better
when the samples really differ.

"Median r across mice" is the Pearson correlation across the 10 mice between a population's
frequency and the frequency of its Hungarian-matched cluster, median over populations. A1
reaches 0.96 and tracks 19 of 24 populations at r above 0.8. C2 reaches 0.75 and tracks 11.
This is the number that says whether a cluster rises and falls with its population from mouse
to mouse, and it orders the arms as F1 does. The per-population table shows the usual
exceptions: a population that no arm recovers (pro B cells, HSC with 29 cells) has no
meaningful r, and a population whose cluster took a neighbor (IgM- IgD- B cells under A2, r
of -0.47 at 24 metaclusters) has a negative one.

Per-mouse mean F1 has a standard deviation of 0.01 for every arm. The accuracy of an arm does
not depend on which mouse is scored, because the clustering ran on the pool.

## What this means for phase 3

The between-sample CV is the right metric only when the true frequencies are the same across
samples. The HIPC Lyoplate set has that property, because one set of samples ran at nine sites.
Ten mice do not. So this run does not replace phase 3, and the result it does give is this: on
samples that differ, the arm that separates populations best also tracks the between-sample
biology best, and a binary arm that lowers the variance does it by blurring. The paper can
state that as the reason it reports r alongside CV.

## Caveats

- One seed on the full data, as in every run. Stability ARI per arm is in the metrics table.
- The clustering pooled all mice. Clustering each mouse alone and matching clusters across
  mice is a different design, closer to how a site-by-site analysis runs, and it was not done.
- The 1 percent floor for the CV median is a choice. At 0.5 percent the medians move by
  less than 0.02 and the order does not change (`results/samusik_all_between_clusters.csv`
  has every cluster).
- Population 25, unassigned, is 39 percent of cells and is left out of the labeled-cell
  frequencies. The all-cell CV column includes those cells in the denominator.
- The Samusik lineage split is the same judgment as in phase 2.
