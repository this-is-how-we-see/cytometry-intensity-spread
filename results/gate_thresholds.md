# Gate thresholds — the upper bound for the binary arms on Levine_32dim

**Run:** 2026-10-18 · **Code:** `code/gate_thresholds.py`, `code/pilot.py --tag gate`, `code/leiden_hamming.py`, tables by `code/gate_report.py` · **Thresholds:** `results/gate_thresholds.json` · **Metrics:** `results/gate_metrics.csv`, `results/gate_perpop.csv`, `results/leiden_metrics.csv`

## The question

The pilot binarized every marker with a two-component Gaussian mixture, and the binary arms
lost to plain arcsinh intensity on Levine_32dim at 40 metaclusters, mean F1 0.41 against 0.79.
The mixture put 148,505 distinct phenotypes among 265,627 cells, so the thresholds were noisy.
This run asks the question the plan names in section 5, method 1: when the thresholds are as
good as the manual labels can make them, does binary input beat arcsinh intensity?

The answer is no. The best binary arm ties arcsinh at the gated number of metaclusters and
loses at the free number, and it loses on the small populations. The paragraph at the end says
what that means for the paper.

## How the thresholds were made

For each of the 32 markers, a truth table names the populations the label set defines as
positive and as negative. The cut on arcsinh(x / 5) is the one that maximizes balanced
accuracy, 0.5 times the sum of the true positive rate and the true negative rate, over 2,000
candidate values between the 0.1st and 99.9th percentile of the labeled cells. Balanced
accuracy counts a 300-cell population the same as a 26,000-cell population. A population with
intermediate expression (CD4 on monocytes, CD7 on T cells, CD20 on pre B cells) is in neither
list, so the cut is not pulled toward it. Ten markers have no positive or no negative
population in the label set and keep the mixture threshold.

The truth table was assigned in this session from the population names and the per-population
median of each marker. It is open to correction. A change to the table changes the cut for that
marker and nothing else, and `code/gate_thresholds.py` reruns in under a minute.

The gate cut sits above the mixture cut on 19 of the 22 defined markers. CD45RA, Flt3 and CD13 are the exceptions. The
mixture fits its lower component to the zero spike of a mass cytometry channel and cuts inside
the noise above it. Balanced accuracy says which markers carry label information at all. Seven
defined markers stay below 0.9 under both methods, and four of them (CD133, Flt3, CD14, CD13)
stay below 0.76. CD14 has a median arcsinh of 0.1 on monocytes in this file. No threshold makes
a dim stain informative.

### Thresholds per marker

| Marker | Method | Gate cut | Mixture cut | Balanced accuracy, gate | Balanced accuracy, mixture | Positive fraction, gate | Positive fraction, mixture | Positive populations | Negative populations | Note |
|---|---|---|---|---|---|---|---|---|---|---|
| CD45RA | gate | 0.67 | 0.83 | 0.80 | 0.79 | 0.43 | 0.34 | Mature B cells, CD16- NK cells, CD16+ NK cells | CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs, Basophils, Pro B cells | T cells left out: naive and memory mix |
| CD133 | gate | 0.05 | 0.02 | 0.56 | 0.56 | 0.47 | 0.51 | CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, Monocytes | stain is dim on HSCs in this file, median arcsinh 0.2 |
| CD19 | gate | 0.91 | 0.14 | 0.99 | 0.87 | 0.20 | 0.42 | Mature B cells, Pre B cells, Pro B cells, Plasma B cells | CD4 T cells, CD8 T cells, CD16- NK cells, CD16+ NK cells, Monocytes, CD34+CD38lo HSCs, Basophils, pDCs |  |
| CD22 | gate | 0.64 | 0.14 | 0.95 | 0.86 | 0.18 | 0.39 | Mature B cells | CD4 T cells, CD8 T cells, CD16- NK cells, CD16+ NK cells, Monocytes, CD34+CD38lo HSCs, pDCs | pre B and pro B left out: dim |
| CD11b | gate | 1.01 | 0.56 | 0.93 | 0.90 | 0.24 | 0.35 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, Pre B cells, Pro B cells, CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs, CD34+CD38+CD123+ HSPCs, pDCs | NK cells and basophils left out: dim |
| CD4 | gate | 0.35 | 0.02 | 0.91 | 0.77 | 0.35 | 0.64 | CD4 T cells | CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, CD34+CD38lo HSCs, Pre B cells | monocytes and pDCs left out: intermediate, the hypothesis 3 case |
| CD8 | gate | 1.01 | 0.25 | 1.00 | 0.92 | 0.19 | 0.34 | CD8 T cells | CD4 T cells, Mature B cells, Monocytes, CD34+CD38lo HSCs, Pre B cells, pDCs | NK cells left out: CD8 dim |
| CD34 | gate | 1.33 | 0.60 | 1.00 | 0.95 | 0.05 | 0.16 | CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs, CD34+CD38+CD123+ HSPCs, Pro B cells | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, Monocytes, Basophils, pDCs, Plasma B cells, Pre B cells |  |
| Flt3 | gate | 0.12 | 0.47 | 0.57 | 0.53 | 0.48 | 0.16 | CD34+CD38+CD123+ HSPCs, pDCs | CD4 T cells, CD8 T cells, CD16- NK cells, CD16+ NK cells, Mature B cells | stain is dim, median arcsinh 0.3 on the positive populations |
| CD20 | gate | 0.41 | 0.13 | 0.90 | 0.85 | 0.22 | 0.37 | Mature B cells | CD4 T cells, CD8 T cells, CD16- NK cells, CD16+ NK cells, Monocytes, CD34+CD38lo HSCs, Pro B cells, Plasma B cells | pre B left out: dim |
| CXCR4 | mixture fallback |  | 0.84 |  |  | 0.41 | 0.41 | none | none | continuous marker, no population in the label set defines it |
| CD235ab | mixture fallback |  | 0.79 |  |  | 0.32 | 0.32 | none | none | erythroid marker, no erythroid population in the label set |
| CD45 | mixture fallback |  | 4.92 |  |  | 0.80 | 0.80 | none | none | every labeled population is positive |
| CD123 | gate | 2.07 | 0.35 | 1.00 | 0.93 | 0.03 | 0.28 | pDCs, Basophils, CD34+CD38+CD123+ HSPCs | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, CD34+CD38+CD123- HSPCs, CD34+CD38lo HSCs | monocytes left out: dim |
| CD321 | mixture fallback |  | 2.04 |  |  | 0.44 | 0.44 | none | none | continuous marker, no population in the label set defines it |
| CD14 | gate | 0.05 | 0.02 | 0.69 | 0.68 | 0.34 | 0.39 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, CD34+CD38lo HSCs, pDCs, Basophils | stain is dim on monocytes in this file, median arcsinh 0.1 |
| CD33 | gate | 0.22 | 0.16 | 0.84 | 0.84 | 0.32 | 0.37 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, Plasma B cells | basophils left out: median arcsinh 0.1 in this file |
| CD47 | mixture fallback |  | 2.36 |  |  | 0.65 | 0.65 | none | none | continuous marker, no population in the label set defines it |
| CD11c | gate | 1.04 | 0.49 | 0.99 | 0.94 | 0.23 | 0.35 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, Pre B cells, Pro B cells, CD34+CD38lo HSCs, pDCs |  |
| CD7 | gate | 1.35 | 0.09 | 1.00 | 0.91 | 0.25 | 0.53 | CD16- NK cells, CD16+ NK cells | Mature B cells, Monocytes, Basophils, pDCs, CD34+CD38lo HSCs, Plasma B cells | T cells left out: CD7 dim in this file |
| CD15 | mixture fallback |  | 0.02 |  |  | 0.62 | 0.62 | none | none | granulocyte marker, no granulocyte population in the label set |
| CD16 | gate | 1.09 | 0.24 | 1.00 | 0.93 | 0.03 | 0.18 | CD16+ NK cells | CD16- NK cells, CD4 T cells, CD8 T cells, Mature B cells, CD34+CD38lo HSCs, pDCs, Basophils | monocytes left out: classical and non-classical mix |
| CD44 | mixture fallback |  | 3.38 |  |  | 0.68 | 0.68 | none | none | continuous marker, no population in the label set defines it |
| CD38 | gate | 1.92 | 1.48 | 1.00 | 0.85 | 0.29 | 0.38 | CD34+CD38+CD123- HSPCs, CD34+CD38+CD123+ HSPCs, Plasma B cells, Pre B cells, Pro B cells | CD34+CD38lo HSCs | the label set defines CD38 only through the CD38+ against CD38lo split of the CD34+ cells |
| CD13 | gate | 0.32 | 0.33 | 0.75 | 0.75 | 0.38 | 0.38 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells | stain is dim on monocytes, median arcsinh 0.6 |
| CD3 | gate | 2.22 | 1.58 | 1.00 | 1.00 | 0.43 | 0.44 | CD4 T cells, CD8 T cells | CD16- NK cells, CD16+ NK cells, Mature B cells, Pre B cells, Pro B cells, Plasma B cells, Monocytes, Basophils, pDCs, CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs, CD34+CD38+CD123+ HSPCs |  |
| CD61 | mixture fallback |  | 0.32 |  |  | 0.30 | 0.30 | none | none | platelet marker, no megakaryocyte population in the label set |
| CD117 | gate | 0.28 | 0.02 | 0.82 | 0.73 | 0.16 | 0.45 | CD34+CD38lo HSCs, CD34+CD38+CD123- HSPCs | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, Monocytes, Pre B cells, Plasma B cells |  |
| CD49d | mixture fallback |  | 0.87 |  |  | 0.40 | 0.40 | none | none | continuous marker, no population in the label set defines it |
| HLA-DR | gate | 0.58 | 0.47 | 0.98 | 0.96 | 0.51 | 0.53 | Mature B cells, Pre B cells, Pro B cells, Monocytes, pDCs, CD34+CD38+CD123- HSPCs, CD34+CD38+CD123+ HSPCs | CD4 T cells, CD8 T cells, CD16- NK cells, CD16+ NK cells, Basophils, Plasma B cells | HSCs left out: intermediate |
| CD64 | gate | 0.75 | 0.35 | 0.99 | 0.95 | 0.23 | 0.34 | Monocytes | CD4 T cells, CD8 T cells, Mature B cells, CD16- NK cells, CD16+ NK cells, pDCs, Basophils, CD34+CD38lo HSCs |  |
| CD41 | mixture fallback |  | 0.13 |  |  | 0.39 | 0.39 | none | none | platelet marker, no megakaryocyte population in the label set |

## FlowSOM, the same clusterer as the pilot

Three new arms. C_gate binarizes the 25 lineage markers with the gate cuts and z-scores the 7
functional markers, as C2 did. D_gate binarizes all 32 markers with the gate cuts. D_gate_clean
binarizes only the 15 markers whose gate cut reaches balanced accuracy 0.9, which drops the 10
fallback markers and the 7 dim ones. D_gate_clean uses the labels twice, once to choose the
markers and once to place the cuts, so it is the upper bound for a binary panel and nothing
more. The pilot arms ran again in the same job. A1, A2 and A3 reproduce the pilot exactly. The
mixture arms moved by up to 0.06 of mean F1, because the mixture fits on a fresh subsample each
run and the mixture threshold is itself unstable.

### 14 metaclusters, the gated number

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | Occupied phenotypes | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|
| A1 arcsinh | 0.47 | 0.66 | 0.83 | 0.94 |  |  |
| A2 z-score | 0.43 | 0.63 | 0.81 | 0.89 |  |  |
| A3 percentile | 0.41 | 0.64 | 0.82 | 0.85 |  |  |
| C2 lineage mixture + z-score | 0.45 | 0.61 | 0.69 | 0.71 | 113,405 | 100,124 |
| C_gate lineage gate + z-score | 0.46 | 0.63 | 0.73 | 0.74 | 56,339 | 43,058 |
| D all mixture | 0.43 | 0.62 | 0.70 | 0.57 | 148,595 | 135,314 |
| D_gate all gate | 0.47 | 0.73 | 0.81 | 0.64 | 120,773 | 107,492 |
| D_gate_clean 15 clean markers | 0.47 | 0.78 | 0.81 | 0.72 | 3,286 | 371 |
| C lineage mixture + percentile | 0.41 | 0.58 | 0.69 | 0.59 | 113,405 | 100,124 |

### 40 metaclusters, the free number

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | Occupied phenotypes | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|
| A1 arcsinh | 0.79 | 0.80 | 0.83 | 0.79 |  |  |
| A2 z-score | 0.63 | 0.68 | 0.76 | 0.70 |  |  |
| A3 percentile | 0.64 | 0.52 | 0.74 | 0.64 |  |  |
| C2 lineage mixture + z-score | 0.47 | 0.32 | 0.58 | 0.36 | 112,505 | 99,224 |
| C_gate lineage gate + z-score | 0.48 | 0.29 | 0.59 | 0.38 | 56,339 | 43,058 |
| D all mixture | 0.47 | 0.25 | 0.57 | 0.40 | 148,294 | 135,013 |
| D_gate all gate | 0.45 | 0.29 | 0.64 | 0.41 | 120,773 | 107,492 |
| D_gate_clean 15 clean markers | 0.68 | 0.74 | 0.79 | 0.68 | 3,286 | 371 |
| C lineage mixture + percentile | 0.35 | 0.23 | 0.55 | 0.39 | 112,505 | 99,224 |

The phenotype count is the diagnostic. With all 32 gate cuts, the cells still occupy 120,773
phenotypes, because 17 of the bits are noise under any threshold. With the 15 clean markers the
cells occupy 3,286 phenotypes, and 371 of them hold 95 percent of the cells. That is the picture
the plan predicted for a 13-marker conventional panel, and on mass cytometry it appears only
after the dim and undefined channels are removed.

### Per-population F1, 40 metaclusters

| Population | Cells | A1 | C2 | C_gate | D | D_gate | D_gate_clean |
|---|---|---|---|---|---|---|---|
| CD4 T cells | 26,366 | 0.96 | 0.51 | 0.41 | 0.35 | 0.42 | 0.93 |
| Monocytes | 21,099 | 0.66 | 0.64 | 0.47 | 0.52 | 0.46 | 0.89 |
| CD8 T cells | 20,108 | 0.96 | 0.46 | 0.53 | 0.47 | 0.55 | 0.83 |
| Mature B cells | 16,520 | 0.96 | 0.58 | 0.59 | 0.54 | 0.52 | 0.81 |
| Pre B cells | 6,135 | 0.91 | 0.49 | 0.53 | 0.54 | 0.79 | 0.60 |
| CD16- NK cells | 3,905 | 0.72 | 0.45 | 0.48 | 0.51 | 0.47 | 0.79 |
| CD34+CD38+CD123- HSPCs | 3,295 | 0.63 | 0.47 | 0.47 | 0.74 | 0.60 | 0.88 |
| CD16+ NK cells | 2,248 | 0.88 | 0.36 | 0.51 | 0.55 | 0.48 | 0.84 |
| pDCs | 1,238 | 0.81 | 0.52 | 0.62 | 0.71 | 0.82 | 0.86 |
| Basophils | 1,207 | 0.90 | 0.57 | 0.63 | 0.72 | 0.81 | 0.84 |
| CD34+CD38lo HSCs | 916 | 0.74 | 0.54 | 0.37 | 0.03 | 0.40 | 0.99 |
| Pro B cells | 513 | 0.66 | 0.15 | 0.11 | 0.21 | 0.00 | 0.01 |
| Plasma B cells | 330 | 0.79 | 0.67 | 0.72 | 0.64 | 0.01 | 0.11 |
| CD34+CD38+CD123+ HSPCs | 304 | 0.55 | 0.21 | 0.23 | 0.05 | 0.01 | 0.09 |

### Per-population F1, 14 metaclusters

| Population | Cells | A1 | C2 | C_gate | D | D_gate | D_gate_clean |
|---|---|---|---|---|---|---|---|
| CD4 T cells | 26,366 | 0.72 | 0.72 | 0.72 | 0.81 | 0.88 | 0.92 |
| Monocytes | 21,099 | 0.88 | 0.94 | 0.92 | 0.74 | 0.85 | 0.92 |
| CD8 T cells | 20,108 | 0.01 | 0.25 | 0.16 | 0.78 | 0.70 | 0.90 |
| Mature B cells | 16,520 | 0.99 | 0.95 | 0.96 | 0.95 | 0.97 | 0.87 |
| Pre B cells | 6,135 | 0.93 | 0.89 | 0.90 | 0.80 | 0.88 | 0.62 |
| CD16- NK cells | 3,905 | 0.76 | 0.02 | 0.62 | 0.66 | 0.73 | 0.65 |
| CD34+CD38+CD123- HSPCs | 3,295 | 0.82 | 0.35 | 0.65 | 0.72 | 0.83 | 0.79 |
| CD16+ NK cells | 2,248 | 0.00 | 0.58 | 0.00 | 0.00 | 0.00 | 0.28 |
| pDCs | 1,238 | 0.63 | 0.06 | 0.30 | 0.04 | 0.60 | 0.00 |
| Basophils | 1,207 | 0.00 | 0.39 | 0.40 | 0.41 | 0.09 | 0.66 |
| CD34+CD38lo HSCs | 916 | 0.01 | 0.43 | 0.02 | 0.04 | 0.00 | 0.00 |
| Pro B cells | 513 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Plasma B cells | 330 | 0.79 | 0.72 | 0.72 | 0.02 | 0.01 | 0.00 |
| CD34+CD38+CD123+ HSPCs | 304 | 0.00 | 0.00 | 0.10 | 0.00 | 0.00 | 0.00 |

The large populations are where the binary arms with noisy bits lose. CD4 T cells go from 0.96
under A1 to 0.42 under D_gate at 40 metaclusters, because FlowSOM splits 26,000 cells along the
noise bits and one-to-one matching gives the population one cluster. D_gate_clean brings them
back to 0.93. The small populations are where the clean binary arm loses. Pro B cells (0.01),
plasma B cells (0.11) and CD34+CD38+CD123+ HSPCs (0.09) each differ from a larger neighbor by
one bit among 15, and FlowSOM merges them. Under A1 they differ by intensity on several markers
at once, and A1 scores 0.66, 0.79 and 0.55. Those three populations are the hypothesis 3 cases
of this dataset: a population that intensity defines and a bit cannot.

## Leiden on a kNN graph, 50,000-cell subsample

Leiden has no cluster count, so the table reports a resolution grid and the number of clusters
each resolution found. Read the rows nearest 14 and 40 clusters. A1 uses Euclidean distance.
C_gate ran twice, Euclidean on the mixed matrix and Hamming on its 25 binary columns alone.
D_gate and D_gate_clean use Hamming. Stability is the mean pairwise ARI over 10 Leiden seeds on
the same graph.

| Arm | Distance | Features | Resolution | Clusters | Mean F1 | ARI | NMI | Stability ARI, 10 seeds |
|---|---|---|---|---|---|---|---|---|
| A1 | euclidean | 32 | 0.1 | 8 | 0.43 | 0.95 | 0.91 | 0.97 |
| A1 | euclidean | 32 | 0.2 | 10 | 0.51 | 0.97 | 0.94 | 0.99 |
| A1 | euclidean | 32 | 0.3 | 11 | 0.51 | 0.97 | 0.94 | 0.99 |
| A1 | euclidean | 32 | 0.5 | 16 | 0.56 | 0.90 | 0.90 | 0.93 |
| A1 | euclidean | 32 | 0.8 | 17 | 0.55 | 0.87 | 0.88 | 0.89 |
| A1 | euclidean | 32 | 1.0 | 19 | 0.59 | 0.69 | 0.83 | 0.87 |
| A1 | euclidean | 32 | 1.5 | 25 | 0.55 | 0.54 | 0.78 | 0.80 |
| C_gate | euclidean | 32 | 0.1 | 4 | 0.23 | 0.58 | 0.74 | 0.94 |
| C_gate | euclidean | 32 | 0.2 | 5 | 0.28 | 0.67 | 0.80 | 0.96 |
| C_gate | euclidean | 32 | 0.3 | 5 | 0.28 | 0.67 | 0.80 | 0.82 |
| C_gate | euclidean | 32 | 0.5 | 8 | 0.50 | 0.91 | 0.89 | 0.96 |
| C_gate | euclidean | 32 | 0.8 | 10 | 0.56 | 0.87 | 0.87 | 0.91 |
| C_gate | euclidean | 32 | 1.0 | 12 | 0.52 | 0.73 | 0.82 | 0.73 |
| C_gate | euclidean | 32 | 1.5 | 20 | 0.47 | 0.45 | 0.72 | 0.64 |
| C_gate | hamming | 25 | 0.1 | 18 | 0.17 | 0.48 | 0.63 | 0.97 |
| C_gate | hamming | 25 | 0.2 | 22 | 0.30 | 0.58 | 0.70 | 0.75 |
| C_gate | hamming | 25 | 0.3 | 30 | 0.34 | 0.66 | 0.72 | 0.80 |
| C_gate | hamming | 25 | 0.5 | 37 | 0.24 | 0.31 | 0.59 | 0.66 |
| C_gate | hamming | 25 | 0.8 | 48 | 0.28 | 0.24 | 0.57 | 0.59 |
| C_gate | hamming | 25 | 1.0 | 57 | 0.31 | 0.19 | 0.57 | 0.58 |
| C_gate | hamming | 25 | 1.5 | 87 | 0.28 | 0.13 | 0.53 | 0.61 |
| D_gate | hamming | 32 | 0.1 | 3 | 0.17 | 0.52 | 0.68 | 0.91 |
| D_gate | hamming | 32 | 0.2 | 5 | 0.28 | 0.67 | 0.79 | 0.98 |
| D_gate | hamming | 32 | 0.3 | 6 | 0.37 | 0.88 | 0.86 | 0.95 |
| D_gate | hamming | 32 | 0.5 | 11 | 0.43 | 0.90 | 0.87 | 0.96 |
| D_gate | hamming | 32 | 0.8 | 10 | 0.42 | 0.79 | 0.83 | 0.81 |
| D_gate | hamming | 32 | 1.0 | 16 | 0.43 | 0.58 | 0.75 | 0.68 |
| D_gate | hamming | 32 | 1.5 | 25 | 0.38 | 0.36 | 0.68 | 0.61 |
| D_gate_clean | hamming | 15 | 0.1 | 393 | 0.38 | 0.07 | 0.47 | 0.68 |
| D_gate_clean | hamming | 15 | 0.2 | 432 | 0.36 | 0.06 | 0.47 | 0.75 |
| D_gate_clean | hamming | 15 | 0.3 | 460 | 0.39 | 0.05 | 0.47 | 0.77 |
| D_gate_clean | hamming | 15 | 0.5 | 507 | 0.38 | 0.04 | 0.46 | 0.81 |
| D_gate_clean | hamming | 15 | 0.8 | 546 | 0.39 | 0.04 | 0.47 | 0.85 |
| D_gate_clean | hamming | 15 | 1.0 | 572 | 0.39 | 0.03 | 0.46 | 0.85 |
| D_gate_clean | hamming | 15 | 1.5 | 618 | 0.38 | 0.03 | 0.46 | 0.87 |

Leiden on a kNN graph is the wrong clusterer for binary data, and the table shows why. On
D_gate_clean, 50,000 cells share 3,286 phenotypes, so most cells have hundreds of neighbors at
distance zero and the graph breaks into islands. Leiden returns 393 to 618 clusters at every
resolution. On D_gate the 17 noise bits break the ties, the graph connects, and at 11 clusters
the result is mean F1 0.43 and ARI 0.90 against 0.51 and 0.97 for A1 at the same cluster count.
The Hamming graph on the 25 lineage bits of C_gate sits between the two cases and loses at every
resolution. The clusterer that suits binary data is a count of phenotypes with a merge rule,
which is flowType, and it is not a kNN graph.

## Does binary input beat A1 when the thresholds are as good as they can be?

No. With gate thresholds on all 32 markers, the binary arms score the same as they did with
mixture thresholds at 40 metaclusters, 0.45 to 0.48 against 0.79, because the label set cannot
define 17 of the 32 markers and those bits are noise under any cut. With the 15 markers it can
define, binary input ties A1 at 14 metaclusters on mean F1 (0.47 against 0.47) and beats it on
ARI (0.78 against 0.66), and at 40 metaclusters it reaches 0.68 against 0.79. That arm used the
labels to choose the markers and to place the cuts, and it still lost at the free cluster count
and on the three smallest populations. Leiden on a Hamming graph did not change the order. The
honest reading is that hypothesis 2 does not hold on Levine_32dim under any threshold the
labels can supply. Two things survive. The phenotype count behaves as the plan predicted once
the bits are clean (3,286 occupied, 371 for 95 percent of cells), which makes it a usable check
on threshold quality. And the binary arm's failures are the hypothesis 3 cases, named per
population above, so the loss has a mechanism and not only a number.

## Caveats

- The truth table is a judgment made in this session. A reviewer with the Levine 2015 gating
  figures may move a population between the lists.
- One seed for every full-data FlowSOM result. The stability column says how far one seed moves
  the result.
- The mixture arms (C, C2, D) do not reproduce the pilot to the third decimal, because the
  mixture fits on a subsample drawn from the module's random state and that state differs when
  the gate arms run in the same job. The gate arms have no such dependence.
- D_gate_clean has 15 columns, so its PCA and two-NN numbers are not comparable to the 32-column
  arms.
- Levine_13dim, Samusik and the conventional flow sets did not run with gate thresholds. The
  conventional panels have fewer dim channels, so the all-marker binary arm has less noise to
  lose to there. That is a reason to run them, not a reason to expect a different verdict.
- Nothing here is in git.
