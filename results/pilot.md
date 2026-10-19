# Pilot results — four representations, FlowSOM, Weber-style scoring

**Run:** 2026-10-18 · **Code:** `code/pilot.py`, one command in `code/run_pilot.sh` · **Data:** see `data/README.md`

## What ran

Three datasets, one clustering method (FlowSOM, 10 by 10 grid, consensus metaclustering), one
seed for the full-data result and 10 seeds on a 50,000-cell subsample for stability.

| Dataset | Cells | Labeled | Markers | Metaclusters | Ground truth |
|---|---|---|---|---|---|
| Levine_32dim | 265,627 | 104,184 | 32 | 14 (the gated number) and 40 (free) | 14 populations, human bone marrow, mass cytometry |
| Nilsson_rare | 44,140 | all | 13 | 40 | 358 HSC, 0.8 percent, conventional flow |
| Mosmann_rare | 396,460 | all | 14 (7 lineage, 7 cytokine) | 40 | 109 activated memory CD4 T cells, 0.03 percent, conventional flow |

The plan defines the arms. Binarization used a two-component Gaussian mixture per marker
on arcsinh values, cut at equal posterior. Arm B ran only on Nilsson_rare, the one set whose
channel names carry fluorochromes. Mosmann_rare names its channels by filter, so its arm B waits
for the panel table in the 2014 paper.

## Hypothesis 1: brightness changes cluster membership

Supported as a statement about membership. Dividing by the conjugate brightness before the
transform moved cells between clusters. It did not improve recovery of the rare population.

| Pair | ARI |
|---|---|
| A1 vs B_geom | 0.69 |
| A1 vs B_lin | 0.59 |
| A2 vs B_geom | 0.53 |
| A1 vs A2 | 0.55 |
| A1 vs A3 | 0.58 |
| B_geom vs B_lin | 0.57 |

Agreement between arm B and the z-score control A2 is 0.53, so arm B is not the same operation as
per-channel scaling. The divisors are the weak, vendor-table source, and I assumed the four Qdot values
because the fluor-ai table has no Qdot entries.

| Marker | Fluorochrome | Divisor (2^rank, geometric mean 1) |
|---|---|---|
| CD38 | FITC | 0.59 |
| CD34 | PerCP-Cy5.5 | 0.59 |
| CD45 | Alexa Fluor 700 | 0.59 |
| CD3 | APC-Cy7 | 0.59 |
| CD45RA | Pacific Blue | 0.59 |
| CD49fpur | Qdot 705 | 0.59 |
| CD90bio | Qdot 800 | 0.59 |
| CD11b | Qdot 605 | 1.17 |
| CD4 | Qdot 655 | 1.17 |
| CD123 | PE | 2.35 |
| CD19 | PE-Cy5 | 2.35 |
| CD10 | PE-Cy7 | 2.35 |
| CD110 | APC | 2.35 |

The same effect appears without any fluorochrome. On Levine_32dim at 40 metaclusters, z-scoring
the channels dropped mean F1 from 0.79 to 0.63 and the ARI between the two clusterings is 0.84 at
14 metaclusters. Per-channel reweighting is not a neutral step on mass cytometry either.

## Hypothesis 2: binary lineage markers recover the gated populations

Not supported in this pilot. With unsupervised thresholds the binary arms lost to the unscaled
intensity arm on every dataset once the number of metaclusters was free.

Levine_32dim, 40 metaclusters:

| Arm | Mean F1 | ARI | NMI | Stability ARI |
|---|---|---|---|---|
| A1 arcsinh | 0.79 | 0.80 | 0.83 | 0.77 |
| A2 z-score | 0.63 | 0.68 | 0.76 | 0.70 |
| A3 percentile | 0.64 | 0.52 | 0.74 | 0.64 |
| C lineage binary + percentile | 0.40 | 0.22 | 0.56 | 0.39 |
| C2 lineage binary + z-score | 0.46 | 0.29 | 0.56 | 0.36 |
| D all binary | 0.41 | 0.24 | 0.58 | 0.40 |

Levine_32dim, 14 metaclusters (forced to the gated number):

| Arm | Mean F1 | ARI | NMI | Stability ARI |
|---|---|---|---|---|
| A1 arcsinh | 0.47 | 0.66 | 0.83 | 0.94 |
| A2 z-score | 0.43 | 0.63 | 0.81 | 0.90 |
| A3 percentile | 0.41 | 0.64 | 0.82 | 0.84 |
| C lineage binary + percentile | 0.40 | 0.52 | 0.66 | 0.55 |
| C2 lineage binary + z-score | 0.46 | 0.62 | 0.71 | 0.71 |
| D all binary | 0.49 | 0.63 | 0.71 | 0.57 |

Nilsson_rare, best cluster for the HSC population:

| Arm | F1 | Precision | Recall | Stability ARI |
|---|---|---|---|---|
| A1 arcsinh | 0.50 | 0.36 | 0.84 | 0.54 |
| A2 z-score | 0.38 | 0.25 | 0.87 | 0.50 |
| A3 percentile | 0.23 | 0.13 | 0.94 | 0.52 |
| B brightness, 2^rank | 0.22 | 0.12 | 0.93 | 0.61 |
| B brightness, linear rank | 0.20 | 0.11 | 0.95 | 0.57 |
| C lineage binary + percentile | 0.41 | 0.29 | 0.72 | 0.79 |
| C2 lineage binary + z-score | 0.40 | 0.44 | 0.37 | 0.76 |
| D all binary | 0.33 | 0.22 | 0.65 | 0.77 |

Mosmann_rare, best cluster for the activated memory CD4 T cells:

| Arm | F1 | Precision | Recall | Stability ARI |
|---|---|---|---|---|
| A1 arcsinh | 0.57 | 0.40 | 0.98 | 0.56 |
| A2 z-score | 0.61 | 0.44 | 0.99 | 0.59 |
| A3 percentile | 0.09 | 0.05 | 0.72 | 0.53 |
| C lineage binary + percentile | 0.10 | 0.06 | 0.67 | 0.96 |
| C2 lineage binary + z-score | 0.55 | 0.38 | 0.96 | 0.66 |
| D all binary | 0.07 | 0.04 | 0.32 | 0.57 |

Three things inside those tables matter more than the headline:

1. **The unscaled arm A1 is the best or near the best everywhere.** Its mean F1 of 0.79 on
   Levine_32dim at 40 metaclusters matches the FlowSOM result in Weber and Robinson, which says
   the pipeline reproduces the benchmark.
2. **Percentile scaling destroys a rare population.** On Mosmann_rare the 0.03 percent
   cytokine-positive cells sit above the 99th percentile of every cytokine channel, so clipping
   puts them in the same bin as the top 1 percent of negatives. F1 fell from 0.57 to 0.09. Arm C
   inherited the fault through its continuous columns, and arm C2, which z-scores them instead,
   recovered to 0.55.
3. **Stability went both ways.** The binary arms were more stable on Nilsson_rare (0.79 against
   0.54) and less stable on Levine_32dim (0.55 against 0.94). On 32 mass cytometry markers the
   Gaussian mixture produced 148,505 distinct phenotypes among 265,627 cells, so the binary data
   is nearly unique per cell and the SOM has nothing to hold on to.

The binary arms did win in one place. At 14 metaclusters every intensity arm merged CD8 T cells
into CD4 T cells (F1 0.01) and arm D kept them apart (F1 0.77), and arm D was the only arm to find
basophils (0.62). That is a merge order effect, not a recovery effect, and it disappears at 40.

## Hypothesis 3: binarization fails on intermediate-expression populations

Supported, and the loss is large. Levine_32dim at 40 metaclusters, F1 per population:

| Population | n | A1 | A2 | A3 | C | C2 | D |
|---|---|---|---|---|---|---|---|
| CD4 T cells | 26,366 | 0.96 | 0.91 | 0.81 | 0.32 | 0.51 | 0.32 |
| Monocytes | 21,099 | 0.66 | 0.59 | 0.38 | 0.37 | 0.58 | 0.32 |
| CD8 T cells | 20,108 | 0.96 | 0.95 | 0.82 | 0.38 | 0.39 | 0.44 |
| Mature B cells | 16,520 | 0.96 | 0.63 | 0.56 | 0.54 | 0.68 | 0.58 |
| Pre B cells | 6,135 | 0.91 | 0.87 | 0.86 | 0.51 | 0.40 | 0.57 |
| CD16- NK cells | 3,905 | 0.72 | 0.94 | 0.86 | 0.59 | 0.21 | 0.41 |
| CD34+CD38+CD123- HSPCs | 3,295 | 0.63 | 0.38 | 0.39 | 0.48 | 0.39 | 0.71 |
| CD16+ NK cells | 2,248 | 0.88 | 0.76 | 0.90 | 0.48 | 0.55 | 0.58 |
| pDCs | 1,238 | 0.81 | 0.69 | 0.89 | 0.65 | 0.56 | 0.72 |
| Basophils | 1,207 | 0.90 | 0.42 | 0.90 | 0.65 | 0.64 | 0.72 |
| CD34+CD38lo HSCs | 916 | 0.74 | 0.43 | 0.51 | 0.34 | 0.45 | 0.09 |
| Pro B cells | 513 | 0.66 | 0.35 | 0.54 | 0.16 | 0.22 | 0.14 |
| Plasma B cells | 330 | 0.79 | 0.78 | 0.01 | 0.01 | 0.68 | 0.01 |
| CD34+CD38+CD123+ HSPCs | 304 | 0.55 | 0.09 | 0.48 | 0.11 | 0.12 | 0.16 |

CD34+CD38lo HSCs go from 0.74 under A1 to 0.09 under D. CD11b-defined monocyte subsets are not in
this label set, so the dim and mid cases named in the plan still need Levine_13dim, which sits in
the mirror and has not run.

## Dimensionality

| Arm | PC1 variance | PC1 to PC5 variance | Two-NN intrinsic dimension | Occupied phenotypes | Possible phenotypes | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|
| A1 arcsinh | 0.37 | 0.68 | 19.11 |  |  |  |
| A2 z-score | 0.17 | 0.49 | 19.60 |  |  |  |
| A3 percentile | 0.23 | 0.56 | 18.84 |  |  |  |
| C lineage binary + percentile | 0.17 | 0.44 | 11.22 | 112988.00 | 33554432.00 | 99707.00 |
| C2 lineage binary + z-score | 0.20 | 0.57 | 17.49 | 112988.00 | 33554432.00 | 99707.00 |
| D all binary | 0.17 | 0.42 | 14.83 | 148395.00 | 4294967296.00 | 135114.00 |

Nilsson_rare:

| Arm | PC1 variance | PC1 to PC5 variance | Two-NN intrinsic dimension | Occupied phenotypes | Possible phenotypes | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|
| A1 arcsinh | 0.52 | 0.92 | 11.63 |  |  |  |
| A2 z-score | 0.39 | 0.88 | 11.57 |  |  |  |
| A3 percentile | 0.46 | 0.90 | 11.95 |  |  |  |
| B brightness, 2^rank | 0.49 | 0.92 | 11.61 |  |  |  |
| B brightness, linear rank | 0.51 | 0.92 | 11.52 |  |  |  |
| C lineage binary + percentile | 0.39 | 0.82 | 1.41 | 596.00 | 4096.00 | 131.00 |
| C2 lineage binary + z-score | 0.48 | 0.84 | 1.37 | 596.00 | 4096.00 | 131.00 |
| D all binary | 0.41 | 0.81 | 14.24 | 851.00 | 8192.00 | 200.00 |

On 13 conventional markers the picture the plan predicted holds: the cells occupy 851 of 8,192
phenotypes, and 200 of them hold 95 percent of the cells. On 32 mass cytometry markers it does not
hold with these thresholds, and the count of occupied phenotypes is the diagnostic that says the
thresholds are noisy.

## Figures

- `fig1_brightness.png` hypothesis 1: the divisors and the agreement between arms
- `fig2_f1.png` hypothesis 2: recovery and stability per arm
- `fig3_failure.png` hypothesis 3: per-population F1 on Levine_32dim at 40 metaclusters
- `fig4_dimensionality.png` occupied phenotypes against 2^k, intrinsic dimension per arm

## Where this leaves the paper

The plan's stop rule said: stop if arm C loses to arm A on every metric and arm B equals A2
everywhere. Arm B does not equal A2, and arm C won on stability for one dataset and on two merged
populations at the forced cluster number, so the rule does not fire. The outlook for hypothesis 2
is still poor. With unsupervised thresholds, binary input lost to plain arcsinh intensity on every
dataset at a free cluster number, and it lost most on the populations that intermediate expression
defines.

Two results from this pilot stand on their own and do not depend on hypothesis 2:

- Per-channel scaling is not neutral. Z-scoring cost 0.16 of mean F1 on Levine_32dim, and
  percentile clipping erased a rare population on Mosmann_rare. Many published pipelines do one of
  these by default.
- Dividing by brightness changes cluster membership on the one conventional flow set that carries
  its fluorochromes, and the change is not the same as scaling.

One test decides the next step. Run the binary arms with the manual gate thresholds, which is the
upper bound the plan names, on Levine_32dim at 40 metaclusters. If binary input still loses with
perfect thresholds, hypothesis 2 ends and the paper is about hypotheses 1 and 3 and the
scaling result. If binary input wins with perfect thresholds, the paper is about thresholds, and
flowType's multi-level thresholds become the method to test.

**That test ran on 2026-10-18 and binary input lost.** The full result is in
`results/gate_thresholds.md`. With gate thresholds on all 32 markers, the binary arms stayed at
mean F1 0.45 to 0.48 against 0.79 at 40 metaclusters, because the label set cannot define 17 of
the markers and those bits are noise under any cut. With only the 15 markers the labels define
cleanly, binary input tied A1 at 14 metaclusters and reached 0.68 against 0.79 at 40, and it
lost the three smallest populations. Leiden on a Hamming graph did not change the order.
Hypothesis 2 ends on this dataset. The paper is about hypotheses 1 and 3 and the scaling
result, and the phenotype count stays as a check on threshold quality.

## Caveats

- One seed for every full-data result. The stability numbers show that one seed moves F1.
- The Gaussian mixture threshold is a poor fit for mass cytometry channels with a zero spike.
- The two-NN intrinsic dimension of 1.4 for arms C and C2 on Nilsson_rare is an artifact. With
  12 binary columns and one continuous column, the nearest neighbors of a cell share its phenotype
  and differ only in CD38, so the estimator sees one dimension. Read that number only for the
  continuous arms.
- Arm B used ordinal vendor ranks, four of them assumed, mapped to divisors by two provisional
  rules. Measured stain index from single stains is the next source up.
- Nothing here is in git.
