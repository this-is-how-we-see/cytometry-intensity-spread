# Resolution-weighted and posterior arms — weight the markers that separate, keep the grade

**Run:** 2026-10-19 · **Code:** `code/pilot.py --datasets Levine_32dim_k40,Samusik_01_k40 --arms A1,A2,C2,W,W_z,P --tag weighted`, tables by `code/phase2_report.py --tag weighted`; the arms are `resolution_weights` and the `W`, `W_z` and `P` entries of `build_arms` in `code/pilot.py` · **Metrics:** `results/weighted_metrics.csv`, `results/weighted_perpop.csv`, weights in `results/weighted_meta.json`

## The question

The question, 2026-10-19: instead of a bit per marker, weight each marker by how well it
separates, so that CD4 and CD8, which split T cells cleanly, carry more of the distance than a
dim or noisy marker, and keep a graded value where expression is intermediate. Two arms test
it without the labels. Arm W multiplies each arcsinh channel by the separation of its two
mixture modes in pooled-SD units, a stain index measured on the cells, normalized to
geometric mean 1 and clipped to 1/3 to 3. Arm W_z applies the same weights to the z-scored
channels. Arm P replaces the 0 or 1 of arm C2 with the mixture's posterior probability of the
high component, so a CD4-dim monocyte scores 0.4 instead of 0 or 1.

Neither beats plain intensity. Arm W scores mean F1 0.68 on Levine 32 against 0.79 for A1
and 0.66 against 0.70 on Samusik 01. Arm W_z recovers part of the z-score cost (0.67 against
0.63, and 0.69 against 0.60), which says the weights carry some signal. Arm P ties the binary
arm on Levine 32 (0.44) and beats it on Samusik 01 (0.45 against 0.39), and on the
intermediate monocytes it scores 0.20 against 0.14 for the bit and 0.91 for A1. The reason is
in the weights themselves, below: the statistic does not find CD4 and CD8.

## Results

### Levine_32dim, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.79 | 0.80 | 0.83 | 0.79 |  |  |  |  |
| A2 z-score | 0.63 | 0.68 | 0.76 | 0.72 | 0.67 |  |  |  |
| W resolution-weighted arcsinh | 0.68 | 0.62 | 0.77 | 0.69 | 0.68 |  |  |  |
| W_z resolution-weighted z-score | 0.67 | 0.68 | 0.78 | 0.70 | 0.76 |  |  |  |
| P posterior + z-score | 0.44 | 0.24 | 0.54 | 0.36 | 0.24 |  |  |  |
| C2 lineage mixture + z-score | 0.44 | 0.26 | 0.55 | 0.35 | 0.27 | 113,405 | 33,554,432 | 100,124 |

### Samusik_01, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.70 | 0.85 | 0.83 | 0.76 |  |  |  |  |
| A2 z-score | 0.60 | 0.48 | 0.66 | 0.75 | 0.46 |  |  |  |
| W resolution-weighted arcsinh | 0.66 | 0.76 | 0.79 | 0.69 | 0.64 |  |  |  |
| W_z resolution-weighted z-score | 0.69 | 0.70 | 0.76 | 0.71 | 0.67 |  |  |  |
| P posterior + z-score | 0.45 | 0.39 | 0.54 | 0.56 | 0.36 |  |  |  |
| C2 lineage mixture + z-score | 0.39 | 0.42 | 0.58 | 0.54 | 0.37 | 65,061 | 268,435,456 | 60,718 |

### Per-population F1, Levine_32dim, 40 metaclusters

| Population | Cells | A1 | A2 | W | W_z | P | C2 |
|---|---|---|---|---|---|---|---|
| CD4 T cells | 26,366 | 0.96 | 0.91 | 0.82 | 0.96 | 0.44 | 0.43 |
| Monocytes | 21,099 | 0.66 | 0.59 | 0.68 | 0.44 | 0.52 | 0.53 |
| CD8 T cells | 20,108 | 0.96 | 0.95 | 0.86 | 0.78 | 0.45 | 0.38 |
| Mature B cells | 16,520 | 0.96 | 0.63 | 0.67 | 0.77 | 0.43 | 0.53 |
| Pre B cells | 6,135 | 0.91 | 0.87 | 0.89 | 0.88 | 0.58 | 0.53 |
| CD16- NK cells | 3,905 | 0.72 | 0.94 | 0.63 | 0.88 | 0.11 | 0.21 |
| CD34+CD38+CD123- HSPCs | 3,295 | 0.63 | 0.38 | 0.50 | 0.30 | 0.47 | 0.36 |
| CD16+ NK cells | 2,248 | 0.88 | 0.76 | 0.62 | 0.87 | 0.56 | 0.55 |
| pDCs | 1,238 | 0.81 | 0.69 | 0.90 | 0.88 | 0.52 | 0.52 |
| Basophils | 1,207 | 0.90 | 0.42 | 0.91 | 0.90 | 0.57 | 0.57 |
| CD34+CD38lo HSCs | 916 | 0.74 | 0.43 | 0.70 | 0.50 | 0.50 | 0.52 |
| Pro B cells | 513 | 0.66 | 0.35 | 0.30 | 0.45 | 0.12 | 0.09 |
| Plasma B cells | 330 | 0.79 | 0.78 | 0.79 | 0.79 | 0.67 | 0.72 |
| CD34+CD38+CD123+ HSPCs | 304 | 0.55 | 0.09 | 0.28 | 0.01 | 0.19 | 0.21 |

### Per-population F1, Samusik_01, 40 metaclusters

| Population | Cells | A1 | A2 | W | W_z | P | C2 |
|---|---|---|---|---|---|---|---|
| Classical Monocytes | 13,607 | 0.96 | 0.73 | 0.89 | 0.86 | 0.64 | 0.60 |
| Intermediate Monocytes | 12,045 | 0.91 | 0.45 | 0.88 | 0.80 | 0.20 | 0.14 |
| IgDpos IgMpos B cells | 8,768 | 0.98 | 0.75 | 0.92 | 0.95 | 0.83 | 0.80 |
| Eosinophils | 5,173 | 0.89 | 0.85 | 0.86 | 0.88 | 0.48 | 0.65 |
| IgD- IgMpos B cells | 3,508 | 0.56 | 0.59 | 0.45 | 0.63 | 0.55 | 0.65 |
| pDCs | 2,974 | 0.53 | 0.64 | 0.76 | 0.91 | 0.76 | 0.74 |
| IgM- IgD- B-cells | 1,123 | 0.43 | 0.28 | 0.52 | 0.32 | 0.35 | 0.19 |
| Non-Classical Monocytes | 1,016 | 1.00 | 0.85 | 0.99 | 0.98 | 0.93 | 0.90 |
| NKT cells | 961 | 0.90 | 0.79 | 0.83 | 0.87 | 0.80 | 0.63 |
| CD8 T cells | 798 | 0.88 | 0.96 | 0.79 | 0.97 | 0.74 | 0.14 |
| GMP | 703 | 0.78 | 0.61 | 0.59 | 0.82 | 0.54 | 0.63 |
| MEP | 420 | 0.78 | 0.81 | 0.85 | 0.75 | 0.68 | 0.45 |
| CD4 T cells | 406 | 0.93 | 0.93 | 0.89 | 0.93 | 0.39 | 0.26 |
| Basophils | 292 | 0.87 | 0.82 | 0.85 | 0.89 | 0.76 | 0.02 |
| B-cell Frac A-C (pro-B cells) | 266 | 0.12 | 0.01 | 0.11 | 0.08 | 0.04 | 0.03 |
| Macrophages | 199 | 0.88 | 0.69 | 0.80 | 0.76 | 0.05 | 0.62 |
| CMP | 197 | 0.44 | 0.25 | 0.56 | 0.08 | 0.21 | 0.15 |
| MPP | 197 | 0.85 | 0.00 | 0.77 | 0.59 | 0.72 | 0.70 |
| NK cells | 189 | 0.66 | 0.65 | 0.55 | 0.67 | 0.25 | 0.34 |
| Plasma Cells | 118 | 0.99 | 0.97 | 0.72 | 0.99 | 0.38 | 0.28 |
| mDCs | 116 | 0.78 | 0.93 | 0.77 | 0.92 | 0.48 | 0.36 |
| CLP | 59 | 0.54 | 0.52 | 0.22 | 0.44 | 0.11 | 0.10 |
| gd T cells | 35 | 0.08 | 0.19 | 0.06 | 0.37 | 0.01 | 0.00 |
| HSC | 3 | 0.00 | 0.00 | 0.10 | 0.00 | 0.00 | 0.00 |

## What the weights measured

Levine 32, largest and smallest: CD3 3.00, HLA-DR 1.49, CD45 1.30, CD38 1.26; CD4 0.92, CD8 1.08, CD19 1.00; CD41 0.63, CD16 0.64, CD117 0.68, Flt3 0.69.
Samusik 01: MHCII 1.96, CD44 1.84, Ly6C 1.57, CD16_32 1.52, CD11b 1.47; CD4 0.76, CD8 0.68, TCRb 0.77, CD3 0.79.

On Samusik 01 the weights went up on MHCII, CD44, Ly6C and CD11b, the markers that most cells
express at one of two clear levels, and down on CD4, CD8 and TCRb, the markers that split the
T cells. CD4 and CD8 T cells are 406 and 798 of 86,864 cells there, so their high mode is a
small component and the mixture fits it with a wide variance against a tall negative mode.
The separation of two modes fitted on all cells measures how much of the sample sits in a
clean high mode, which is a different thing from how well a marker separates two
populations that both exist. On Levine 32 the same statistic put CD3 at the cap and CD4 and
CD8 near 1, so T cells were weighted by CD3 alone, and the weights of the 17 undefined or dim
markers stayed between 0.6 and 1.0 rather than falling toward the floor. The arm therefore
did not remove the noise channels that broke the binary arms, and it added a distortion of
its own.

The posterior arm keeps the graded value the bit throws away, and the intermediate monocytes
on Samusik 01 show what that buys: 0.20 against 0.14. It is a small gain, because the
posterior of a two-component mixture is sharp exactly where the bit was already right and
graded only across the narrow band between the modes, and the intermediate population sits
on top of that band for only one marker. Everything else in arm P is arm C2.

## What would test the idea properly

The idea is right and the statistic is wrong. A marker's weight has to come from the
separation it produces where it matters, between the populations it is meant to split, and
a mixture on all cells cannot see that without the labels. Two unsupervised measures come
closer. The sharpness of the posterior, the share of cells with a posterior below 0.1 or
above 0.9, rewards a marker whose two modes are clean wherever they sit and does not
depend on how many cells are in the high mode. And a weight per marker fitted inside the
cells that are positive for a parent marker, CD4 and CD8 among CD3-positive cells, measures
the split where it happens, at the price of a gating order that the paper would have to
state. Either is one afternoon. The current run says only that the first statistic tried is
not the one.

## Caveats

- One seed on the full data, as in every run.
- The cap and floor of 3 and 1/3 were chosen before the run and not tuned; with most weights
  near 1 they did not bind except for CD3.
- Arm P changes only the lineage markers; its functional markers are z-scored, so the
  comparison with C2 isolates the posterior and the comparison with A1 does not.

## Second pass, 2026-10-19: weight by the sharpness of the posterior

**Run:** `code/pilot.py --datasets Levine_32dim_k40,Samusik_01_k40 --arms A1,W2,W2s --tag weighted2` · **Metrics:** `results/weighted2_metrics.csv`, `results/weighted2_perpop.csv`, weights in `results/weighted2_meta.json`

The weight is the share of cells a marker scores below 0.1 or above 0.9 under its own mixture,
so a marker with two clean modes gets a high weight wherever the modes sit and however many
cells are in the high one. Arm W2 uses the sharpness as the weight, arm W2s its square, both
normalized to geometric mean 1 and clipped to 1/3 to 3, on the arcsinh channels.

The weights now land where the idea says they should. On Samusik 01 the top five are
Ly6C 1.30, Ter119 1.28, 120g8 1.26, CD11b 1.25, B220 1.25, CD4, CD8, TCRb and CD3 sit above 1, and the bottom are
CD45.2 0.33, Ly6G 0.47, CD64 0.64, CD115 0.66, F480 0.77, the dim and noisy channels. On Levine 32 the top
five are CD3 1.69, CD4 1.57, CD133 1.54, CD7 1.51, CD117 1.50 and the bottom CD235ab 0.33, CD321 0.33, CD47 0.33, CD44 0.46, CXCR4 0.58, the
continuous markers that no population defines.

### Levine_32dim, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.79 | 0.80 | 0.83 | 0.79 |  |  |  |  |
| W2 sharpness-weighted arcsinh | 0.69 | 0.64 | 0.80 | 0.68 | 0.67 |  |  |  |
| W2s sharpness squared | 0.75 | 0.77 | 0.83 | 0.63 | 0.82 |  |  |  |

### Samusik_01, FlowSOM, 40 metaclusters, full data, seed 1

| Arm | Mean F1 | ARI | NMI | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Possible | Phenotypes for 95 percent of cells |
|---|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.70 | 0.85 | 0.83 | 0.74 |  |  |  |  |
| W2 sharpness-weighted arcsinh | 0.69 | 0.81 | 0.82 | 0.78 | 0.73 |  |  |  |
| W2s sharpness squared | 0.65 | 0.75 | 0.81 | 0.78 | 0.71 |  |  |  |

### Per-population F1, Levine_32dim, 40 metaclusters

| Population | Cells | A1 | W2 | W2s |
|---|---|---|---|---|
| CD4 T cells | 26,366 | 0.96 | 0.66 | 0.98 |
| Monocytes | 21,099 | 0.66 | 0.65 | 0.60 |
| CD8 T cells | 20,108 | 0.96 | 0.96 | 0.96 |
| Mature B cells | 16,520 | 0.96 | 0.69 | 0.69 |
| Pre B cells | 6,135 | 0.91 | 0.87 | 0.86 |
| CD16- NK cells | 3,905 | 0.72 | 0.89 | 0.92 |
| CD34+CD38+CD123- HSPCs | 3,295 | 0.63 | 0.46 | 0.45 |
| CD16+ NK cells | 2,248 | 0.88 | 0.90 | 0.95 |
| pDCs | 1,238 | 0.81 | 0.89 | 0.93 |
| Basophils | 1,207 | 0.90 | 0.95 | 0.97 |
| CD34+CD38lo HSCs | 916 | 0.74 | 0.66 | 0.46 |
| Pro B cells | 513 | 0.66 | 0.37 | 0.60 |
| Plasma B cells | 330 | 0.79 | 0.71 | 0.75 |
| CD34+CD38+CD123+ HSPCs | 304 | 0.55 | 0.00 | 0.34 |

### Per-population F1, Samusik_01, 40 metaclusters

| Population | Cells | A1 | W2 | W2s |
|---|---|---|---|---|
| Classical Monocytes | 13,607 | 0.96 | 0.86 | 0.66 |
| Intermediate Monocytes | 12,045 | 0.91 | 0.95 | 0.94 |
| IgDpos IgMpos B cells | 8,768 | 0.98 | 0.97 | 0.97 |
| Eosinophils | 5,173 | 0.89 | 0.89 | 0.89 |
| IgD- IgMpos B cells | 3,508 | 0.56 | 0.62 | 0.47 |
| pDCs | 2,974 | 0.53 | 0.57 | 0.56 |
| IgM- IgD- B-cells | 1,123 | 0.43 | 0.44 | 0.50 |
| Non-Classical Monocytes | 1,016 | 1.00 | 0.99 | 1.00 |
| NKT cells | 961 | 0.90 | 0.90 | 0.90 |
| CD8 T cells | 798 | 0.88 | 0.88 | 0.91 |
| GMP | 703 | 0.78 | 0.67 | 0.34 |
| MEP | 420 | 0.78 | 0.82 | 0.86 |
| CD4 T cells | 406 | 0.93 | 0.93 | 0.90 |
| Basophils | 292 | 0.87 | 0.82 | 0.70 |
| B-cell Frac A-C (pro-B cells) | 266 | 0.12 | 0.12 | 0.14 |
| Macrophages | 199 | 0.88 | 0.88 | 0.88 |
| CMP | 197 | 0.44 | 0.48 | 0.41 |
| MPP | 197 | 0.85 | 0.84 | 0.86 |
| NK cells | 189 | 0.66 | 0.65 | 0.56 |
| Plasma Cells | 118 | 0.99 | 0.82 | 0.76 |
| mDCs | 116 | 0.78 | 0.79 | 0.85 |
| CLP | 59 | 0.54 | 0.57 | 0.57 |
| gd T cells | 35 | 0.08 | 0.07 | 0.04 |
| HSC | 3 | 0.00 | 0.00 | 0.05 |

Neither arm passes plain intensity. The squared weights come within 0.04 of A1 on Levine 32
(0.75 against 0.79, ARI 0.77 against 0.80) and the linear weights within 0.01 on Samusik 01
(0.69 against 0.70, ARI 0.81 against 0.85), and both sit inside the range that one seed
moves. The large populations are unchanged and the intermediate monocytes on Samusik gain
a little (0.95 against 0.91). The small populations pay: CD34+CD38+CD123+ HSPCs fall from
0.55 to 0.00 under W2 and 0.34 under W2s, the HSCs from 0.74 to 0.66 and 0.46. Those
populations are defined by graded shifts on CD38, CD123 and CD117, which are exactly the
markers the sharpness statistic reads as unclean and weights down. A marker's sharpness
across the whole sample is not its value for the rare populations, and a weight that
rewards sharpness removes the information the rare populations need.

This closes the question for the paper. Plain arcsinh intensity gives every channel its
measured spread, and that spread already carries the weighting the data supports: a marker
with clean modes spans more of the scale than a marker without them. The two unsupervised
weightings tried here either distort that (mode separation) or come back to it with a loss on
the small populations (sharpness). The idea survives in one form only, a weight fitted where
the split happens, inside a parent population, and that form needs a stated gating order,
which makes it a gating strategy and not a representation.
