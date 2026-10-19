# Mass cytometry against conventional flow on shared markers — the no-dye control

**Run:** 2026-10-19 · **Code:** `code/pilot.py --datasets Nilsson_shared,Levine_13dim_shared --arms A1,A2,A3,B_geom,B_lin,C2,D --tag massflow`, then `--arms A1,B_meas --tag massflow_meas` on Levine 13 · **Metrics:** `results/massflow_metrics.csv`, `results/massflow_meas_metrics.csv` · **Loaders:** `load_levine13_shared`, `load_nilsson_shared` in `code/pilot.py`

## The question

Arcsinh intensity wins on every labeled set. Mass cytometry has no fluorochrome, so a mass set
is the nearest thing to a control where the dye term is absent. If arm B removes the dye term
from a flow set, the flow set under arm B behaves more like the mass set than the flow set
under A1 does. Can the two platforms be compared on the same markers?

They can, on one pair, and the pair answers a narrower question than the one asked. Levine 13
(CyTOF, human bone marrow) and Nilsson rare (conventional flow, human bone marrow) share 10
markers: CD45, CD45RA, CD19, CD11b, CD4, CD34, CD123, CD38, CD90 and CD3. Both label the HSC.
On those 10 markers the flow set recovers its HSC at F1 0.50 under plain intensity and at 0.48
to 0.50 under arm B, so the brightness division costs nothing once the three unshared channels
are out, and arm B agrees with A1 at ARI 0.59 and 0.73. The mass set recovers its HSC at 0.07
under every arm, and the shared markers do not resolve it there: a box on CD34, CD38 and CD90
around the HSC is 17 percent pure on the mass set against 64 percent on the flow set. The
comparison therefore does not reach the dye term. It says that the two HSC problems differ in
their parent population and in how CD38 stains, and both of those are bigger than any dye
effect.

## What ran

| Set | Platform | Cells | HSC | Parent population | Arms |
|---|---|---|---|---|---|
| Nilsson rare | Conventional flow, 13 colors, 10 used | 44,140 | 358, 0.81 percent | CD34-enriched bone marrow (the file's median CD34 is 3.3 on the arcsinh scale, so nearly every cell is CD34-positive) | A1, A2, A3, B (vendor ranks), C2, D |
| Levine 13 | CyTOF, 13 metals, 10 used | 167,044 | 261, 0.16 percent | Whole bone marrow, one donor | A1, A2, A3, C2, D, and B_meas with the CyTOF sensitivity per mass as the divisor |

Same FlowSOM settings as every run: 10 by 10 grid, 40 metaclusters, seed 1 on the full file,
10 seeds on a subsample for stability. The HSC is scored as a rare population, best-cluster
F1, on both sets, so the metric is the same. CD38 stays continuous in the binary arms.

The Levine 13 metals come from the HDCytoData channel names: CD45 on indium 115, CD45RA on
lanthanum 139, CD19 on 142, CD11b on 144, CD4 on 145, CD34 on 148, CD123 on 160, CD38 on
erbium 167, CD90 on 170, and CD3 on the cadmium channel (110 to 114). The sensitivity per mass
that arm B_meas divides by is in the section below.

## Results

### Nilsson rare, conventional flow, HSC 358 of 44,140 cells

| Arm | HSC F1 | Precision | Recall | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Phenotypes for 95 percent |
|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.50 | 0.45 | 0.56 | 0.59 |  |  |  |
| A2 z-score | 0.36 | 0.23 | 0.87 | 0.54 | 0.46 |  |  |
| A3 percentile | 0.37 | 0.23 | 0.90 | 0.61 | 0.54 |  |  |
| B brightness, 2^rank | 0.48 | 0.45 | 0.53 | 0.71 | 0.59 |  |  |
| B brightness, linear rank | 0.50 | 0.44 | 0.58 | 0.68 | 0.73 |  |  |
| C2 lineage binary + z-score | 0.37 | 0.32 | 0.44 | 0.80 | 0.45 | 247 of 512 | 65 |
| D all binary | 0.36 | 0.27 | 0.58 | 0.87 | 0.47 | 401 of 1024 | 100 |

### Levine 13, mass cytometry, HSC 261 of 167,044 cells

| Arm | HSC F1 | Precision | Recall | Stability ARI, 10 seeds | ARI against A1 | Occupied phenotypes | Phenotypes for 95 percent |
|---|---|---|---|---|---|---|---|
| A1 arcsinh | 0.07 | 0.04 | 0.21 | 0.70 |  |  |  |
| A2 z-score | 0.07 | 0.04 | 0.19 | 0.58 | 0.52 |  |  |
| A3 percentile | 0.11 | 0.07 | 0.26 | 0.62 | 0.46 |  |  |
| C2 lineage binary + z-score | 0.02 | 0.01 | 0.16 | 0.41 | 0.21 | 512 of 512 | 297 |
| D all binary | 0.02 | 0.01 | 0.15 | 0.42 | 0.17 | 1021 of 1024 | 545 |

### Levine 13 with the CyTOF sensitivity as the divisor, arm B_meas

The divisor per marker, geometric mean 1, from the sensitivity of each marker's mass channel
(Tricot et al. 2015, Cytometry A 87:357, Figure 2B, mean of three CyTOF 1 instruments; indium 115
extrapolated, the cadmium channel set neutral):

| Marker | Mass | Divisor |
|---|---|---|
| CD45 |  | 0.51 |
| CD45RA |  | 0.64 |
| CD19 |  | 0.72 |
| CD11b |  | 0.73 |
| CD4 |  | 0.79 |
| CD34 |  | 1.00 |
| CD123 |  | 2.29 |
| CD38 |  | 1.87 |
| CD90 |  | 1.84 |
| CD3 |  | 0.93 |

| Arm | HSC F1 | Precision | Recall | Stability ARI, 10 seeds | ARI against A1 |
|---|---|---|---|---|---|
| A1 arcsinh | 0.07 | 0.04 | 0.21 | 0.69 | |
| B measured divisor | 0.07 | 0.04 | 0.20 | 0.69 | 0.71 |

The divisor spans 4.5-fold, about the spread of the 2^rank divisors on the flow side, and it
moves the clustering (ARI 0.71 against A1) without moving the HSC. On the mass set as on the
flow set, the per-channel instrument term changes membership and does not change whether the
rare population is found. That is the same statement hypothesis 1 makes on the flow sets, now
with a measured divisor on a platform that has no dye.

## How the HSC differs between the two sets

Medians on the arcsinh scale of the HSC against all cells, with the separation as the
difference over the robust SD of all cells:

| Marker | Levine 13 HSC median | Levine 13 all cells | Separation | Nilsson HSC median | Nilsson all cells | Separation |
|---|---|---|---|---|---|---|
| CD34 | 2.44 | 0.27 | 3.95 | 2.98 | 3.33 | -0.49 |
| CD38 | 3.39 | 2.86 | 0.32 | 2.24 | 4.38 | -3.03 |
| CD90 | 1.25 | 0.51 | 0.93 | 3.06 | 2.98 | 0.12 |
| CD45RA | 0.32 | 1.26 | -0.52 | 0.99 | 2.88 | -2.33 |
| CD45 | 2.28 | 4.39 | -2.55 | 3.18 | 2.47 | 1.38 |

On the flow set the HSC is the CD38-low, CD45RA-low fraction of CD34-positive cells, and CD38
separates it by three robust SDs. On the mass set CD34 separates the HSC from whole marrow,
but CD38 does not separate it from the other CD34-positive cells: the HSC's CD38 median sits
above the all-cell median, and the 10th to 90th percentile range of the HSC on CD38 spans
2.4 to 4.4. Whether that is the erbium 167 stain, the gate the authors drew, or the biology of
this donor cannot be told from the file. The consequence for clustering is the same under
every arm: the mass HSC shares its position with 650 non-HSC cells on the three markers that
define it, and no representation separates cells that sit on top of each other.

## What this means for the paper

- The no-dye control is a weaker argument than it looks, because mass cytometry has its own
  per-channel term, the detector's sensitivity by mass, and because a matched pair on public
  data differs in donor, clone, enrichment and gate as well as in dye.
- On the one matched pair, the dye term is not what limits recovery on either platform. On
  flow the HSC is recovered at 0.50 with or without brightness division; on mass it is not
  recovered because the markers do not resolve it.
- The earlier arm B cost on Nilsson (0.50 to 0.22 on 13 markers) came through the three
  channels the comparison drops: CD10 on PE-Cy7, CD110 on APC and CD49f on Qdot 705. With
  those out, arm B and A1 agree. That is worth one sentence in the Discussion, because it
  says where arm B's damage comes from: the bright channels that the division pulls down.
- The stronger design is one sample stained once with the same clones and acquired on both
  platforms. No public set in the plan's tiers has that.

## Caveats

- One seed on the full file, as in every run.
- The flow divisors are vendor ranks with four Qdot values assumed, as in the pilot.
- The mass divisor is a published sensitivity curve for the instrument generation, not a
  measurement on the instrument that acquired this file.
- Nilsson rare is a CD34-enriched file and Levine 13 is whole marrow, so the two HSC
  populations have different parents and different frequencies, 0.81 and 0.16 percent.
