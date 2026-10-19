## Between-site variance, results/labels_hipc_site_HIPC_tcell_k30.npz, 7 centers, 3 lots, 63 files (0 files without manual statistics)

CD3 cut on the deposited scale for the T cell denominator: 1859

| Arm | Clusters | At 1 percent or more | Median between-site CV, all events | Lot 12828 | Lot 1349 | Lot 1369 | Median within-site CV, replicates | T cell fraction | Clusters at 1 percent of T cells | Median between-site CV, of T cells | Median within-site CV, of T cells | Median r to best cluster |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Central manual gates | 14 | 11 | 0.25 | 0.33 | 0.24 | 0.18 | 0.05 |  | 12 | 0.14 | 0.04 |  |
| A1 logicle as deposited | 30 | 20 | 0.47 | 0.50 | 0.49 | 0.25 | 0.08 | 0.46 | 13 | 0.30 | 0.06 | 0.77 |
| B brightness, 2^rank | 30 | 19 | 0.55 | 0.55 | 0.54 | 0.55 | 0.07 | 0.47 | 15 | 0.50 | 0.06 | 0.66 |
| B brightness, stain index pooled over centers | 30 | 19 | 0.44 | 0.49 | 0.44 | 0.37 | 0.08 | 0.47 | 11 | 0.30 | 0.05 | 0.69 |
| B brightness, stain index per center | 30 | 21 | 0.52 | 0.51 | 0.53 | 0.46 | 0.08 | 0.49 | 12 | 0.48 | 0.06 | 0.71 |

## Manual populations: between-site CV per lot (the reference)

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

## Best-matching cluster per manual population, per arm (r across files; between-site CV of the cluster against the population)

| Population | Manual CV | A1 r / CV | B_geom r / CV | B_pool r / CV | B_site r / CV |
|---|---|---|---|---|---|
| Lymphocytes | 0.18 | 0.74 / 0.18 | 0.65 / 0.55 | 0.67 / 0.40 | 0.67 / 0.58 |
| CD3 | 0.23 | 0.79 / 0.18 | 0.65 / 0.73 | 0.74 / 0.14 | 0.72 / 0.16 |
| CD4 | 0.24 | 0.90 / 0.18 | 0.77 / 0.16 | 0.87 / 0.14 | 0.80 / 0.16 |
| CD8 | 0.22 | 0.85 / 0.19 | 0.84 / 0.32 | 0.84 / 0.15 | 0.75 / 0.24 |
| CD4 Central Memory | 0.32 | 0.87 / 0.18 | 0.66 / 0.19 | 0.75 / 0.14 | 0.71 / 0.16 |
| CD4 Naive | 0.22 | 0.90 / 0.22 | 0.88 / 0.19 | 0.90 / 0.19 | 0.78 / 0.34 |
| CD8 Naive | 0.24 | 0.96 / 0.19 | 0.96 / 0.32 | 0.74 / 0.37 | 0.74 / 0.40 |
| CD4 Effector Memory | 0.26 | 0.75 / 0.23 | 0.41 / 0.73 | 0.58 / 0.14 | 0.86 / 0.26 |
| CD8 Effector | 0.47 | 0.89 / 0.50 | 0.86 / 0.47 | 0.61 / 0.15 | 0.53 / 0.32 |
| CD8 Effector Memory | 0.27 | 0.56 / 0.25 | 0.67 / 0.29 | 0.68 / 0.52 | 0.63 / 0.22 |
| CD8 Central Memory | 0.42 | 0.68 / 2.12 | 0.65 / 2.60 | 0.71 / 0.37 | 0.70 / 0.40 |
| CD4 Activated | 0.36 | 0.41 / 0.47 | 0.52 / 0.47 | 0.41 / 0.28 | 0.44 / 0.40 |
| CD8 Activated | 0.34 | 0.58 / 0.19 | 0.59 / 0.27 | 0.61 / 0.25 | 0.70 / 0.32 |
| CD4 Effector | 0.53 | 0.71 / 0.23 | 0.28 / 0.73 | 0.58 / 0.52 | 0.76 / 0.26 |
