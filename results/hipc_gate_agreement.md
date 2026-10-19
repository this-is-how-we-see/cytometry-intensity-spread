# HIPC Lyoplate T cell panel: agreement of the Python gating with the central manual popstats

**Run:** 2026-10-19 · **Code:** `code/flowjo9_gates.py . --write --report results/hipc_gate_agreement.csv` · **Table:** `results/hipc_gate_agreement.csv` (63 files, 15 populations each) · **Labels:** `data/HIPC/labels/<Center>/<file>.npz`

Relative difference of the Python count from the popstats count, in percent, as the range over the nine files of a center (min to max). A positive value means the Python gate holds more events. Live is the `not dead` population, the parent of CD3 in popstats.

| Population | Baylor | CIMR | Miami | NHLBI | Stanford | UCLA | Yale |
|---|---|---|---|---|---|---|---|
| Lymphocytes | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Live | +0.11 to +1.04 | 0 | +0.47 to +1.40 | +0.28 to +0.96 | 0 | +0.17 to +0.73 | -0.08 to +0.00 |
| CD3 | +0.11 to +0.89 | 0 | +0.45 to +1.18 | +0.22 to +0.43 | 0 | +0.10 to +0.62 | -0.08 to +0.00 |
| CD4 | +0.12 to +0.65 | 0 | +0.40 to +1.30 | +0.17 to +0.37 | 0 | +0.09 to +0.55 | -0.10 to +0.00 |
| CD8 | +0.06 to +0.67 | 0 | +0.31 to +0.71 | +0.18 to +0.36 | 0 | +0.09 to +0.65 | -0.03 to +0.01 |
| CD4 Naive | +0.23 to +1.04 | 0 | +1.25 to +3.24 | +0.24 to +0.57 | 0 | +0.08 to +0.89 | -0.08 to +0.00 |
| CD4 Central Memory | +0.05 to +0.24 | 0 | +0.09 to +0.30 | +0.10 to +0.18 | +0.00 to +0.01 | +0.06 to +0.47 | -0.18 to +0.00 |
| CD4 Effector Memory | +0.01 to +0.89 | 0 | +0.00 to +0.21 | +0.02 to +0.21 | +0.00 to +0.05 | +0.11 to +0.33 | 0 |
| CD4 Effector | +0.00 to +22.22 | 0 | +0.00 to +10.00 | +0.00 to +4.76 | 0 | +0.19 to +2.48 | 0 |
| CD8 Naive | +0.11 to +0.99 | 0 | +0.59 to +1.02 | +0.29 to +0.80 | 0 | +0.19 to +1.12 | -0.04 to +0.02 |
| CD8 Central Memory | +0.00 to +0.36 | +0.00 to +0.05 | +0.00 to +0.35 | +0.00 to +0.37 | 0 | +0.09 to +0.55 | -0.08 to +0.00 |
| CD8 Effector Memory | +0.00 to +0.20 | +0.00 to +0.02 | +0.02 to +0.12 | +0.03 to +0.29 | 0 | +0.00 to +0.53 | 0 |
| CD8 Effector | +0.00 to +0.73 | 0 | +0.00 to +0.43 | +0.04 to +0.53 | 0 | +0.03 to +0.71 | -0.06 to +0.00 |
| CD4 Activated | +1.98 to +9.68 | 0 | +6.63 to +19.09 | +3.04 to +9.76 | 0 | +0.90 to +9.72 | -0.35 to +0.00 |
| CD8 Activated | +0.75 to +14.29 | 0 | +3.88 to +8.87 | +1.69 to +7.36 | 0 | +0.79 to +10.49 | 0 |

Singlets gate: applied at Stanford (FSC-H from the raw files), none in the workspace at CIMR and Yale, skipped at Baylor, Miami, NHLBI and UCLA because the deposited files carry no FSC-H or FSC-W.

## Per file, CD3

| Center | File | Python | popstats | FlowJo | Difference |
|---|---|---|---|---|---|
| Baylor | 12828_1_T CELL.fcs | 11635 | 11532 | 11522 | +0.89 % |
| Baylor | 12828_2_T CELL.fcs | 106512 | 106252 | 106307 | +0.24 % |
| Baylor | 12828_3_T CELL.fcs | 92043 | 91744 | 91757 | +0.33 % |
| Baylor | 1349_1_T CELL.fcs | 107259 | 106912 | 106920 | +0.32 % |
| Baylor | 1349_2_T CELL.fcs | 121316 | 121013 | 120982 | +0.25 % |
| Baylor | 1349_3_T CELL.fcs | 110319 | 110090 | 110034 | +0.21 % |
| Baylor | 1369_1_T CELL.fcs | 87747 | 87562 | 87559 | +0.21 % |
| Baylor | 1369_2_T CELL.fcs | 95340 | 95182 | 95225 | +0.17 % |
| Baylor | 1369_3_T CELL.fcs | 91676 | 91578 | 91660 | +0.11 % |
| CIMR | T_CELL_1369_P1.fcs | 50180 | 50180 | 50182 | +0.00 % |
| CIMR | T_CELL_1369_001_P1.fcs | 50630 | 50630 | 50620 | +0.00 % |
| CIMR | T_CELL_1369_002_P1.fcs | 50291 | 50291 | 50286 | +0.00 % |
| CIMR | T_CELL_1349_P1.fcs | 54995 | 54995 | 54999 | +0.00 % |
| CIMR | T_CELL_1349_001_P1.fcs | 54728 | 54728 | 54728 | +0.00 % |
| CIMR | T_CELL_1349_002_P1.fcs | 54730 | 54730 | 54733 | +0.00 % |
| CIMR | T_CELL_12828_P1.fcs | 68071 | 68071 | 68062 | +0.00 % |
| CIMR | T_CELL_12828_001_P1.fcs | 57674 | 57674 | 57721 | +0.00 % |
| CIMR | T_CELL_12828_002_P1.fcs | 41952 | 41952 | 41955 | +0.00 % |
| Miami | lot 12828_A1_A01.fcs | 18062 | 17851 | 17851 | +1.18 % |
| Miami | lot 12828_A2_A02.fcs | 37875 | 37573 | 37572 | +0.80 % |
| Miami | lot 12828_A3_A03.fcs | 37603 | 37327 | 37323 | +0.74 % |
| Miami | lot 1349_A4_A04.fcs | 48330 | 48077 | 48078 | +0.53 % |
| Miami | lot 1349_A5_A05.fcs | 47688 | 47464 | 47466 | +0.47 % |
| Miami | lot 1349_A6_A06.fcs | 47079 | 46869 | 46870 | +0.45 % |
| Miami | lot 1369_A7_A07.fcs | 49404 | 49092 | 49090 | +0.64 % |
| Miami | lot 1369_A8_A08.fcs | 49598 | 49243 | 49241 | +0.72 % |
| Miami | lot 1369_A9_A09.fcs | 48074 | 47734 | 47737 | +0.71 % |
| NHLBI | 12828_1_Tcell_A01.fcs | 29265 | 29200 | 29194 | +0.22 % |
| NHLBI | 12828_2_Tcell_A02.fcs | 51222 | 51002 | 51005 | +0.43 % |
| NHLBI | 12828_3_Tcell_A03.fcs | 51468 | 51268 | 51263 | +0.39 % |
| NHLBI | 1349_1_Tcell_A04.fcs | 44310 | 44214 | 44208 | +0.22 % |
| NHLBI | 1349_2_Tcell_A05.fcs | 42807 | 42709 | 42693 | +0.23 % |
| NHLBI | 1349_3_Tcell_A06.fcs | 42514 | 42399 | 42385 | +0.27 % |
| NHLBI | 1369_1_Tcell_A07.fcs | 51583 | 51419 | 51410 | +0.32 % |
| NHLBI | 1369_2_Tcell_A08.fcs | 51373 | 51214 | 51200 | +0.31 % |
| NHLBI | 1369_3_Tcell_A09.fcs | 51482 | 51338 | 51331 | +0.28 % |
| Stanford | 1228-1_A1_A01.fcs | 32238 | 32238 | 32318 | +0.00 % |
| Stanford | 1228-2_A2_A02.fcs | 38180 | 38180 | 38301 | +0.00 % |
| Stanford | 1228-3_A3_A03.fcs | 36363 | 36363 | 36436 | +0.00 % |
| Stanford | 1349-1_A4_A04.fcs | 34836 | 34836 | 34903 | +0.00 % |
| Stanford | 1349-2_A5_A05.fcs | 34239 | 34239 | 34323 | +0.00 % |
| Stanford | 1349-3_A6_A06.fcs | 40991 | 40991 | 41084 | +0.00 % |
| Stanford | 1369-1_A7_A07.fcs | 38361 | 38361 | 38410 | +0.00 % |
| Stanford | 1369-2_A8_A08.fcs | 38332 | 38332 | 38378 | +0.00 % |
| Stanford | 1369-3_A9_A09.fcs | 38546 | 38546 | 38602 | +0.00 % |
| UCLA | TCELL 22013_12828_001.fcs | 22174 | 22151 | 22140 | +0.10 % |
| UCLA | TCELL 22013_12828_002.fcs | 45631 | 45541 | 45501 | +0.20 % |
| UCLA | TCELL 22013_12828_003.fcs | 46509 | 46406 | 46396 | +0.22 % |
| UCLA | TCELL 22013_1349_001.fcs | 48291 | 48065 | 48096 | +0.47 % |
| UCLA | TCELL 22013_1349_002.fcs | 48163 | 47949 | 47969 | +0.45 % |
| UCLA | TCELL 22013_1349_003.fcs | 47187 | 46988 | 46993 | +0.42 % |
| UCLA | TCELL 22013_1369_001.fcs | 49720 | 49482 | 49499 | +0.48 % |
| UCLA | TCELL 22013_1369_002.fcs | 49191 | 48948 | 48978 | +0.50 % |
| UCLA | TCELL 22013_1369_003.fcs | 48686 | 48388 | 48400 | +0.62 % |
| Yale | 12828_1_A1_A01.fcs | 22378 | 22390 | 22361 | -0.05 % |
| Yale | 12828_2_A2_A02.fcs | 22716 | 22727 | 22701 | -0.05 % |
| Yale | 12828_3_A3_A03.fcs | 23422 | 23437 | 23418 | -0.06 % |
| Yale | 1349_1_A4_A04.fcs | 19472 | 19488 | 19448 | -0.08 % |
| Yale | 1349_2_A5_A05.fcs | 20593 | 20603 | 20583 | -0.05 % |
| Yale | 1349_3_A6_A06.fcs | 16717 | 16717 | 16695 | +0.00 % |
| Yale | 1369_1_A7_A07.fcs | 20732 | 20743 | 20724 | -0.05 % |
| Yale | 1369_2_A8_A08.fcs | 17751 | 17751 | 17740 | +0.00 % |
| Yale | 1369_3_A9_A09.fcs | 17317 | 17317 | 17300 | +0.00 % |
