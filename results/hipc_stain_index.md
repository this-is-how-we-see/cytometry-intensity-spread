# Stain index per center from the HIPC Lyoplate single-stain controls

**Run:** 2026-10-19 · **Code:** `code/hipc_stain_index.py` · **Data:** `data/HIPC/COMP/<Center>/`, the raw FACSDiva compensation controls of the seven centers with central gates · **Output:** `results/hipc_stain_index.csv` (one row per center and dye), `results/hipc_stain_index_skipped.csv` (the other controls on the same dye)

## What ran

The antibody controls are BD CompBeads, about half positive and half negative beads in one tube, 5,000 to 11,000 events. The Live/Dead FITC controls are cells. For each center and each of the eight T cell panel dyes the script takes the control that carries the panel's own reagent when the center has it, else another antibody on the same dye, and splits the tube on its primary detector with a two-component Gaussian mixture on arcsinh(x / 150). The stain index is the median positive minus the median negative over twice the robust SD of the negative, on the linear scale, uncompensated. On the cell tubes the positive is the brightest of three components, because the two-component cut sits inside the negative when the positive fraction is small. When a tube holds one population, the unstained tube of the center is the negative.

The index is a property of one conjugate on one instrument at one setting. On beads it carries the conjugate and the detector gain, and not the antigen. The control antibodies differ by center, so part of a between-center difference on a dye is the antibody and not the instrument. The FITC value is a Live/Dead stain on cells; Live is not a feature of the clustering and the value is here for completeness.

## Stain index by dye and center

The control marker is in parentheses. "none, pooled" marks a center without a control for that dye; arm B_site takes the mean over the other centers there.

| Dye | Baylor | CIMR | Miami | NHLBI | Stanford | UCLA | Yale | Mean over centers | Brightest over dimmest |
|---|---|---|---|---|---|---|---|---|---|
| FITC | none, pooled 144 | 91 (Live, cells) | 98 (Live, cells) | 196 (no marker named, cells) | 128 (Live, cells) | 138 (Live, cells) | 213 (Live, cells) | 144 | 2.3 |
| PE | 618 (CCR7) | 487 (CCR7) | 600 (CCR7) | none, pooled 579 | 281 (CCR7) | 960 (CCR7) | 530 (no marker named, cells) | 579 | 3.4 |
| PerCP-Cy5.5 | 30 (CD4) | 25 (CD4) | 19 (CD4) | 26 (no marker named) | 19 (CD4) | none, pooled 23 | 21 (CD4) | 23 | 1.6 |
| PE-Cy7 | 469 (CD45RA) | 504 (CD45RA) | 158 (CD45RA) | 590 (CD45RA) | 255 (CD194) | 619 (CD194) | 608 (CD45RA) | 457 | 3.9 |
| V450 / Pacific Blue | 112 (CD3) | 59 (CD3) | 98 (CD3) | 110 (no marker named) | 142 (CD3) | none, pooled 107 | 121 (CD3) | 107 | 2.4 |
| V500 / AmCyan | none, pooled 90 | 44 (HLA-DR) | 130 (HLA-DR) | 97 (no marker named) | 115 (HLA-DR) | 135 (HLA-DR) | 18 (HLA-DR) | 90 | 7.7 |
| APC | 253 (CD38) | 150 (CD38) | 139 (CD38) | 298 (CD38) | 283 (CD38) | none, pooled 228 | 244 (CD38) | 228 | 2.1 |
| APC-H7 / APC-Cy7 | 98 (CD20) | 67 (CD8) | 56 (CD8) | 83 (CD20) | 51 (CD20) | 138 (CD8) | 113 (CD8) | 87 | 2.7 |

## Every control used

| Center | Cytometer | Dye | Marker | File | Detector | Volts | Control | Events | Positive fraction | Negative | Median neg | Robust SD neg | Median pos | Stain index | Separation | Note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Baylor | FORTESSA (LSRII) H1773014 | PE-Cy7 | CD45RA | Compensation Controls_CD45RA PE-Cy7 Stained Control.fcs | PE-Cy7-A | 548 | beads | 10035 | 0.47 | in-tube | 3 | 20 | 19,069 | 469 | 6,485 |  |
| Baylor | FORTESSA (LSRII) H1773014 | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-H7 Stained Control.fcs | APC-H7-A | 735 | beads | 10015 | 0.46 | in-tube | 11 | 130 | 25,310 | 98 | 2,378 |  |
| Baylor | FORTESSA (LSRII) H1773014 | PerCP-Cy5.5 | CD4 | Compensation Controls_CD4 PerCP-Cy5-5 Stained Control.fcs | PerCP-Cy5-5-A | 518 | beads | 10034 | 0.49 | in-tube | 14 | 42 | 2,485 | 30 | 181 |  |
| Baylor | FORTESSA (LSRII) H1773014 | APC | CD38 | Compensation Controls_CD38 APC Stained Control.fcs | APC-A | 630 | beads | 10026 | 0.46 | in-tube | 66 | 47 | 24,018 | 253 | 362 |  |
| Baylor | FORTESSA (LSRII) H1773014 | PE | CCR7 | Compensation Controls_CD197 PE Stained Control.fcs | PE-A | 552 | beads | 10811 | 0.49 | in-tube | 41 | 49 | 61,112 | 618 | 1,484 |  |
| Baylor | FORTESSA (LSRII) H1773014 | FITC |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| Baylor | FORTESSA (LSRII) H1773014 | V450 / Pacific Blue | CD3 | Compensation Controls_CD3 V450 Stained Control.fcs | V450-A | 380 | beads | 11479 | 0.46 | in-tube | 81 | 34 | 7,758 | 112 | 95 |  |
| Baylor | FORTESSA (LSRII) H1773014 | V500 / AmCyan |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| CIMR | LSRFortessa H7B200001 | PE-Cy7 | CD45RA | Compensation Controls_CD45RA PE-Cy7 YG Stained Control.fcs | PE-Cy7 YG-A | 540 | beads | 5000 | 0.48 | in-tube | -1 | 25 | 24,760 | 504 |  |  |
| CIMR | LSRFortessa H7B200001 | APC-H7 / APC-Cy7 | CD8 | Compensation Controls_CD8 APC-Cy7 Stained Control.fcs | APC-Cy7-A | 620 | beads | 5000 | 0.50 | in-tube | -25 | 224 | 30,047 | 67 |  |  |
| CIMR | LSRFortessa H7B200001 | PerCP-Cy5.5 | CD4 | Compensation Controls_CD4 PerCP-Cy5-5 Stained Control.fcs | PerCP-Cy5-5-A | 555 | beads | 5000 | 0.48 | in-tube | 10 | 46 | 2,346 | 25 | 237 |  |
| CIMR | LSRFortessa H7B200001 | APC | CD38 | Compensation Controls_CD38 APC Stained Control.fcs | APC-A | 520 | beads | 5000 | 0.47 | in-tube | 275 | 85 | 25,688 | 150 | 93 |  |
| CIMR | LSRFortessa H7B200001 | PE | CCR7 | Compensation Controls_CCR7 PE YG Stained Control.fcs | PE YG-A | 600 | beads | 5000 | 0.50 | in-tube | 52 | 81 | 78,673 | 487 | 1,516 |  |
| CIMR | LSRFortessa H7B200001 | FITC | Live | Compensation Controls_LIVE_GREEN Alexa Fluor 488 Stained Control.fcs | Alexa Fluor 488-A | 480 | cells | 2428 | 0.20 | in-tube, 3 components | 88 | 107 | 19,628 | 91 | 223 |  |
| CIMR | LSRFortessa H7B200001 | V450 / Pacific Blue | CD3 | Compensation Controls_CD3 Pacific Blue Stained Control.fcs | Pacific Blue-A | 400 | beads | 5000 | 0.47 | in-tube | 76 | 63 | 7,481 | 59 | 98 |  |
| CIMR | LSRFortessa H7B200001 | V500 / AmCyan | HLA-DR | Compensation Controls_HLADR AmCyan Stained Control.fcs | AmCyan-A | 460 | beads | 5000 | 0.47 | in-tube | 155 | 170 | 15,127 | 44 | 97 |  |
| Miami | LSRII 1 | PE-Cy7 | CD45RA | Compensation Controls_CD45RA PE Cy7 YG Stained Control_F05.fcs | PE Cy7 YG-A | 540 | beads | 5000 | 0.53 | in-tube | 2 | 63 | 19,956 | 158 | 11,211 |  |
| Miami | LSRII 1 | APC-H7 / APC-Cy7 | CD8 | Compensation Controls_CD8 APC-Cy7 Stained Control_F12.fcs | APC-Cy7-A | 685 | beads | 5000 | 0.51 | in-tube | 6 | 226 | 25,461 | 56 | 4,359 |  |
| Miami | LSRII 1 | PerCP-Cy5.5 | CD4 | Compensation Controls_CD4 PerCP-Cy5-5 Stained Control_F04.fcs | PerCP-Cy5-5-A | 593 | beads | 5000 | 0.52 | in-tube | 15 | 65 | 2,434 | 19 | 156 |  |
| Miami | LSRII 1 | APC | CD38 | Compensation Controls_CD38 APC Stained Control_F10.fcs | APC-A | 599 | beads | 5000 | 0.52 | in-tube | 211 | 89 | 24,887 | 139 | 117 |  |
| Miami | LSRII 1 | PE | CCR7 | Compensation Controls_CD197 PE Stained Control_F03.fcs | PE-A | 462 | beads | 5000 | 0.53 | in-tube | 84 | 51 | 61,823 | 600 | 738 |  |
| Miami | LSRII 1 | FITC | Live | Compensation Controls_Live Green FITC Stained Control_F02.fcs | FITC-A | 517 | cells | 5000 | 0.23 | in-tube, 3 components | 160 | 113 | 22,323 | 98 | 138 |  |
| Miami | LSRII 1 | V450 / Pacific Blue | CD3 | Compensation Controls_CD3 Pacific Blue Stained Control_G04.fcs | Pacific Blue-A | 498 | beads | 5000 | 0.51 | in-tube | 77 | 39 | 7,722 | 98 | 100 |  |
| Miami | LSRII 1 | V500 / AmCyan | HLA-DR | Compensation Controls_HLA-DR AmCyan Stained Control_G05.fcs | AmCyan-A | 549 | beads | 5000 | 0.50 | in-tube | 186 | 58 | 15,252 | 130 | 81 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | PE-Cy7 | CD45RA | Compensation Controls_CD45RA PE-Cy7 G780 Stained Control.fcs | G780-A | 520 | beads | 5000 | 0.49 | in-tube | 1 | 17 | 20,227 | 590 | 14,872 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-H7 R780 Stained Control.fcs | R780-A | 660 | beads | 5000 | 0.50 | in-tube | 10 | 157 | 26,093 | 83 | 2,570 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | PerCP-Cy5.5 |  | Compensation Controls_B710 Stained Control.fcs | B710-A | 705 | beads | 5000 | 0.50 | in-tube | 37 | 47 | 2,511 | 26 | 66 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | APC | CD38 | Compensation Controls_CD38 APC R660 Stained Control.fcs | R660-A | 610 | beads | 5000 | 0.49 | in-tube | 121 | 45 | 26,750 | 298 | 220 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | PE |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | FITC |  | Compensation Controls_B515 Stained Control.fcs | B515-A | 500 | cells | 2027 | 0.36 | in-tube, 3 components | 139 | 54 | 21,519 | 196 | 154 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | V450 / Pacific Blue |  | Compensation Controls_V450 Stained Control.fcs | V450-A | 367 | beads | 5000 | 0.49 | in-tube | 113 | 33 | 7,351 | 110 | 64 |  |
| NHLBI | LSRFortessa 2012 (LSRII) H64717700044 | V500 / AmCyan |  | Compensation Controls_V545 Stained Control.fcs | V545-A | 527 | beads | 5000 | 0.48 | in-tube | 188 | 78 | 15,323 | 97 | 81 |  |
| Stanford | LSRII H47100031 | PE-Cy7 | CD194 | Compensation Controls_CD194 PE-Cy7 Stained Control.fcs | PE-Cy7-A | 590 | beads | 5000 | 0.48 | in-tube | 11 | 55 | 28,085 | 255 | 2,478 |  |
| Stanford | LSRII H47100031 | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-Cy7 Stained Control.fcs | APC-Cy7-A | 660 | beads | 5000 | 0.47 | in-tube | -6 | 216 | 22,001 | 51 |  |  |
| Stanford | LSRII H47100031 | PerCP-Cy5.5 | CD4 | Compensation Controls_CD4 PerCP-Cy5-5 Stained Control.fcs | PerCP-Cy5-5-A | 587 | beads | 5000 | 0.48 | in-tube | 4 | 55 | 2,131 | 19 | 516 |  |
| Stanford | LSRII H47100031 | APC | CD38 | Compensation Controls_CD38 APC Stained Control.fcs | APC-A | 572 | beads | 5000 | 0.46 | in-tube | 33 | 39 | 21,915 | 283 | 660 |  |
| Stanford | LSRII H47100031 | PE | CCR7 | Compensation Controls_CD197 PE Stained Control.fcs | PE-A | 430 | beads | 5000 | 0.49 | in-tube | 22 | 95 | 53,740 | 281 | 2,426 |  |
| Stanford | LSRII H47100031 | FITC | Live | Compensation Controls_Live,2f,Dead FITC Stained Control.fcs | FITC-A | 427 | cells | 3532 | 0.25 | in-tube, 3 components | 72 | 79 | 20,146 | 128 | 278 |  |
| Stanford | LSRII H47100031 | V450 / Pacific Blue | CD3 | Compensation Controls_CD3 Pacific Blue Stained Control.fcs | Pacific Blue-A | 339 | beads | 5000 | 0.43 | in-tube | 114 | 31 | 8,892 | 142 | 77 |  |
| Stanford | LSRII H47100031 | V500 / AmCyan | HLA-DR | Compensation Controls_HLA DR Am Cyan Stained Control.fcs | Am Cyan-A | 509 | beads | 5000 | 0.48 | in-tube | 175 | 71 | 16,536 | 115 | 94 |  |
| UCLA | LSRII H64717700004 | PE-Cy7 | CD194 | Compensation Controls_CD194 PE-Cy7 Stained Control.fcs | PE-Cy7-A | 555 | beads | 5000 | 0.51 | in-tube | 5 | 31 | 38,320 | 619 | 7,602 |  |
| UCLA | LSRII H64717700004 | APC-H7 / APC-Cy7 | CD8 | Compensation Controls_CD8 APC-H7 Stained Control.fcs | APC-H7-A | 540 | beads | 5000 | 0.48 | in-tube | 6 | 89 | 24,567 | 138 | 3,886 |  |
| UCLA | LSRII H64717700004 | PerCP-Cy5.5 |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| UCLA | LSRII H64717700004 | APC |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| UCLA | LSRII H64717700004 | PE | CCR7 | Compensation Controls_CD197 PE Stained Control.fcs | PE-A | 443 | beads | 5000 | 0.49 | in-tube | 104 | 31 | 59,534 | 960 | 569 |  |
| UCLA | LSRII H64717700004 | FITC | Live | Compensation Controls_LIVE DEAD FITC Stained Control.fcs | FITC-A | 366 | cells | 5000 | 0.21 | in-tube, 3 components | 58 | 77 | 21,422 | 138 | 367 |  |
| UCLA | LSRII H64717700004 | V450 / Pacific Blue |  |  | nan | nan | nan | 0 |  | nan |  |  |  |  |  | no control for this dye at this center |
| UCLA | LSRII H64717700004 | V500 / AmCyan | HLA-DR | Compensation Controls_HLA-DR BD Horizon V500 Stained Control.fcs | BD Horizon V500-A | 397 | beads | 5000 | 0.46 | in-tube | 115 | 52 | 14,252 | 135 | 123 |  |
| Yale | LSRII 1 | PE-Cy7 | CD45RA | Compensation Controls_CD45RA PE-Cy7 Stained Control.fcs | PE-Cy7-A | 505 | beads | 5000 | 0.48 | in-tube | 2 | 25 | 30,208 | 608 | 15,028 |  |
| Yale | LSRII 1 | APC-H7 / APC-Cy7 | CD8 | Compensation Controls_CD8 APC-Cy7 Stained Control.fcs | APC-Cy7-A | 550 | beads | 5000 | 0.51 | in-tube | 0 | 113 | 25,615 | 113 |  |  |
| Yale | LSRII 1 | PerCP-Cy5.5 | CD4 | Specimen_001_CD4 PerCP5-5.fcs | PerCP-Cy5-5-A | 560 | beads | 5000 | 0.50 | in-tube | 18 | 59 | 2,537 | 21 | 138 |  |
| Yale | LSRII 1 | APC | CD38 | Specimen_001_CD38 APC.fcs | APC-A | 475 | beads | 5000 | 0.51 | in-tube | 166 | 51 | 25,320 | 244 | 151 |  |
| Yale | LSRII 1 | PE |  | Compensation Controls_PE Stained Control.fcs | PE-A | 533 | cells | 5000 | 0.06 | in-tube, 3 components | 52 | 52 | 55,340 | 530 | 1,065 |  |
| Yale | LSRII 1 | FITC | Live | Specimen_001_Live Green.fcs | FITC-A | 475 | cells | 6060 | 0.03 | in-tube, 3 components | 36 | 44 | 18,739 | 213 | 513 |  |
| Yale | LSRII 1 | V450 / Pacific Blue | CD3 | Specimen_001_CD3 V450.fcs | V450-A | 322 | beads | 5000 | 0.53 | in-tube | 85 | 31 | 7,635 | 121 | 89 |  |
| Yale | LSRII 1 | V500 / AmCyan | HLA-DR | Specimen_001_HLA-DR V500.fcs | V500-A | 415 | beads | 5000 | 0.48 | in-tube | 863 | 408 | 15,239 | 18 | 17 |  |

## Controls on the same dye that the script did not use

| Center | Dye | Marker | File |
|---|---|---|---|
| Baylor | PE-Cy7 | CD11c | Compensation Controls_CD11c PE-Cy7 Stained Control.fcs |
| Baylor | PE-Cy7 | CD194 | Compensation Controls_CD194 PE-Cy7 Stained Control.fcs |
| Baylor | PE-Cy7 | CD196 | Compensation Controls_CD196 PE-Cy7 Stained Control.fcs |
| Baylor | PE-Cy7 | CD27 | Compensation Controls_CD27 PE-Cy7 Stained Control.fcs |
| Baylor | APC-H7 / APC-Cy7 | CD3+CD19+CD20 | Compensation Controls_CD3+19+20 APC-H7 Stained Control.fcs |
| Baylor | APC-H7 / APC-Cy7 | CD45RO | Compensation Controls_CD45RO APC-H7 Stained Control.fcs |
| Baylor | APC | CD127 | Compensation Controls_CD127 APC Stained Control.fcs |
| CIMR | PE-Cy7 | CCR4 | Compensation Controls_CCR4 PE-Cy7 YG Stained Control.fcs |
| CIMR | PE-Cy7 | CCR6 | Compensation Controls_CCR6 PE-Cy7 YG Stained Control.fcs |
| CIMR | PE-Cy7 | CD11c | Compensation Controls_CD11c PE-Cy7 YG Stained Control.fcs |
| CIMR | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-Cy7 Stained Control.fcs |
| CIMR | APC-H7 / APC-Cy7 | CD3+CD19+CD20 | Compensation Controls_CD3_CD19_CD20 APC-Cy7 Stained Control.fcs |
| CIMR | APC-H7 / APC-Cy7 | CD45RO | Compensation Controls_CD45RO APC-Cy7 Stained Control.fcs |
| CIMR | PerCP-Cy5.5 | CD123 | Compensation Controls_CD123 PerCP-Cy5-5 Stained Control.fcs |
| CIMR | PerCP-Cy5.5 | CD19 | Compensation Controls_CD19 PerCP-Cy5-5 Stained Control.fcs |
| CIMR | APC | CD127 | Compensation Controls_CD127 APC Stained Control.fcs |
| CIMR | APC | CD16 | Compensation Controls_CD16 APC Stained Control.fcs |
| CIMR | PE | CD24 | Compensation Controls_CD24 PE YG Stained Control.fcs |
| CIMR | PE | CD25 | Compensation Controls_CD25 PE YG Stained Control.fcs |
| CIMR | PE | CD56 | Compensation Controls_CD56 PE YG Stained Control.fcs |
| CIMR | PE | CXCR3 | Compensation Controls_CXCR3 PE YG Stained Control.fcs |
| CIMR | V450 / Pacific Blue | CD14 | Compensation Controls_CD14 Pacific Blue Stained Control.fcs |
| CIMR | V500 / AmCyan | IGD | Compensation Controls_IGD AmCyan Stained Control.fcs |
| Miami | PE-Cy7 | CD11c | Compensation Controls_CD11c PE Cy7 YG Stained Control_F08.fcs |
| Miami | PE-Cy7 | CD194 | Compensation Controls_CD194 PE Cy7 YG Stained Control_F06.fcs |
| Miami | PE-Cy7 | CD196 | Compensation Controls_CD196 PE Cy7 YG Stained Control_F09.fcs |
| Miami | PE-Cy7 | CD27 | Compensation Controls_CD27 PE Cy7 YG Stained Control_F07.fcs |
| Miami | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-Cy7 Stained Control_G02.fcs |
| Miami | APC-H7 / APC-Cy7 | CD3+CD19+CD20 | Compensation Controls_CD3+CD19+CD20+ APC-Cy7 Stained Control_G03.fcs |
| Miami | APC-H7 / APC-Cy7 | CD45RO | Compensation Controls_CD45RO APC-Cy7 Stained Control_G01.fcs |
| Miami | APC | CD127 | Compensation Controls_CD127 APC Stained Control_F11.fcs |
| NHLBI | PE-Cy7 | CCR4 | Compensation Controls_CCR4 PE-Cy7 G780 Stained Control.fcs |
| NHLBI | PE-Cy7 | CCR6 | Compensation Controls_CCR6 PE-Cy7 G780 Stained Control.fcs |
| NHLBI | PE-Cy7 | CD11c | Compensation Controls_CD11c PE-Cy7 G780 Stained Control.fcs |
| NHLBI | PE-Cy7 | CD27 | Compensation Controls_CD27 PE-Cy7 G780 Stained Control.fcs |
| NHLBI | APC-H7 / APC-Cy7 | CD3+CD19+CD20 | Compensation Controls_Lineage APC-H7 R780 Stained Control.fcs |
| NHLBI | APC | CD127 | Compensation Controls_CD127 AF647 R660 Stained Control.fcs |
| UCLA | PE-Cy7 | CD196 | Compensation Controls_CD196 PE-Cy7 Stained Control.fcs |
| UCLA | PE-Cy7 | CD27 | Compensation Controls_CD27 PE-Cy7 Stained Control.fcs |
| UCLA | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-H7 Stained Control.fcs |
| UCLA | APC-H7 / APC-Cy7 | CD45RO | Compensation Controls_CD45RO APC-H7 Stained Control.fcs |
| Yale | PE-Cy7 | CD45RA | Specimen_001_CD45RA PE-Cy7.fcs |
| Yale | PE-Cy7 | CD11c | Compensation Controls_CD11c PE-Cy7 Stained Control.fcs |
| Yale | PE-Cy7 | CD194 | Compensation Controls_CD194 PE-Cy7 Stained Control.fcs |
| Yale | PE-Cy7 | CD196 | Compensation Controls_CD196 PE-Cy7 Stained Control.fcs |
| Yale | PE-Cy7 | CD27 | Compensation Controls_CD27 PE-Cy7 Stained Control.fcs |
| Yale | APC-H7 / APC-Cy7 | CD20 | Compensation Controls_CD20 APC-Cy7 Stained Control.fcs |
| Yale | APC-H7 / APC-Cy7 | CD20 | Specimen_001_CD20 APC-H7.fcs |
| Yale | APC-H7 / APC-Cy7 | CD3+CD19+CD20 | Compensation Controls_CD3+CD19+CD20 APC-Cy7 Stained Control.fcs |
| Yale | APC-H7 / APC-Cy7 | CD45RO | Compensation Controls_CD45RO APC-Cy7 Stained Control.fcs |
| Yale | PerCP-Cy5.5 |  | Compensation Controls_PerCP-Cy5-5 Stained Control.fcs |
| Yale | APC | CD127 | Compensation Controls_CD127 Alexa-647 APC Stained Control.fcs |
| Yale | APC |  | Compensation Controls_APC Stained Control.fcs |
| Yale | FITC |  | Compensation Controls_FITC Stained Control.fcs |
| Yale | V450 / Pacific Blue |  | Compensation Controls_V450 Stained Control.fcs |
