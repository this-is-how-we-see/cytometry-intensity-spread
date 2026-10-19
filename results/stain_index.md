# Measured stain index against the vendor brightness rank

**Run:** 2026-10-18 · **Code:** `code/stain_index.py` · **Data:** the DURAClone IM Basic run of 10-DEC-2024 on the LSRFortessa, 12 compensation controls and one unstained control, read in place from `<DURAClone folder>/` (permission to publish that set stays open, and the script copied nothing) · **Inventory:** `data/README.md`

## What ran

Each compensation control is a tube of capture beads, 5,000 events, one antibody conjugate per tube. The tubes hold 95 to 99 percent positive beads, so a tube has no negative population of its own. The negative for every channel is the unstained bead tube. The positive beads are the events above the unstained median plus 3 robust SD on arcsinh(x / 150). The stain index is the median positive minus the median negative, over twice the robust SD of the negative. The script computes it on the linear scale, on the uncompensated primary detector.

A single stain carries one dye, so its own detector needs no compensation. The column "SI compensated" applies the stored matrix of the stained sample files as a check. "Separation" is the median positive minus the median negative, over the median negative. The unstained median sits between 0 and 40 on nine channels, so that ratio is unstable there. The table keeps it only because the plan asked for it.

The two-component Gaussian mixture that the plan named splits the positive mode in two on these tubes, because the tubes have no negative beads. The script reports the in-tube index only when the smaller component holds at least 10 percent of events. That happens on three tubes (BUV661, Krome Orange, APC). Those three values are splits inside the positive mode and are no stain index. They are in `stain_index.csv` and out of the table.

## Table

Sorted by stain index. The spill pattern of a tube is the median of the tube minus the median of the unstained tube, in every detector. The script divides that by the same difference in the tube's own detector. "Spill check" is the largest difference between that pattern and the stored SPILL row of the stained sample files, with the detector where it occurs. A difference above 0.5 means the stored matrix did not come from this tube as exported.

| Fluorochrome | Marker | Detector | Volts | Positive beads | Median neg | Robust SD neg | Median pos | Stain index | SI compensated | Separation | Spill check | Vendor rank | Rank source |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BUV737 | CD69 | BUV737-A | 687 | 4,983 | 186 | 188 | 1,602,383 | 4,265 | 4,255 | 8,635 | 0.09 in Alexa 700-A | 4 | fluor-ai |
| BUV661 | PD-L1 | BUV661/Qdot655-A | 495 | 4,913 | 17 | 24 | 112,287 | 2,324 | 1,606 | 6,579 | **3.37 in BUV737-A** | 4 | fluor-ai |
| BV786 | HLA-DR | BV786-A | 416 | 4,981 | 1 | 11 | 41,788 | 1,834 | 1,815 | 74,083 | **1.60 in BUV737-A** | 4 | fluor-ai |
| FITC | CD16 | BB515/FITC-A | 360 | 4,961 | 2 | 13 | 33,271 | 1,299 | 1,306 | 19,669 | 0.02 in BUV737-A | 3 | fluor-ai |
| APC-Alexa Fluor 750 | CD3 | APC-Cy7-A | 386 | 4,936 | 0 | 10 | 20,977 | 1,075 | 816 | 417,406 | **0.84 in BUV737-A** | 3 | assumed, not in fluor-ai; rated with the other APC tandems |
| PE | CD56 | PE-A | 410 | 4,946 | 108 | 34 | 65,685 | 964 | 880 | 608 | 0.06 in BV570-A | 5 | fluor-ai |
| Alexa Fluor 700 | CD8 | Alexa 700-A | 394 | 4,975 | 19 | 17 | 29,678 | 871 | 370 | 1,523 | **1.16 in BUV737-A** | 3 | fluor-ai |
| Viability dye, violet | LD | BV480-A | 288 | 4,971 | 36 | 27 | 43,427 | 790 | 797 | 1,205 | 0.11 in BV570-A | 3 | assumed; the FCS keywords name no product, fluor-ai rates every violet viability dye 3 |
| PE-Texas Red | CD19 | PE-CF594/ECD-A | 380 | 4,892 | 30 | 28 | 30,951 | 543 | -26 | 1,028 | **2.07 in PE-A** | 4 | fluor-ai |
| PE-Cy7 | CD14 | BB790-A | 522 | 4,933 | 4 | 33 | 32,215 | 481 | 419 | 8,312 | 0.07 in BUV737-A | 5 | fluor-ai |
| Krome Orange | CD45 | BV570-A | 418 | 4,617 | 20 | 25 | 14,033 | 276 | 278 | 691 | 0.28 in BUV737-A | 3 | assumed, not in fluor-ai |
| APC | CD4 | APC-A | 389 | 4,820 | 568 | 62 | 10,796 | 83 | 68 | 18 | **2.12 in BUV737-A** | 5 | fluor-ai |

## Comparison with the vendor rank

| Set | n | Spearman rho | p | Brightest / dimmest, measured | Brightest / dimmest, 2^rank | Brightest / dimmest, linear rank |
|---|---|---|---|---|---|---|
| all 12 tubes | 12 | -0.14 | 0.67 | 52 (BUV737 / APC) | 4 (rank 5 / 3) | 1.67 |
| 11 antibody conjugates | 11 | -0.19 | 0.57 | 52 (BUV737 / APC) | 4 (rank 5 / 3) | 1.67 |
| tubes consistent with the matrix | 6 | 0.09 | 0.86 | 15 (BUV737 / Krome Orange) | 4 (rank 5 / 3) | 1.67 |
| antibodies consistent with the matrix | 5 | 0.00 | 1.00 | 15 (BUV737 / Krome Orange) | 4 (rank 5 / 3) | 1.67 |

Members of the two consistent sets: BUV737, FITC, PE, Viability dye, violet, PE-Cy7, Krome Orange.

Divisors with geometric mean 1, as arm B applies them before arcsinh, for the 11 antibody conjugates. Arm B divides the channel by the divisor. A value above 1 pulls a bright channel down, and a value below 1 lifts a dim channel up.

| Fluorochrome | Marker | Measured SI | Divisor, measured SI | Divisor, 2^rank | Divisor, linear rank | Spill check |
|---|---|---|---|---|---|---|
| BUV737 | CD69 | 4,265 | 5.14 | 1.07 | 1.04 | passes |
| BUV661 | PD-L1 | 2,324 | 2.80 | 1.07 | 1.04 | fails |
| BV786 | HLA-DR | 1,834 | 2.21 | 1.07 | 1.04 | fails |
| FITC | CD16 | 1,299 | 1.57 | 0.53 | 0.78 | passes |
| APC-Alexa Fluor 750 | CD3 | 1,075 | 1.30 | 0.53 | 0.78 | fails |
| PE | CD56 | 964 | 1.16 | 2.13 | 1.31 | passes |
| Alexa Fluor 700 | CD8 | 871 | 1.05 | 0.53 | 0.78 | fails |
| PE-Texas Red | CD19 | 543 | 0.66 | 1.07 | 1.04 | fails |
| PE-Cy7 | CD14 | 481 | 0.58 | 2.13 | 1.31 | passes |
| Krome Orange | CD45 | 276 | 0.33 | 0.53 | 0.78 | passes |
| APC | CD4 | 83 | 0.10 | 2.13 | 1.31 | fails |

## What the numbers say

The measured stain index does not agree with the vendor rank on this run. Spearman rho is
-0.19 on the 11 antibody conjugates and 0.00 on the five antibodies whose tubes pass the
spill check. The vendor table puts PE, APC and PE-Cy7 at rank 5 above everything else. The
beads put BUV737 first at 4,265, then BUV661, BV786 and FITC, with PE sixth at 964, PE-Cy7
tenth at 481 and APC last at 83.

The measured spread is also far wider than either mapping. Brightest over dimmest is 52 on the 11 antibodies and 15 on the five that pass the spill
check. The 2^rank mapping spans 4 and the linear mapping spans 1.67.

The 2^rank mapping in pilot.py is therefore not close to the measured spread, and the
direction of the divisors is also different. On this panel the rank mapping divides PE, APC
and PE-Cy7 by 2.13 and lifts FITC, Krome Orange and the APC tandem by 0.53. The measured
divisors divide BUV737 by 5.14 and the three Brilliant dyes by 2.2 to 5.1, leave PE at 1.16,
and lift PE-Cy7 to 0.58 and APC to 0.10. A paper that runs arm B with vendor ranks and a paper
that runs it with these single stains make opposite moves on three of the eight DURAClone
channels.

Part of the gap is the instrument and not the dye. The stain index of a bead is the dye
times the conjugate times the detector setting. The BUV737 detector ran at 687 volts and the
PE detector at 410, and on the LSRFortessa the UV detectors need the high voltage, so the
BUV737 index carries the gain. That is the correct quantity for arm B, because the stained
samples carry the same gain, and arm B is meant to remove the scale the acquisition put on
each channel. The vendor rank cannot see the voltage at all.

Part of the gap is this run. Six tubes show a spill pattern that the stored SPILL matrix does
not describe (spill check above 0.5): the ECD tube shows 2.3 times more signal in the PE
detector than in its own, which is the pattern of a PE conjugate or of a broken PE tandem,
and the APC, Alexa Fluor 700, APC-Alexa Fluor 750, BV786 and BUV661 tubes show 0.9 to 3.6
times their own signal in the BUV737 detector, where the stored matrix says 0.05 to 0.27.
Those five all emit in the far red, and a UV detector that reads a far-red emitter at that
level points to a laser delay or window setting, which the files do not record. The stored
matrix came from some other acquisition of the controls, and this export cannot say which one
is right. The ECD, APC and BUV737 values
in the table are therefore the least trustworthy, and the five-antibody set is the one to
quote. On that set the rank correlation is still zero and the spread is still 15.

## Divisor rule for arm B when single stains exist

Use the measured stain index as the divisor on the linear scale before arcsinh, with the
geometric mean of the divisors set to 1 so that the cofactor keeps its meaning:

    divisor_j = SI_j / geomean(SI)        x'_j = x_j / divisor_j        arcsinh(x'_j / cofactor)

Compute SI_j on the uncompensated primary detector of each single-stain control, with the
unstained tube as the negative and the same positive threshold as this script. Run the spill
check first, and drop a control whose pattern fails it; for that channel fall back to the
vendor rank and say so in the table. Report the rank correlation between the measured index
and the vendor rank as a quality figure, because a value near zero, as here, says the rank
mapping was not a usable stand-in for that run. Fall back to the 2^rank mapping only when no
controls exist, and say in the Methods that the rank mapping spans a factor of 4 where
measured values span a factor of 15 to 50.

## Limits

A stain index is a property of one conjugate on one instrument at one setting, and it is not
a property of the dye alone. Each control here is a different antibody, so the index mixes
the dye brightness with the fluorochrome-to-protein ratio of that lot and with the amount of
antibody the capture bead binds. On cells the index also mixes in the antigen density, and the
stained samples of this run are not single stains, so this run cannot separate the two. A
divisor from bead controls corrects for the conjugate and the detector setting, and a divisor
from cell controls corrects for the antigen density as well. The paper has to say which one it
uses and why, and on a cell control the Krome Orange CD45 channel, the dimmest conjugate
here, resolves well because CD45 is among the densest antigens on a leukocyte.

Two smaller limits. The viability dye is not an antibody, and its bead control (BV480
detector, index 790) is a stain of the bead surface, so it is in the table and out of the
antibody comparisons. And the unstained tube reads 568 on the APC detector with a robust SD of
62, against 0 to 40 on nine other detectors, so the APC index of 83 rests on a negative that
is itself unusual.
