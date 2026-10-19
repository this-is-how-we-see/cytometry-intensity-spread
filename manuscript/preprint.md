# Intensity spread is information: what FlowSOM loses when cytometry data are scaled, binarized or brightness-normalized

**Howie Seay, SCYM(ASCP)** · HOW WE SEE LLC · howie@thisishowwesee.com · ORCID [0000-0001-7002-0187](https://orcid.org/0000-0001-7002-0187)

**Preprint, draft of 2026-10-19.** Code and derived results: https://github.com/this-is-how-we-see/cytometry-intensity-spread. Data: public benchmark sets as cited, and the HIPC Lyoplate files on ImmuneSpace.

## Abstract

Clustering methods for high-dimensional cytometry take transformed intensity as their input, and most pipelines offer optional steps that change that input before clustering: per-channel scaling, binarization, or a correction for the brightness of the conjugate. This study asks what each of those steps does to the populations the clustering finds. Six public benchmark datasets with manual gates (three mass cytometry sets from Levine et al. and Samusik et al., two conventional flow sets with a rare population, and the seven-center HIPC Lyoplate T cell panel) were clustered with FlowSOM from seven representations of the same cells: arcsinh intensity, z-scored and percentile-scaled intensity, intensity divided by conjugate brightness, lineage markers as bits, every marker as a bit, and intensity reweighted by each marker's resolution. Recovery was scored by mean F1 with Hungarian matching against the manual gates, as in the Weber and Robinson benchmark. Plain arcsinh intensity recovered the gated populations best or near best on every set. Z-scoring cost 0.05 to 0.16 of mean F1 at the free cluster number and up to 0.25 of adjusted Rand index; percentile clipping erased a population of 0.03 percent. Binary input lost on every mass cytometry set under unsupervised thresholds and under thresholds taken from the gates themselves, and the loss fell on populations that a marker defines at an intermediate level or that differ from a neighbor by one marker. Division by conjugate brightness changed cluster membership and did not improve recovery. Across seven centers that stained one set of samples with the same lyophilized reagents, the between-site spread of cluster frequency under plain intensity matched the central manual gates, and no representation lowered it. The measured spread of each channel carries information that every one of these steps removes, and the cost is paid by the intermediate and rare populations. As a secondary finding, the deposited HIPC T cell files carry FlowJo's per-channel transform rather than the one their README states, and the reconstruction is given.

## 1. Introduction

Automated analysis of high-dimensional cytometry runs on transformed intensity. FlowSOM [2], PhenoGraph [3], X-shift [4] and the pipelines built on them [13, 14] take a matrix of cells by markers, apply an inverse hyperbolic sine with a fixed cofactor, and cluster on Euclidean distance in that space. Weber and Robinson [1] benchmarked eighteen clustering methods on six public datasets with manual gates and fixed this input for every method. The benchmark settled which method to use. It did not ask what the input should be.

The input is a choice, and the choice has a physical reading. On a conventional flow cytometer the intensity of a marker is the product of the antigen density on the cell and the brightness of the fluorochrome that the panel designer assigned to that antigen. The second term is a design decision, and it differs by more than tenfold across the dyes of a large panel. A clustering on Euclidean distance lets a bright conjugate pull harder than a dim one, so the panel design is a hidden weight on every published clustering. Manual gating never had this problem, because a gate is a threshold and a threshold is invariant to brightness once it is set.

Several common practices change the input to remove or reduce that weight, each for its own reason. Many pipelines offer per-channel scaling, z-score or percentile, as an optional step before clustering. Single-cell protein analysis in CITE-seq, where no fluorochrome exists, calls each protein positive or negative with a mixture model and clusters the calls. flowType [10] enumerated every positive and negative combination of markers as a phenotype and found correlates of outcome among them. Clinical immunophenotype scores are binary marker calls by design. Each of these treats the intensity of a marker as a nuisance and keeps only its sign or its rank.

This study tests those practices head to head on the same cells. The hypothesis was that a representation which removes the brightness term, by division or by binarization, recovers the manually gated populations at least as well as intensity does, with fewer clusters and better run-to-run stability. The study also asked where such a representation fails, and named the cases in advance: populations that a marker defines at an intermediate level, such as CD4 on monocytes. The scope is the input to an existing algorithm. No new clustering method is proposed, and no new software beyond the scripts that reproduce the results.

The hypothesis did not hold. Plain arcsinh intensity recovered the gated populations best on every set, and every representation that removed or reweighted the measured spread of a channel lost, in a pattern that repeated across datasets and platforms. The paper reports that pattern, its mechanism, and two results that stand on their own: the cost of the optional scaling step, and the between-site behavior of each representation on a multi-center dataset with central gates.

## 2. Methods

### 2.1 Datasets

Table 1 lists the datasets. Four are mass cytometry sets with multi-population manual gates from the Weber and Robinson benchmark [1]: Levine_32dim and Levine_13dim, human bone marrow from Levine et al. [3] on the panels of Bendall et al. [5], and Samusik_01 and Samusik_all, mouse bone marrow from Samusik et al. [4], one mouse and all ten. Two are conventional flow sets from the same benchmark with one gated rare population each: Nilsson_rare, hematopoietic stem cells at 0.8 percent of a CD34-enriched human marrow sample [6], and Mosmann_rare, activated memory CD4 T cells at 0.03 percent of influenza-stimulated human PBMC [7]. The Levine files came from the benchmark's GitHub mirror and the rest from the HDCytoData copies on Bioconductor ExperimentHub, which keep the original channel names. The HDCytoData copy of Nilsson_rare carries its fluorochromes in the channel names; the Mosmann_rare panel came from the MIFlowCyt checklist published as the first supplementary file of Mosmann et al. [7], experiment JMW034, which is the experiment FlowRepository FR-FCM-ZZ8J holds.

The cross-site dataset is the T cell panel of the HIPC Lyoplate standardization study [8], from ImmuneSpace project HIPC/Lyoplate: three lots of SeraCare PBMC stained in triplicate at seven centers with the same lyophilized BD reagents, acquired on BD LSR II and LSRFortessa instruments, and gated centrally by one analyst in FlowJo 9.6.3. Eight markers: CD3, CD4, CD8, CD45RA, CCR7, CD38, HLA-DR and a viability dye. The deposited files, 63 for the seven centers with central gates, are not raw, and section 2.6 says what they are.

One dataset, provided in kind by the University of Florida, contributed a measurement and no clustering result. It holds twelve single-stain bead compensation controls from a BD LSRFortessa, acquired in December 2024 with the Beckman Coulter DURAClone IM Basic tube and drop-in reagents, used to measure a stain index per conjugate.

### 2.2 The representations

Every representation starts from the deposited intensity, compensated where the depositor compensated it. Table 2 defines the arms.

- **A1, arcsinh intensity.** The inverse hyperbolic sine of intensity over a cofactor, 5 for mass cytometry and 150 for fluorescence, as in Weber and Robinson [1]. This is the control arm and the field's default.
- **A2, z-score.** A1, then each channel centered on its mean and divided by its standard deviation.
- **A3, percentile.** A1, then each channel scaled from its first to its 99th percentile onto the interval from 0 to 1, with values outside clipped.
- **B, brightness-normalized.** Intensity divided by the relative brightness of the conjugate on the linear scale, then transformed as A1. Division before the transform is the only order with a physical meaning. Two divisor rules mapped a vendor brightness rank from 1 to 5 onto a multiplicative factor, 2 to the power of the rank minus 3 and the rank over 3, each normalized to a geometric mean of 1 so the cofactor keeps its meaning. The ranks came from vendor brightness tables; four Qdot values on Nilsson_rare and five on Mosmann_rare were assumed. On the HIPC panel a third divisor was measured on each center's own single-stain controls.
- **C and C2, lineage markers as bits.** Each lineage marker cut to 0 or 1 at the equal-posterior point of a two-component Gaussian mixture fitted to its arcsinh values on a subsample of 50,000 cells; functional markers kept continuous, percentile-scaled in C and z-scored in C2. The lineage and functional split is stated per dataset in the repository.
- **D, every marker as a bit.** All markers cut as in C, in the manner of flowType [10].
- **Gate thresholds.** On Levine_32dim the binary arms also ran with a supervised threshold per marker, the cut that maximizes balanced accuracy between the populations the label set defines as positive and as negative, as the upper bound a binary representation can reach on that set. This arm uses the labels and is reported as an upper bound only.
- **W and P, reweighted intensity.** Two unsupervised weightings tested whether a graded weight does what a bit cannot. W multiplies each A1 channel by the separation of its two mixture modes in pooled-SD units; W2 multiplies by the sharpness of the mixture posterior, the share of cells scored below 0.1 or above 0.9; both normalized to a geometric mean of 1 and clipped to the interval from 1/3 to 3. P replaces each bit of C2 with the posterior probability of the high component.

### 2.3 Clustering

FlowSOM [2], in the Python port maintained by the Saeys lab, with a 10 by 10 grid and consensus metaclustering to a fixed number of metaclusters. Every dataset ran at the number of populations in its manual gates, the gated number, and at a free number, 40 for the benchmark sets and 30 for the HIPC panel, which is how FlowSOM is used in practice. The claims in this paper are made at the free number; the gated number is reported because merge-order effects at a forced count reverse some comparisons, and the Results say where. One seed on the full data gave the scored clustering. Ten seeds on a subsample of 50,000 cells gave the run-to-run stability as the mean pairwise adjusted Rand index. Leiden on a shared-nearest-neighbor graph, with Hamming distance for the binary arms, served as a second clustering method on Levine_32dim.

### 2.4 Scoring

Recovery of the gated populations is the mean F1 under Hungarian matching [20], exactly as in Weber and Robinson [1]: each manual population is matched to one cluster so that the summed F1 is as large as possible, and the F1 values are averaged with every population weighing the same. The adjusted Rand index [19] and normalized mutual information are secondary. For the two rare-population sets the score is the F1 of the best single cluster for the rare population. Per-population F1 is reported for the populations the hypothesis named. For the binary arms the number of occupied Boolean phenotypes, against the 2 to the k possible, is reported as a diagnostic of threshold quality. Agreement between two arms is the adjusted Rand index between their full clusterings.

### 2.5 Cross-site variance

For the HIPC panel the measurement is the spread of each cluster's frequency across the seven centers, for a given lot, after averaging the three replicates within each center. It is the coefficient of variation, with its median taken over clusters that hold at least 1 percent of cells. Two denominators are reported, all events and the arm's own T cell clusters, those whose median CD3 lies above the mixture cut, because the all-events denominator inherits the lymphocyte gate and the debris fraction differs by center. The central manual gates give the reference, from their population counts per file, with the same two denominators. The replicate coefficient of variation within a center is the noise floor. The same analysis ran on the ten mice of Samusik_all as a contrast, where the samples really differ and the manual labels set the floor.

Per-cell manual labels for the HIPC files were produced by applying the central gating trees from the FlowJo 9 workspaces in Python, after the transform of section 2.6 was reproduced. Counts matched the published population statistics to the event at two centers and within 1.2 percent on CD3 at the others, and these labels gave Hungarian F1 on the HIPC panel.

### 2.6 The deposited HIPC transform

The HIPC T cell files deposited for the seven centers with central gates are flowCore exports from flowWorkspace, compensated and transformed, with the fluorescence columns reordered into the marker order of the HIPC scripts and the fluorochrome keyword left stale from the original column order. The README with the files describes an iplogicle transform (T 16777215, M 7.22, W 1.16). That transform does not reproduce the files. Stanford deposited its raw Diva exports beside the processed ones, with the same events in the same order, so the transform could be fitted event by event. It is a logicle [9] with T 262144, M 4.5 decades, no extra negative decades, scaled to 4096, and a width chosen per channel, which is FlowJo's biexponential as flowWorkspace reproduced it from each center's workspace; the compensation matrix is the one the analyst edited in FlowJo rather than the acquisition matrix. The fitted widths and the reconstruction error are in Results 3.6. Arm B on this panel inverts that transform per channel, divides, and applies it again. Only Stanford deposited raw files, so the other six centers take Stanford's widths, and the sensitivity of arm B to that assumption was measured. About 0.6 percent of events carry values below minus 4096, an artifact of the export, and were dropped.

### 2.7 Measured brightness

The stain index of each conjugate on the LSRFortessa was computed on the bead compensation controls, uncompensated, on the primary detector, as the median of the positive beads minus the median of the unstained beads, over twice the robust standard deviation of the unstained beads. Spearman's rank correlation between the measured index and the vendor rank is reported. On the HIPC panel a stain index per dye per center was computed the same way from each center's single-stain controls and used as the per-site divisor.

### 2.8 Software and reproducibility

Python 3.14 with flowsom, scikit-learn, flowio, flowutils and rdata. Every run is one script with the random seeds recorded; every table in this paper is written by a script from the metrics files; and the repository holds the plan, the per-run reports and the data notes with sources and checksums. The MIFlowCyt items for the DURAClone controls are in the data notes.

## 3. Results

### 3.1 Per-channel scaling costs accuracy

Figure 1 and Table 3 give mean F1 and adjusted Rand index for the three intensity arms at the free cluster number. Plain arcsinh intensity scored highest or within 0.01 of highest on every labeled set: 0.79 on Levine_32dim, 0.53 on Levine_13dim, 0.70 on Samusik_01 and 0.67 on Samusik_all. The Levine_32dim value matches the FlowSOM result in Weber and Robinson [1], which says the pipeline reproduces the benchmark.

Z-scoring cost 0.16, 0.12, 0.10 and 0.05 of mean F1 on the four sets, and 0.12, 0.20, 0.37 and 0.25 of adjusted Rand index. The z-scored clustering is a different partition of the same cells and not a small perturbation, since its agreement with the arcsinh clustering was 0.67, 0.53, 0.46 and 0.57. Percentile scaling cost less on mean F1, 0.15, 0.05, 0.04 and nothing on Samusik_all, and 0.10 to 0.28 on adjusted Rand index. At the gated cluster number the costs were smaller and in two cases reversed, because at a forced count the merge order decides which populations share a cluster; on Levine_13dim at 24 metaclusters arcsinh merged the naive and mature CD4 T cells and percentile scaling kept them apart.

On the two rare-population sets the cost fell on the population that mattered. On Nilsson_rare the best cluster for the HSC scored F1 0.50 under arcsinh, 0.38 under z-score and 0.23 under percentile. On Mosmann_rare the 109 activated memory CD4 T cells, 0.03 percent of the file, sit above the 99th percentile of every cytokine channel, so percentile clipping put them in one bin with the top 1 percent of negative cells, and F1 fell from 0.57 to 0.09. Z-scoring did not hurt that population (0.61).

The cofactor is the one preprocessing choice the field argues about, so the three intensity arms also ran across cofactors (Figure 6). On the mass cytometry sets the benchmark's value of 5 was the best choice for plain intensity, with mean F1 0.75, 0.79, 0.67 and 0.56 on Levine_32dim at cofactors 1, 5, 15 and 50, and the same shape on Levine_13dim and Samusik_01. The ranking of the three arms held at every cofactor, so the scaling cost is a property of the scaling and not of the cofactor it was measured at. The cofactor moved the result as much as scaling did, 0.79 to 0.56 from 5 to 50 on Levine_32dim. On the two flow sets the picture differed. The rare HSC on Nilsson_rare was recovered at 0.48 to 0.50 under plain intensity with cofactors of 50 and 150 and fell to 0.26 and 0.18 at 500 and 1,500, while under z-score and percentile scaling it rose to 0.54 and 0.61 at 500. The rare cytokine-positive cells on Mosmann_rare were recovered at 0.12 at a cofactor of 50 and at 0.69 at 1,500 under plain intensity, while z-score scaling recovered them at 0.43 at a cofactor of 50 and at 0.72 at 1,500. A large cofactor flattens the low range of a channel and leaves the clustering to the bright channels, and that is the condition under which scaling helps. Scaling is a correction for a cofactor that is wrong for the data, and at a cofactor that is right for the data it costs accuracy.

### 3.2 The loss falls on intermediate and small populations

Figure 2 gives per-population F1 for the populations the hypothesis named, across the intensity arms and the binary arms. Three patterns repeat.

A marker that the gates read at three levels cannot be cut once. On Levine_13dim the monocytes are labeled CD11b-negative, mid and high; the mid level scored 0.32 under arcsinh and 0.13 to 0.16 under the binary arms, while the high and negative levels held up. On Samusik_all the 82,384 intermediate monocytes, the second largest population in the set and defined by Ly6C at a middle level, scored 0.97 under arcsinh, 0.23 under z-score and 0.17 under lineage bits. Z-scoring Ly6C cost almost as much as binarizing it.

A population that differs from a larger neighbor by one marker is merged. On Levine_32dim, pro B cells (513 cells), plasma B cells (330) and CD34+CD38+CD123+ progenitors (304) each differ from a neighbor by one bit among the fifteen markers the labels define cleanly, and the binary arms merged them: 0.66 to 0.14, 0.79 to 0.01 and 0.55 to 0.16. Under arcsinh they differ by intensity on several markers at once.

A dim shift on several markers is what the rare populations are made of. The CD34+CD38lo HSCs on Levine_32dim fell from 0.74 under arcsinh to 0.09 under all-marker bits and to 0.43 under z-score.

### 3.3 Binary and reweighted inputs do not pass intensity

Table 4 gives the binary arms on the four mass cytometry sets at the free cluster number. With unsupervised mixture thresholds they scored 0.40 to 0.46 on Levine_32dim against 0.79, 0.25 to 0.32 on Levine_13dim against 0.53, 0.38 to 0.48 on Samusik_01 against 0.70 and 0.39 to 0.55 on Samusik_all against 0.67. The phenotype count says why. On Levine_32dim the mixture thresholds put 148,395 distinct phenotypes among 265,627 cells, so the binary data was nearly unique per cell and the self-organizing map had nothing to hold on to. The mixture fits its lower component to the zero spike of a mass cytometry channel and cuts inside the noise above it.

The question was then whether binary input fails because of the thresholds or because of the bits. With a threshold per marker chosen to maximize balanced accuracy against the labels, the upper bound, the binary arms on Levine_32dim scored 0.45 to 0.48, because the label set cannot define 17 of the 32 markers and those bits are noise under any cut. With only the 15 markers whose supervised cut reached balanced accuracy 0.9, binary input tied arcsinh at the gated cluster number (0.47 against 0.47, with a higher adjusted Rand index, 0.78 against 0.66) and reached 0.68 against 0.79 at the free number. That arm used the labels twice, to choose the markers and to place the cuts, and it still lost the three smallest populations. Leiden on a Hamming graph did not change the order. The phenotype count behaved as a clean panel should once the noise bits were removed, 3,286 occupied with 371 covering 95 percent of cells, which makes the count a usable diagnostic of threshold quality.

The seven-marker HIPC panel is the one conventional flow panel in the study where a bit per marker is not noise, because the cells occupy all 32 and 128 possible phenotypes, with 18 and 45 covering 95 percent of cells. The binary arms were still the least stable there and did not pass intensity on recovery.

Reweighting intensity rather than discretizing it did not pass intensity either. Weighting each channel by the separation of its mixture modes scored 0.68 against 0.79 on Levine_32dim and 0.66 against 0.70 on Samusik_01, and the weights went to the broadly expressed markers rather than to the ones that split populations. Weighting by the sharpness of the posterior put the weights where a cytometrist would, on CD3, CD4, CD8 and TCRb and away from the dim channels, and scored 0.75 against 0.79 and 0.69 against 0.70, with the small populations paying. The CD34+CD38+CD123+ progenitors fell from 0.55 to 0.34, because the graded markers that define them are the ones sharpness weights down. The posterior in place of the bit matched the bit. In every case the measured spread of a channel already carried the weighting the data supports, and a representation that changed it lost.

### 3.4 Brightness division changes membership and does not improve recovery

On the two conventional flow sets with known fluorochromes, dividing by conjugate brightness moved cells between clusters (Figure 3b). On Nilsson_rare the agreement between the arcsinh clustering and arm B was 0.69 and 0.59 for the two divisor rules, and between arm B and the z-score control 0.53; on Mosmann_rare 0.56 and 0.78 against arcsinh and 0.47 and 0.60 against z-score. Arm B is therefore not per-channel scaling in disguise, which was the first objection the study had to answer.

It also cost the rare population on both full panels: 0.50 to 0.22 on Nilsson_rare and 0.57 to 0.13 on Mosmann_rare, with recall intact and precision gone. The three cytokines that define the Mosmann population are on APC, PE and PE-Cy7, the three brightest dyes in the panel, and the division pulled them down by 2.6 while multiplying the nine dim channels by 1.6. The panel designer had put the bright dyes on the markers the question depended on. On the ten markers that Nilsson_rare shares with Levine_13dim, with the three channels that carried the cost removed, arm B and arcsinh recovered the HSC alike, 0.48 to 0.50 against 0.50.

The divisor was the weak point of arm B, and a measurement settled how weak. On the University of Florida LSRFortessa the bead stain index of the twelve DURAClone conjugates spanned from 83 to 4,265, a factor of about 50, where the 2-to-the-rank rule spans 4, and it had no rank correlation with the vendor ranks. Six of the twelve control tubes did not match the stored spillover matrix, which limits three of the values. A vendor rank is the wrong divisor, and a measured one belongs to one instrument on one day.

Mass cytometry has no dye, but its detector sensitivity varies with isotope mass by several fold [11], so it is not a clean control for the brightness term. With the published sensitivity curve of the first-generation CyTOF as a measured divisor on Levine_13dim, on the ten markers shared with Nilsson_rare, arm B moved membership (agreement 0.71 with arcsinh) and left the HSC where it was, at F1 0.07 under every arm, because the shared markers do not resolve that population in that file. A box on CD34, CD38 and CD90 around the labeled HSC is 17 percent pure there against 64 percent on the flow set. The matched pair repeats the finding on a platform without dye and does not sharpen it.

### 3.5 No representation removes the instrument

Figure 4 gives the between-site spread on the HIPC panel. The central manual gates themselves vary across the seven centers at a median coefficient of variation of 0.14 among T cells, with CD4 and CD8 at 0.04 and 0.06 and the memory subsets at 0.09 to 0.18, as Finak et al. reported [8]; relative to all events the median is 0.25, and most of that is the lymphocyte gate, because the debris fraction differs by center. The replicate spread within a center is 0.04.

At 12 metaclusters, the granularity of the 14 gated populations, plain arcsinh intensity matched the gates: 0.15 among T cells. Z-score and percentile scaling gave 0.17, the binary arms 0.17 to 0.19, and arm B 0.23 and 0.27. A divisor measured on each center's own single stains did not lower it (0.30 per center, 0.28 pooled). At 30 metaclusters every arm spread more, 0.26 to 0.50, because finer clusters cut the CD45RA, CCR7, CD38 and HLA-DR axes finer than the analyst did, and the positions along those axes shift by instrument; arm B with rank divisors was the widest at 0.50. The between-site component is three to five times the replicate noise under every representation, and none of them removes it. Hungarian F1 against the per-cell central gates was 0.37 to 0.42 across arms at 12 metaclusters and 0.55 to 0.74 at 30, the same per center as pooled, so the site split does not explain the low F1 at the coarse count.

The ten mice of Samusik_all show the other side of this measurement. There the samples really differ, the manual labels vary across mice at a median coefficient of variation of 0.19, and every arm that reports the biology sits at 0.20 to 0.21. The one arm that went lower, lineage bits at 0.15, was the least accurate, because its clusters mix populations whose frequencies move in different directions, and the mixture averages out. A lower between-sample spread is better only when the samples are the same, and the correlation across mice between a population's frequency and its matched cluster's frequency, 0.96 for arcsinh and 0.75 for lineage bits, orders the arms as F1 does.

### 3.6 The deposited HIPC files carry FlowJo's transform

Figure 5 gives the reconstruction. The iplogicle in the files' README leaves a residual that grows with the value, a ratio of 1.19 at the low end and 0.83 at the top on CD3. A logicle with T 262144, M 4.5 and no extra negative decades, scaled to 4096, reproduces the deposited values to a median error of 1 to 6 units when the width is fitted per channel: 0.40 on the viability dye, 0.90 on CCR7, 0.80 on CD4 and CD45RA, 0.15 on CD3, 0.75 on HLA-DR, 0.60 on CD38 and 1.10 on CD8, the same on two files from the same workspace. The compensation matrix that reproduces the files is the one the analyst edited in FlowJo, which at Stanford differs from the acquisition matrix in the FITC to PE cell (0.43 against 0.97) and the PE-Cy7 to APC-Cy7 cell (0.16 against 0.31). The marker keyword on each column is right and the fluorochrome keyword is stale. Among CD3-bright cells the columns named CD4 and CD8 exclude each other at r of minus 0.62 and minus 0.65, and no other pair does.

The widths moved arm B on Stanford by an adjusted Rand index of 0.80 against a flat width and 0.60 against doubled widths, while arm B against arcsinh was 0.82 and the seed-to-seed stability of arcsinh 0.80. On this panel, at this granularity, arm B is within seed noise of intensity, and for the six centers without raw files the inverse transform is as large an effect as the arm. The arm B numbers on the HIPC panel carry that caveat.

## 4. Discussion

The study set out to show that the fluorochrome is a hidden parameter of every clustering, and it found that the fixes for that parameter do more harm than the parameter. Every representation tested here removes or reweights the measured spread of a channel, and every one of them lost to plain intensity in the same place, on the populations that an intermediate level or a dim shift across several markers defines. The spread is where those populations live, and scaling or binarization cleans it away.

The practical result is the cost of the optional scaling step. Z-scoring before FlowSOM costs 0.05 to 0.16 of mean F1 and up to 0.25 of adjusted Rand index at the cluster number the field uses, on four benchmark sets, and percentile clipping erases a rare population that sits above the 99th percentile of its own channels. Many pipelines offer one of these as a default or a recommended option. The recommendation from this study is plain. Keep the arcsinh intensity, and treat scaling as a change to the input that needs the same justification as a change to the algorithm. This is the result that changes what an analyst does, and it is the one a short report could carry alone.

The brightness term is real and it is small. Dividing by conjugate brightness changes which cluster a cell joins, by an amount a reviewer can see and distinct from scaling, and it does not change whether a population is found. Where it changed recovery it did harm, because panel designers put the bright dyes on the markers that matter, and the division removes that weight. A reader who wants to make a clustering independent of the panel should know that the panel was designed to make the clustering good. The vendor rank is also the wrong number, because a measured stain index spans fifty-fold where the rank spans four and the two do not agree in order. That measurement is of one instrument on one day, which is the point.

On the multi-site question the paper has a clear answer and a limit. Across seven centers with central gates, plain intensity reproduces the manual gates' between-site spread and nothing goes lower, under either a panel-level or a per-instrument divisor. The spread is a property of where the populations sit on each instrument, and a representation moves it only if it moves those positions. Finer clustering doubles the spread for every arm, which means a between-site number without its cluster count cannot be read. The limit is the per-cell scoring. Hungarian F1 on the HIPC panel is low at the coarse count for every arm, and the per-center runs say that comes from the gates' resolution against FlowSOM's and does not come from the sites.

The study has limits of its own. The manual gates are thresholds, so a binary representation scored against them can win by construction; this paper leaned on the two measurements that do not use the gates, stability and between-site spread, and the binary arms lost on those too. Every full-data result is one seed, and the stability column says how far one seed moves a number. The cofactor sweep in 3.1 shows that the cofactor moves a result as much as scaling does, and that the two flow sets want different cofactors; the paper's main comparisons were run at the benchmark's values and the sweep bounds what another choice would change. The brightness divisors on the two public flow sets are vendor ranks with assumed values, and the one measured set is one instrument. The gate thresholds and the lineage splits are judgments made for each dataset and stated in the repository. No full-spectrum dataset with manual gates was available, so the panel on which the brightness argument would bite hardest, forty conjugates from PerCP to PE, is untested.

Two results stand apart from the hypothesis. The occupied phenotype count is a cheap diagnostic of whether a binary representation carries signal, 148,395 phenotypes among 265,627 cells against 3,286 with 371 covering most cells. And the HIPC T cell files on ImmuneSpace are not what their README says; the reconstruction in 3.6 is the one that reproduces them, and anyone who needs the linear scale from those files needs it.

## Acknowledgments

The DURAClone bead compensation controls were provided in kind by the University of Florida. The author thanks the authors of the benchmark datasets and the Human ImmunoPhenotyping Consortium for making their data and gates public, and Lukas Weber and Mark Robinson for the benchmark mirrors and HDCytoData.

## Funding

This work received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

## Competing interests

The author has been employed by BD Biosciences, which manufactures the Lyoplate reagents and the cytometers used in the HIPC study, and by Cellares. The author declares no other competing interests.

## Author contributions

H.S. conceived the study, designed and ran the analyses, and wrote the manuscript.

## Data and code availability

The benchmark datasets are public: Levine_32dim and Levine_13dim from the Weber and Robinson GitHub mirrors, Nilsson_rare, Mosmann_rare, Samusik_01 and Samusik_all from HDCytoData on Bioconductor ExperimentHub, all originally FlowRepository FR-FCM-ZZPH. The HIPC Lyoplate files are on ImmuneSpace, project HIPC/Lyoplate, with a free account. The DURAClone bead controls were provided in kind by the University of Florida and are available from the author on request; they are capture beads with no cells. The code, the per-run reports, the data notes with sources and checksums, the FlowJo 9 gate application, and the per-cell HIPC labels are in the repository at https://github.com/this-is-how-we-see/cytometry-intensity-spread.

## Tables

**Table 1. Datasets.** The six benchmark sets are FlowRepository FR-FCM-ZZPH as distributed by Weber and Robinson [1] and HDCytoData; the HIPC files are on ImmuneSpace, project HIPC/Lyoplate.

| Dataset | Cytometry | Sample | Cells | Labeled | Markers | Populations | Source |
|---|---|---|---|---|---|---|---|
| Levine_32dim | Mass (CyTOF) | Human bone marrow, 2 donors | 265,627 | 104,184 | 32 | 14 | Levine et al. [3]; Bendall et al. [5] |
| Levine_13dim | Mass (CyTOF) | Human bone marrow, 1 donor | 167,044 | 81,747 | 13 | 24 | Levine et al. [3]; Bendall et al. [5] |
| Samusik_01 | Mass (CyTOF) | Mouse bone marrow, 1 mouse | 86,864 | 53,173 | 39 | 24 | Samusik et al. [4] |
| Samusik_all | Mass (CyTOF) | Mouse bone marrow, 10 mice | 841,644 | 514,386 | 39 | 24 | Samusik et al. [4] |
| Nilsson_rare | Flow, 13 colors | Human bone marrow, CD34-enriched | 44,140 | 358 rare | 13 | 1 | Rundberg Nilsson et al. [6] |
| Mosmann_rare | Flow, BD LSR-II, 15 colors | Human PBMC, influenza-stimulated | 396,460 | 109 rare | 14 | 1 | Mosmann et al. [7] |
| HIPC Lyoplate T cell | Flow, BD LSR II and LSRFortessa, 7 centers | SeraCare PBMC, 3 lots, triplicates | 2,520,000 sampled from 63 files | central gates | 7 plus viability | 14 | Finak et al. [8] |
| DURAClone controls | Flow, BD LSRFortessa | Capture beads, one conjugate per tube | 12 tubes | none | 12 | none | Provided in kind by the University of Florida |

**Table 2. Representations.** Every arm starts from the compensated intensity x as deposited.

| Arm | Definition | What it removes | Distance |
|---|---|---|---|
| A1 | arcsinh(x / c), c = 5 for mass and 150 for fluorescence | Nothing; the control | Euclidean |
| A2 | A1, then each channel centered and divided by its SD | Differences in spread between channels | Euclidean |
| A3 | A1, then each channel scaled from its 1st to 99th percentile onto 0 to 1, clipped | Spread and the tails | Euclidean |
| B | arcsinh((x / b) / c), b the relative brightness of the conjugate (vendor rank or measured) | The conjugate term | Euclidean |
| C, C2 | Lineage markers 0 or 1 at the mixture cut; functional markers percentile (C) or z-score (C2) | Lineage intensity | Euclidean on the mixed matrix |
| D | Every marker 0 or 1 at the mixture cut | All intensity | Euclidean; Hamming for Leiden |
| W, W2 | A1 times a weight per marker from mode separation (W) or posterior sharpness (W2) | Nothing; reweights | Euclidean |
| P | C2 with the mixture posterior in place of each bit | Lineage intensity outside the band between modes | Euclidean |

**Table 3. Scaling cost.** FlowSOM, 10 by 10 grid, one seed on the full data; stability is the mean pairwise ARI over 10 seeds on 50,000 cells.

| Dataset | Metaclusters | Arm | Mean F1 | ARI | NMI | Stability |
|---|---|---|---|---|---|---|
| Levine_32dim | 14 | A1 arcsinh | 0.47 | 0.66 | 0.83 | 0.94 |
| Levine_32dim | 14 | A2 z-score | 0.43 | 0.63 | 0.81 | 0.90 |
| Levine_32dim | 14 | A3 percentile | 0.41 | 0.64 | 0.82 | 0.84 |
| Levine_32dim | 40 | A1 arcsinh | 0.79 | 0.80 | 0.83 | 0.77 |
| Levine_32dim | 40 | A2 z-score | 0.63 | 0.68 | 0.76 | 0.70 |
| Levine_32dim | 40 | A3 percentile | 0.64 | 0.52 | 0.74 | 0.64 |
| Levine_13dim | 24 | A1 arcsinh | 0.36 | 0.67 | 0.80 | 0.83 |
| Levine_13dim | 24 | A2 z-score | 0.33 | 0.65 | 0.78 | 0.82 |
| Levine_13dim | 24 | A3 percentile | 0.39 | 0.74 | 0.82 | 0.81 |
| Levine_13dim | 40 | A1 arcsinh | 0.53 | 0.79 | 0.83 | 0.70 |
| Levine_13dim | 40 | A2 z-score | 0.41 | 0.59 | 0.75 | 0.65 |
| Levine_13dim | 40 | A3 percentile | 0.48 | 0.61 | 0.77 | 0.69 |
| Samusik_01 | 24 | A1 arcsinh | 0.62 | 0.60 | 0.78 | 0.73 |
| Samusik_01 | 24 | A2 z-score | 0.54 | 0.56 | 0.73 | 0.78 |
| Samusik_01 | 24 | A3 percentile | 0.62 | 0.61 | 0.78 | 0.75 |
| Samusik_01 | 40 | A1 arcsinh | 0.70 | 0.85 | 0.83 | 0.75 |
| Samusik_01 | 40 | A2 z-score | 0.60 | 0.48 | 0.66 | 0.71 |
| Samusik_01 | 40 | A3 percentile | 0.66 | 0.60 | 0.76 | 0.67 |
| Samusik_all | 24 | A1 arcsinh | 0.58 | 0.64 | 0.78 | 0.72 |
| Samusik_all | 24 | A2 z-score | 0.59 | 0.66 | 0.77 | 0.80 |
| Samusik_all | 24 | A3 percentile | 0.62 | 0.69 | 0.81 | 0.75 |
| Samusik_all | 40 | A1 arcsinh | 0.67 | 0.88 | 0.85 | 0.77 |
| Samusik_all | 40 | A2 z-score | 0.62 | 0.63 | 0.75 | 0.75 |
| Samusik_all | 40 | A3 percentile | 0.68 | 0.78 | 0.81 | 0.69 |

| Rare-population set | Metaclusters | Arm | Best-cluster F1 | Precision | Recall | Stability |
|---|---|---|---|---|---|---|
| Nilsson_rare | 40 | A1 arcsinh | 0.50 | 0.36 | 0.84 | 0.54 |
| Nilsson_rare | 40 | A2 z-score | 0.38 | 0.25 | 0.87 | 0.50 |
| Nilsson_rare | 40 | A3 percentile | 0.23 | 0.13 | 0.94 | 0.52 |
| Mosmann_rare | 40 | A1 arcsinh | 0.57 | 0.40 | 0.98 | 0.52 |
| Mosmann_rare | 40 | A2 z-score | 0.61 | 0.44 | 0.99 | 0.57 |
| Mosmann_rare | 40 | A3 percentile | 0.09 | 0.05 | 0.72 | 0.59 |

**Table 4. Binary arms at the free cluster number (40).** Occupied phenotypes are the distinct bit patterns among all cells; possible is 2 to the number of bits.

| Dataset | Arm | Mean F1 | ARI | Stability | Occupied phenotypes | Possible | For 95 percent of cells |
|---|---|---|---|---|---|---|---|
| Levine_32dim | A1 arcsinh | 0.79 | 0.80 | 0.77 |  |  |  |
| Levine_32dim | C lineage bits + percentile | 0.40 | 0.22 | 0.39 | 112,988 | 33,554,432 | 99,707 |
| Levine_32dim | C2 lineage bits + z-score | 0.46 | 0.29 | 0.36 | 112,988 | 33,554,432 | 99,707 |
| Levine_32dim | D all bits | 0.41 | 0.24 | 0.40 | 148,395 | 4,294,967,296 | 135,114 |
| Levine_32dim | C2 with gate thresholds | 0.48 | 0.29 | 0.38 | 56,339 | 33,554,432 | 43,058 |
| Levine_32dim | D with gate thresholds | 0.45 | 0.29 | 0.41 | 120,773 | 4,294,967,296 | 107,492 |
| Levine_32dim | D with gate thresholds, 15 clean markers | 0.68 | 0.74 | 0.68 | 3,286 | 32,768 | 371 |
| Levine_13dim | A1 arcsinh | 0.53 | 0.79 | 0.70 |  |  |  |
| Levine_13dim | C lineage bits + percentile | 0.25 | 0.30 | 0.39 | 3,226 | 4,096 | 1,480 |
| Levine_13dim | C2 lineage bits + z-score | 0.32 | 0.37 | 0.45 | 3,226 | 4,096 | 1,480 |
| Levine_13dim | D all bits | 0.29 | 0.32 | 0.40 | 5,562 | 8,192 | 2,640 |
| Samusik_01 | A1 arcsinh | 0.70 | 0.85 | 0.75 |  |  |  |
| Samusik_01 | C lineage bits + percentile | 0.48 | 0.31 | 0.45 | 65,032 | 268,435,456 | 60,689 |
| Samusik_01 | C2 lineage bits + z-score | 0.38 | 0.38 | 0.55 | 65,032 | 268,435,456 | 60,689 |
| Samusik_01 | D all bits | 0.48 | 0.29 | 0.44 | 82,394 | 549,755,813,888 | 78,051 |
| Samusik_all | A1 arcsinh | 0.67 | 0.88 | 0.77 |  |  |  |
| Samusik_all | C2 lineage bits + z-score | 0.39 | 0.35 | 0.58 | 450,047 | 268,435,456 | 407,965 |
| Samusik_all | D all bits | 0.55 | 0.32 | 0.45 | 732,570 | 549,755,813,888 | 690,488 |

The gate-threshold arms use the labels to place the cuts, and the 15-marker arm also uses them to choose the markers; they are upper bounds and not methods.

## Figure legends

**Figure 1. Per-channel scaling costs accuracy.** (a) Recovery of the gated populations at the free cluster number under arcsinh intensity, z-score and percentile scaling: mean F1 by Hungarian matching for the four labeled mass cytometry sets and best-cluster F1 for the rare population in the two flow sets. (b) Adjusted Rand index for the labeled sets. FlowSOM, 10 by 10 grid, 40 metaclusters, one seed.

**Figure 2. Where the loss falls.** Per-population F1 at the free cluster number for the populations the hypothesis named: populations a marker defines at an intermediate level (CD11b on Levine_13dim monocytes, Ly6C on Samusik monocytes), populations one marker from a neighbor (pro B, plasma B and CD123-positive progenitors on Levine_32dim), and the rare populations of the two flow sets.

**Figure 3. Every representation, and brightness against scaling.** (a) Recovery under all seven arms per dataset. (b) Agreement, as adjusted Rand index, between the arcsinh clustering, the z-score clustering and the two brightness-divided clusterings on the two flow sets with known fluorochromes.

**Figure 4. Between-site spread on the HIPC Lyoplate T cell panel.** (a) Median coefficient of variation of cluster frequency across the seven centers, among T cells, for clusters at 1 percent or more, per representation at 12 metaclusters (full color) and 30 (light). The dotted line is the replicate spread within a center. (b) The same measurement for the central manual gates, per population, relative to CD3.

**Figure 5. The transform the deposited HIPC files carry.** (a) CD3 reconstructed from the Stanford raw file with the README's iplogicle against the deposited value. (b) The same with a logicle, T 262144, M 4.5, width 0.15. (c) The width fitted per channel.

**Figure 6. The cofactor.** Mean F1 under the three intensity arms across arcsinh cofactors, 1 to 50 on the mass cytometry sets and 50 to 1,500 on the flow sets, at the free cluster number.

## References

1. Weber LM, Robinson MD. Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data. Cytometry A. 2016;89:1084-1096. doi:10.1002/cyto.a.23030
2. Van Gassen S, Callebaut B, Van Helden MJ, et al. FlowSOM: using self-organizing maps for visualization and interpretation of cytometry data. Cytometry A. 2015;87:636-645. doi:10.1002/cyto.a.22625
3. Levine JH, Simonds EF, Bendall SC, et al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell. 2015;162:184-197. doi:10.1016/j.cell.2015.05.047
4. Samusik N, Good Z, Spitzer MH, Davis KL, Nolan GP. Automated mapping of phenotype space with single-cell data. Nat Methods. 2016;13:493-496. doi:10.1038/nmeth.3863
5. Bendall SC, Simonds EF, Qiu P, et al. Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. Science. 2011;332:687-696. doi:10.1126/science.1198704
6. Rundberg Nilsson A, Bryder D, Pronk CJ. Frequency determination of rare populations by flow cytometry: a hematopoietic stem cell perspective. Cytometry A. 2013;83:721-727. doi:10.1002/cyto.a.22324
7. Mosmann TR, Naim I, Rebhahn J, et al. SWIFT-scalable clustering for automated identification of rare cell populations in large, high-dimensional flow cytometry datasets, part 2: biological evaluation. Cytometry A. 2014;85:422-433. doi:10.1002/cyto.a.22445
8. Finak G, Langweiler M, Jaimes M, et al. Standardizing flow cytometry immunophenotyping analysis from the Human ImmunoPhenotyping Consortium. Sci Rep. 2016;6:20686. doi:10.1038/srep20686
9. Parks DR, Roederer M, Moore WA. A new "Logicle" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. Cytometry A. 2006;69:541-551. doi:10.1002/cyto.a.20258
10. Aghaeepour N, Chattopadhyay PK, Ganesan A, et al. Early immunologic correlates of HIV protection can be identified from computational analysis of complex multivariate T-cell flow cytometry assays. Bioinformatics. 2012;28:1009-1016. doi:10.1093/bioinformatics/bts082
11. Tricot S, Meyrand M, Sammicheli C, et al. Evaluating the efficiency of isotope transmission for improved panel design and a comparison of the detection sensitivities of mass cytometer instruments. Cytometry A. 2015;87:357-368. doi:10.1002/cyto.a.22648
12. Maecker HT, Frey T, Nomura LE, Trotter J. Selecting fluorochrome conjugates for maximum sensitivity. Cytometry A. 2004;62:169-173. doi:10.1002/cyto.a.20092
13. Nowicka M, Krieg C, Crowell HL, et al. CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets. F1000Res. 2019;6:748, version 3. doi:10.12688/f1000research.11622.3
14. Weber LM, Nowicka M, Soneson C, Robinson MD. diffcyt: differential discovery in high-dimensional cytometry via high-resolution clustering. Commun Biol. 2019;2:183. doi:10.1038/s42003-019-0415-5
15. Finak G, Perez JM, Weng A, Gottardo R. Optimizing transformations for automated, high throughput analysis of flow cytometry data. BMC Bioinformatics. 2010;11:546. doi:10.1186/1471-2105-11-546
16. Diggins KE, Greenplate AR, Leelatian N, Wogsland CE, Irish JM. Characterizing cell subsets using marker enrichment modeling. Nat Methods. 2017;14:275-278. doi:10.1038/nmeth.4149
17. Hahne F, LeMeur N, Brinkman RR, et al. flowCore: a Bioconductor package for high throughput flow cytometry. BMC Bioinformatics. 2009;10:106. doi:10.1186/1471-2105-10-106
18. Finak G, Jiang W, Gottardo R. CytoML for cross-platform cytometry data sharing. Cytometry A. 2018;93:1189-1196. doi:10.1002/cyto.a.23663
19. Hubert L, Arabie P. Comparing partitions. J Classif. 1985;2:193-218. doi:10.1007/BF01908075
20. Kuhn HW. The Hungarian method for the assignment problem. Naval Res Logist Q. 1955;2:83-97. doi:10.1002/nav.3800020109
