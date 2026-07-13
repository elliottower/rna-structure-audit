# Do RNA Foundation Models Learn Secondary Structure? A Preregistered Cross-Architecture Comparison

**Paper C in the portfolio** (standalone, grant-facing). Cites Paper B for the theoretical principle that geometric metrics need construct validity checks.

**Prereg SHA**: `694b43b` (causal-rna), `a7d10f5` (factorization-unified)

---

## Abstract

We evaluate whether five RNA foundation models — spanning four architectures, three tokenization strategies, and three orders of magnitude in parameter count — encode secondary structure in their learned representations. An initial Grassmannian geodesic distance metric (d_g) collapsed under a zero-parameter 3-mer frequency baseline, revealing that the signal reflected nucleotide composition rather than learned biology. We designed a replacement metric — per-position mutation sensitivity under complement swap — with a nucleotide-stratified permutation null that controls for the 24.7 percentage-point GC enrichment in stems over loops. Across 12 RNA families, the composition control reduced the apparent structure signal from 5.5x to 1.8x. The preregistered decision criterion (trained/untrained ratio $\geq$ 2.0) was not met (1.74x). Among 60 model-family pairs, only 11 exceeded the nucleotide-stratified null. Evo (7B parameters, DNA-pretrained) showed the weakest sensitivity despite being 500x larger than Caduceus (14M), indicating that pretraining domain dominates scale. The sole strong positive result is NT v2 attention: trained attention-contact correlation of 0.278 versus exactly 0.000 for randomized weights, confirming that ESM+GLU attention learns base-pairing structure. The paper's methodological contribution is the composition-controlled evaluation framework and the honest report that preregistered hypotheses failed.

## 1. Introduction

RNA secondary structure — the pattern of intramolecular base pairing that folds a linear nucleotide chain into stems, loops, bulges, and junctions — governs RNA stability, function, and interaction with binding partners. Several foundation models pretrained on nucleotide sequences have been proposed for RNA biology, but whether these models encode structural information in their learned representations remains an open question.

Prior work has evaluated RNA foundation models on downstream benchmarks such as splice-site prediction and RNA family classification, but these tasks test general sequence discrimination rather than secondary structure encoding specifically. Geometric metrics such as Grassmannian subspace distance have been proposed as representation-quality measures across domains, but recent work has shown that such metrics require construct validity checks before deployment — they can be confounded by properties any random linear map preserves \textcolor{red}{[CITE: Paper B boundary conditions]}.

We designed a preregistered evaluation targeting structure sensitivity directly: does a model's representation change more when a mutation disrupts a base pair (stem position) than when it occurs in an unpaired region (loop position)? The critical design choice was a nucleotide-stratified permutation null, motivated by the observation that stems are 68.5\% GC while loops are 43.9\% GC. Complement swaps at G and C positions produce inherently larger embedding perturbations than at A and U positions, creating a composition confound that inflates naive stem-vs-loop comparisons.

## 2. Models

We compare five models spanning the major RNA/DNA foundation model architectures:

| Model | Params | Architecture | Tokenization | Pretraining data |
|-------|--------|-------------|-------------|------------------|
| RNA-FM | 99M | BERT encoder | Character (A/U/G/C) | RNA sequences |
| NT v2 | 56M | ESM + GLU attention | 6-mer | DNA sequences |
| HyenaDNA | 5.4M | Hyena SSM (long-range conv) | Character | DNA sequences |
| Caduceus | 14M | BiMamba SSM (bidirectional) | Character | DNA sequences |
| Evo | 7B | StripedHyena (hybrid Hyena + attention) | Byte-level | DNA sequences |

RNA-FM is the only RNA-pretrained model. The remaining four were pretrained on DNA. Evo is two orders of magnitude larger than the next-largest model.

## 3. Methods

### 3.1 RNA structure dataset

We selected 12 RNA families spanning tRNAs, rRNAs, ribozymes, riboswitches, and non-coding RNAs, each with an experimentally determined secondary structure from PDB or Rfam. Each position is labeled as stem (base-paired) or loop (unpaired) from the dot-bracket structure annotation.

### 3.2 Metric 1: Grassmannian geodesic distance (d_g) — falsified

We initially computed d_g between wild-type and mutant embedding subspaces as a structure-sensitivity metric. A zero-parameter 3-mer frequency baseline matched the trained models on d_g, falsifying the metric as a measure of learned structure. The same failure mode has been reported for Grassmannian subspace alignment in single-cell foundation models, where random projections reproduce learned-model scores \textcolor{red}{[CITE: Paper B]}.

### 3.3 Metric 2: Mutation sensitivity ratio (complement swap)

For each position in each RNA sequence, we perform a complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) and measure the cosine distance between the wild-type and mutant embeddings at the mutated position. The structure sensitivity ratio is the mean cosine distance at stem positions divided by the mean cosine distance at loop positions. A ratio greater than 1.0 indicates the model is more sensitive to structure-disrupting mutations (stems) than structure-neutral mutations (loops).

We report the best ratio across all layers for each model-family pair, as structure information may be encoded at different depths.

### 3.4 Nucleotide-stratified permutation null

Stems are enriched for G and C (68.5\% vs 43.9\% in loops). Because complement swaps at G/C positions produce larger embedding perturbations than at A/U positions regardless of structural context, a naive comparison inflates the stem-vs-loop ratio. Our null shuffles stem/loop labels independently within each nucleotide type (A, U, G, C), preserving per-nucleotide composition while destroying the structure-label association. We compute 100 permutations per model-family-layer and report the 95th percentile as the null threshold.

### 3.5 Attention-contact correlation

For models with accessible attention matrices (RNA-FM, NT v2), we compute the Pearson correlation between attention weights and a binary contact matrix derived from the secondary structure (1 for base-paired positions, 0 otherwise). We compare trained attention to untrained baselines — randomized weights for NT v2, untrained initialization for RNA-FM.

### 3.6 Structure probing

We train a logistic regression probe on per-position embeddings to predict stem vs loop labels, using GroupKFold cross-validation (by RNA family) to test cross-family generalization.

### 3.7 Preregistered hypotheses

All hypotheses were registered before examining results (SHA `694b43b`):

- **H1**: RNA-FM trained/untrained mutation sensitivity ratio $\geq$ 2.0
- **H1b**: RNA-FM exceeds nucleotide-stratified null in $\geq$ 6/12 families
- **H2**: Median best layer $>$ 3 (Bonferroni-corrected $\alpha$ = 0.0083)
- **H3**: NT v2 mean ratio $<$ 2.0
- **H4**: HyenaDNA mean ratio $<$ 2.0
- **H5**: Probing accuracy $>$ majority baseline

## 4. Results

### 4.1 d_g falsification

The Grassmannian distance metric was matched by a zero-parameter 3-mer frequency baseline, confirming that d_g captures sequence composition rather than learned structure. This motivated the switch to mutation sensitivity with composition controls.

### 4.2 Mutation sensitivity

| Family | RNA-FM | RNA-FM (untrained) | NT v2 | HyenaDNA | Caduceus | Evo |
|--------|--------|---------------------|-------|----------|----------|-----|
| tRNA_Phe_yeast | 2.704 | 1.050 | 1.180 | 1.185 | 1.374 | 0.927 |
| tRNA_Ala_human | 3.210 | 1.091 | 1.207 | 1.284 | 1.328 | 0.899 |
| 5S_rRNA_ecoli | 1.157 | 1.008 | 0.963 | 1.012 | 1.075 | 1.048 |
| hammerhead | 1.329 | 1.076 | 1.118 | 1.292 | 1.456 | 1.147 |
| SAM_riboswitch | 2.515 | 1.056 | 1.042 | 1.074 | 1.269 | 1.283 |
| TPP_riboswitch | 1.441 | 1.093 | 1.058 | 1.227 | 1.324 | 1.888 |
| SRP_RNA_helix8 | 1.350 | 1.035 | 1.168 | 1.029 | 1.136 | 0.925 |
| mir_21_precursor | 1.312 | 1.049 | 1.086 | 1.868 | 1.341 | 1.156 |
| IRES_HCV_domII | 1.804 | 1.045 | 0.928 | 1.190 | 1.225 | 0.994 |
| HDV_ribozyme | 1.318 | 1.021 | 1.175 | 1.150 | 1.291 | 1.316 |
| RNaseP_spec | 1.510 | 1.034 | 0.998 | 1.217 | 1.086 | 1.217 |
| U2_snRNA_stem | 2.200 | 1.024 | 1.203 | 1.144 | 1.173 | 1.234 |
| **Mean** | **1.821** | **1.048** | **1.094** | **1.223** | **1.257** | **1.170** |
| **Families > null** | **2/12** | **—** | **6/12** | **0/12** | **3/12** | **0/12** |

The composition control reduced the apparent RNA-FM signal from ~5.5x (naive) to 1.8x (composition-controlled). The nucleotide-stratified null thresholds are typically 1.1–1.7x, meaning the residual signal is modest relative to within-nucleotide variation.

### 4.3 Hypothesis results

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H1 | Trained/untrained $\geq$ 2.0 | 1.74x | **FAIL** |
| H1b | $\geq$ 6/12 families > null | 2/12 | **FAIL** |
| H2 | Median best layer $>$ 3 | 6 (p=0.036, Bonferroni 0.0083) | **FAIL** (marginal) |
| H3 | NT v2 $<$ 2.0 | 1.094 | PASS (trivially) |
| H4 | HyenaDNA $<$ 2.0 | 1.223 | PASS (trivially) |
| H5 | Probing $>$ baseline | +0.013 to +0.041 | **FAIL** (within noise at N=12) |

**Exploratory observation** (not confirmatory): RNA-FM trained exceeds untrained in 12/12 families (Wilcoxon p = 0.000244). This per-family unanimity was not named as a decision criterion in the preregistration, so it is reported as exploratory evidence of a genuine learned effect that falls below the preregistered magnitude threshold.

### 4.4 Attention-contact correlation

| Model | Trained | Untrained | Interpretation |
|-------|---------|-----------|---------------|
| NT v2 | 0.278 | 0.000 | Fully learned |
| RNA-FM | 0.067 | 0.061 | Zero learned |

NT v2 trained attention correlates with base-pairing contacts (mean 0.278 across 12 families, max 0.509 on hammerhead ribozyme). Randomized NT v2 weights produce uniform attention with zero contact correlation. The entire 0.278 signal is learned during pretraining. This is the strongest positive result.

RNA-FM trained attention (0.067) is indistinguishable from untrained (0.061). RNA-FM's attention encodes no learned structural information.

### 4.5 Structure probing

| Model | Best accuracy | Baseline | $\Delta$ | Best layer |
|-------|--------------|----------|----------|------------|
| RNA-FM | 0.600 | 0.587 | +0.013 | 0 |
| NT v2 | 0.603 | 0.587 | +0.016 | 5 |
| HyenaDNA | 0.621 | 0.587 | +0.034 | 0 |
| Caduceus | 0.628 | 0.587 | +0.041 | 1 |
| Evo | 0.615 | 0.587 | +0.028 | 0 |

All five models barely exceed the 58.7\% majority baseline. Four of five peak at layer 0 or 1 (embedding), indicating the probe relies on input features rather than learned representations. GroupKFold cross-validation with 12 families (2–3 per holdout fold) yields high variance; the margins are consistent with noise.

### 4.6 Scale vs pretraining domain

Evo (7B parameters, DNA-pretrained, byte tokenization) shows the weakest mutation sensitivity of all trained models: mean 1.170, 0/12 families exceed the null. Two families (tRNA-Phe: 0.927, tRNA-Ala: 0.899) show ratios below 1.0 — Evo is less sensitive to stem mutations than loop mutations. Despite being 500x larger than Caduceus and 70x larger than RNA-FM, Evo does not learn RNA secondary structure. Scale does not compensate for pretraining domain mismatch (DNA, not RNA) or tokenization mismatch (byte-level, destroying nucleotide identity).

## 5. Discussion

### 5.1 The composition confound was the dominant signal

The original d_g metric and the naive mutation sensitivity both captured GC enrichment in stems rather than learned structure. The nucleotide-stratified null absorbs this confound, leaving a residual signal that is real but modest (1.8x vs 5.5x). The 24.7 percentage-point GC enrichment in stems over loops is a fundamental property of RNA secondary structure — Watson-Crick base pairing is enriched for G-C pairs (three hydrogen bonds vs two for A-U) — and any embedding metric that is sensitive to nucleotide identity will inherit this enrichment as a confound.

This echoes findings in geometric transfer metrics for single-cell foundation models, where Johnson-Lindenstrauss preservation causes random projections to match learned models on subspace alignment and graph curvature metrics \textcolor{red}{[CITE: Paper B]}. The shared principle: geometric evaluation metrics require construct validity checks — null-model discrimination, composition controls, or untrained baselines — before deployment.

### 5.2 NT v2 attention as the positive result

NT v2 (ESM+GLU architecture, 6-mer tokenization) learns attention patterns that correlate with base-pairing contacts. The 0.278 trained vs 0.000 untrained comparison is the most conservative possible baseline — randomized weights eliminate all learned structure. The hammerhead ribozyme (simple stem-loop, 0.509) shows the strongest correlation, while complex multi-branch structures (5S rRNA, RNaseP) show weaker correlation, consistent with the prediction that simple stem-loops produce cleaner contact matrices.

### 5.3 Preregistered hypotheses failed — and that is the result

The preregistered criterion (trained/untrained $\geq$ 2.0) was chosen to demand a substantial, practically meaningful effect. 1.74x is genuine — 12/12 families show trained $>$ untrained — but falls below the threshold that would justify strong claims about structure encoding. Reporting the per-family unanimity as confirmatory evidence would be a post-hoc endpoint switch. The honest conclusion is that RNA-FM's structure sensitivity exceeds untrained weights, but the effect size does not meet the preregistered bar.

### 5.4 Limitations

The evaluation uses 12 RNA families, providing limited statistical power for cross-family generalization tests. Expanding to 50+ Rfam families would either confirm or kill the probing and mutation sensitivity signals. The nucleotide-stratified null controls for first-order composition (per-nucleotide GC enrichment) but does not control for dinucleotide context (stem-G in GC pairs vs loop-G in GA contexts), which could contribute to the two families (tRNA-Ala, SAM riboswitch) that survive the null. The evaluation targets secondary structure only, and cannot assess whether models encode tertiary contacts, pseudoknots, or other higher-order structural features.

## 6. Conclusion

Across five RNA foundation models, three evaluation approaches (mutation sensitivity, probing, attention analysis), and 12 RNA families, the strongest signal is a composition confound — GC enrichment in stems inflates stem-vs-loop comparisons by approximately 3x. After controlling for composition, the residual structure-dependent signal is real but modest, falling below preregistered thresholds. The sole unambiguous positive result is NT v2 attention-contact correlation, which is entirely learned (0.278 trained vs 0.000 untrained). The paper's contribution is methodological: a composition-controlled evaluation framework for RNA structure awareness, an honest report of preregistered hypothesis failure, and the finding that model scale does not compensate for domain mismatch in structure sensitivity.

## References

\textcolor{red}{[CITE: Paper B — boundary conditions / geometric metrics need construct validity]}

\textcolor{red}{[CITE: RNA-FM — Chen et al. 2022]}

\textcolor{red}{[CITE: Nucleotide Transformer v2 — Dalla-Torre et al. 2023]}

\textcolor{red}{[CITE: HyenaDNA — Nguyen et al. 2024]}

\textcolor{red}{[CITE: Caduceus — Schiff et al. 2024]}

\textcolor{red}{[CITE: Evo — Nguyen et al. 2024]}

\textcolor{red}{[CITE: Johnson-Lindenstrauss lemma]}

---

## Data availability

All results are saved as JSON files in `data/gpu_results/batch3_structure_metrics/rerun_stratified/`. RNA structures and code are in the `causal-rna` and `factorization-unified` repositories. Preregistration SHA: `694b43b`.
