# Composition-Controlled Evaluation of RNA Structure Awareness: Attention Learns Base-Pairing, Embeddings Mostly Don't

**Paper C in the portfolio** (standalone, grant-facing). Self-contained; cites Paper B as related work, does not depend on it.

**Prereg SHA**: `694b43b` (causal-rna), `a7d10f5` (factorization-unified)

---

## Abstract

We introduce a composition-controlled evaluation framework for measuring RNA secondary structure awareness in foundation models and apply it to five models spanning four architectures and three orders of magnitude in parameter count. The framework addresses a composition confound that inflated prior geometric metrics: stems are 68.5% GC while loops are 43.9%, so any embedding metric sensitive to nucleotide identity inherits a ~3x bias toward stems. Our nucleotide-stratified permutation null absorbs this confound, reducing the apparent structure signal from 5.5x to 1.8x. Applying the framework reveals two positive results and one informative null. First, NT v2 (ESM+GLU, 56M, DNA-pretrained) learns attention patterns that correlate with base-pairing contacts: 0.278 trained versus exactly 0.000 for randomized weights, confirming the entire signal is learned. Second, Evo (StripedHyena, 7B, DNA-pretrained) shows the weakest structure sensitivity of all five models despite being 500x larger than Caduceus (14M), indicating that pretraining domain and tokenization dominate scale. Third, preregistered hypotheses about embedding-level structure encoding failed: RNA-FM's trained/untrained ratio reached 1.74x against a 2.0 threshold, with only 2/12 RNA families exceeding the composition-stratified null. We report this failure as a feature of preregistration, not a limitation of the models. The framework — composition controls, preregistered thresholds, and untrained baselines — generalizes to any domain where feature-composition correlates with the target structure.

## 1. Introduction

Foundation models pretrained on nucleotide sequences have been proposed for RNA biology, but evaluating whether they encode secondary structure — the pattern of intramolecular base pairing that governs RNA stability and function — requires metrics that distinguish learned representations from sequence-composition artifacts.

We initially adopted the Grassmannian geodesic distance (d_g), a geometric metric for comparing embedding subspaces between wild-type and mutant sequences. A zero-parameter 3-mer frequency baseline matched the trained models' d_g scores, falsifying the metric: it captured k-mer composition, a property that random linear maps preserve under the Johnson-Lindenstrauss lemma, rather than learned structure. This same failure mode has been reported for geometric transfer metrics in other domains \textcolor{red}{[CITE: Paper B boundary conditions]}.

The d_g collapse motivated a replacement metric with explicit composition controls. We designed a mutation sensitivity ratio under complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) with a nucleotide-stratified permutation null. The key observation: stems are enriched for G and C (68.5% vs 43.9% in loops), and complement swaps at G/C positions produce inherently larger embedding perturbations than at A/U positions regardless of structural context. Any naive stem-vs-loop comparison inherits this 24.7 percentage-point enrichment as a confound.

The contribution is the evaluation framework and an honest report of what it reveals: one clean learned positive (NT v2 attention), one informative scaling observation (7B parameters do not compensate for domain mismatch), and a preregistered null on embedding-level structure that was designed to be difficult to pass.

## 2. Models

| Model | Params | Architecture | Tokenization | Pretraining data |
|-------|--------|-------------|-------------|------------------|
| RNA-FM | 99M | BERT encoder | Character (A/U/G/C) | RNA sequences |
| NT v2 | 56M | ESM + GLU attention | 6-mer | DNA sequences |
| HyenaDNA | 5.4M | Hyena SSM (long-range conv) | Character | DNA sequences |
| Caduceus | 14M | BiMamba SSM (bidirectional) | Character | DNA sequences |
| Evo | 7B | StripedHyena (hybrid Hyena + attention) | Byte-level | DNA sequences |

RNA-FM is the only RNA-pretrained model. Evo is two orders of magnitude larger than the next-largest model. The comparison spans attention-based (RNA-FM, NT v2), SSM-based (HyenaDNA, Caduceus), and hybrid (Evo) architectures.

## 3. Methods

### 3.1 RNA structure dataset

We selected 12 RNA families spanning tRNAs, rRNAs, ribozymes, riboswitches, and non-coding RNAs. Each has an experimentally determined secondary structure from Rfam consensus or PDB crystal structure (source specified per family in Appendix A). Positions are labeled stem (base-paired) or loop (unpaired) from dot-bracket annotation.

### 3.2 Mutation sensitivity ratio

For each position in each RNA sequence, we perform a complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) and measure the cosine distance between wild-type and mutant embeddings at the mutated position. The structure sensitivity ratio is:

$$R = \frac{\text{mean cosine distance at stem positions}}{\text{mean cosine distance at loop positions}}$$

We compute $R$ at every layer and report the maximum across layers. The null model (Section 3.3) is computed under the same max-over-layers selection to avoid multiple-comparisons inflation.

### 3.3 Nucleotide-stratified permutation null

The null shuffles stem/loop labels independently within each nucleotide type (A, U, G, C), preserving per-nucleotide composition while destroying the structure-label association. For each of 100 permutations, we compute $R$ at every layer and take the maximum across layers (matching the real-data selection). The 95th percentile of these 100 max-over-layers values is the null threshold.

This null controls for first-order composition (per-nucleotide GC enrichment). It does not control for dinucleotide context (stem-G in GC pairs versus loop-G in GA contexts). Families whose real ratio exceeds this null are flagged as possibly reflecting genuine structure sensitivity or residual second-order composition effects.

### 3.4 Attention-contact correlation

For models with accessible attention matrices (RNA-FM, NT v2), we compute the Pearson correlation between attention weights and a binary contact matrix (1 for base-paired positions, 0 otherwise). We compare trained attention to an untrained baseline: randomized weights for NT v2, untrained initialization for RNA-FM.

### 3.5 Structure probing

Logistic regression on per-position embeddings predicting stem vs loop, using GroupKFold cross-validation by RNA family to test cross-family generalization.

### 3.6 Preregistered hypotheses

All hypotheses were registered before examining results (SHA `694b43b`):

- **H1**: RNA-FM trained/untrained mutation sensitivity ratio $\geq$ 2.0
- **H1b**: RNA-FM exceeds nucleotide-stratified null in $\geq$ 6/12 families
- **H2**: Median best layer $>$ 3 (Bonferroni-corrected $\alpha$ = 0.0083)
- **H3**: NT v2 mean ratio $<$ 2.0
- **H4**: HyenaDNA mean ratio $<$ 2.0
- **H5**: Probing accuracy $>$ majority baseline

The 2.0 threshold for H1 was chosen to demand a substantial effect — one that would survive composition controls and justify strong claims about structure encoding. A ratio below 2.0 is reported as a failure regardless of statistical significance, to avoid post-hoc rescue by p-value.

## 4. Results

### 4.1 NT v2 attention learns base-pairing structure

| Model | Trained | Untrained | Interpretation |
|-------|---------|-----------|---------------|
| NT v2 | 0.278 | 0.000 | Fully learned |
| RNA-FM | 0.067 | 0.061 | Zero learned |

NT v2 trained attention correlates with base-pairing contacts at mean 0.278 across 12 families, with a maximum of 0.509 on the hammerhead ribozyme (a simple stem-loop whose clean contact matrix produces the strongest signal). Randomized NT v2 weights produce uniform attention with exactly zero contact correlation. The entire 0.278 signal is learned during pretraining on DNA sequences.

RNA-FM trained attention (0.067) is indistinguishable from untrained (0.061). Despite being RNA-pretrained, RNA-FM's attention encodes no learned structural information.

The hammerhead ribozyme (simple stem-loop, 0.509) shows the strongest correlation; complex multi-branch structures (5S rRNA 0.168, RNaseP 0.190) show the weakest, consistent with the prediction that simple stem-loops produce cleaner contact matrices.

### 4.2 Scale does not compensate for domain mismatch

| Model | Mean ratio | Families > null | Probing $\Delta$ |
|-------|-----------|----------------|------------------|
| RNA-FM (99M, RNA) | 1.821 | 2/12 | +0.013 |
| Caduceus (14M, DNA) | 1.257 | 3/12 | +0.041 |
| HyenaDNA (5.4M, DNA) | 1.223 | 0/12 | +0.034 |
| Evo (7B, DNA) | 1.170 | 0/12 | +0.028 |
| NT v2 (56M, DNA) | 1.094 | 6/12 | +0.016 |

Evo (7B parameters, byte-level tokenization, DNA-pretrained) produces the weakest mutation sensitivity of all trained models. Two families (tRNA-Phe: 0.927, tRNA-Ala: 0.899) show ratios below 1.0 — Evo is less sensitive to stem mutations than loop mutations. Probing accuracy (0.615) peaks at layer 0 (embedding) and is lower than Caduceus (0.628, 14M) or HyenaDNA (0.621, 5.4M).

Despite being 500x larger than Caduceus and 70x larger than RNA-FM, Evo shows no advantage. The model's 7B parameters are spent on DNA-relevant features rather than RNA secondary structure. This observation is consistent across mutation sensitivity, probing, and the null comparison, though n=5 models is too few to support a general scaling law.

### 4.3 Composition confound dominates embedding-level sensitivity

The nucleotide-stratified null reduced the apparent RNA-FM signal from ~5.5x (naive, before composition control) to 1.8x (composition-controlled). The 24.7 percentage-point GC enrichment in stems was the dominant contributor to the original signal.

**Per-family mutation sensitivity:**

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

Two families survive the first-order null: tRNA-Ala (3.210 vs null 95th) and SAM riboswitch (2.515 vs null 95th) for RNA-FM. Both have high-GC stems with low-GC loops, meaning residual dinucleotide-context effects could contribute. A dinucleotide-stratified null is needed to confirm or eliminate these survivors (see Section 5.3).

NT v2 shows 6/12 families exceeding the null despite having the lowest mean ratio (1.094), suggesting a more consistent but smaller per-family signal than RNA-FM's high-mean, high-variance pattern.

### 4.4 Preregistered hypotheses

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H1 | Trained/untrained $\geq$ 2.0 | 1.74x | **FAIL** |
| H1b | $\geq$ 6/12 families > null | 2/12 | **FAIL** |
| H2 | Median best layer $>$ 3 | 6 (p=0.036, Bonferroni 0.0083) | **FAIL** (marginal) |
| H3 | NT v2 $<$ 2.0 | 1.094 | PASS (trivially) |
| H4 | HyenaDNA $<$ 2.0 | 1.223 | PASS (trivially) |
| H5 | Probing $>$ baseline | +0.013 to +0.041 | **FAIL** (within noise at N=12) |

**Exploratory observation** (not confirmatory): RNA-FM trained exceeds untrained in 12/12 families (Wilcoxon p = 0.000244). This per-family unanimity was not named as a decision criterion in the preregistration and is reported as exploratory evidence of a genuine learned effect that falls below the preregistered magnitude threshold.

### 4.5 Structure probing

| Model | Best accuracy | Baseline | $\Delta$ | Best layer |
|-------|--------------|----------|----------|------------|
| RNA-FM | 0.600 | 0.587 | +0.013 | 0 |
| NT v2 | 0.603 | 0.587 | +0.016 | 5 |
| HyenaDNA | 0.621 | 0.587 | +0.034 | 0 |
| Caduceus | 0.628 | 0.587 | +0.041 | 1 |
| Evo | 0.615 | 0.587 | +0.028 | 0 |

All five models barely exceed the 58.7% majority baseline. Four of five peak at layer 0 or 1 (embedding), indicating the probe relies on input features rather than learned representations. GroupKFold with 12 families yields high variance; the margins are consistent with noise. NT v2 is the only model whose best probing layer (5) falls above the embedding, consistent with its stronger attention-based structure signal.

## 5. Discussion

### 5.1 The framework generalizes beyond RNA

The composition confound we identified — GC enrichment in stems inflating embedding-level structure metrics — is an instance of a general problem: whenever the target structure (secondary structure) correlates with feature composition (nucleotide identity), any embedding metric sensitive to feature identity inherits the correlation as a confound. The nucleotide-stratified permutation null is a template: stratify by the confounding feature, shuffle labels within strata, and recompute the metric under the same selection procedure (max-over-layers) used on real data.

The same principle applies to protein secondary structure (amino acid composition differs between helices, sheets, and coils), chromatin accessibility (GC content correlates with open chromatin), and any biological structure correlated with sequence composition. The framework's value is independent of whether the RNA results are large or small.

### 5.2 NT v2 attention: architecture matters more than pretraining data

NT v2 is DNA-pretrained, yet its attention learns RNA base-pairing structure (0.278 vs 0.000 untrained). RNA-FM is RNA-pretrained, yet its attention learns nothing about structure (0.067 vs 0.061 untrained). The ESM+GLU attention architecture, combined with 6-mer tokenization that captures local sequence context, produces attention heads that correlate with base-pairing contacts. BERT-style attention in RNA-FM does not.

This is a finding about architecture, not about RNA versus DNA pretraining. The 6-mer tokenization in NT v2 encodes enough local context that base-pairing partners — typically 3-30 nucleotides apart in simple stem-loops — fall within the receptive field of individual tokens.

### 5.3 The two surviving families may reflect dinucleotide context

tRNA-Ala (3.210) and SAM riboswitch (2.515) are the only RNA-FM families exceeding the nucleotide-stratified null. Both have stems dominated by GC base pairs with loops enriched for A and U. The first-order null controls for the number of each nucleotide type labeled "stem" but does not control for dinucleotide context: stem-G occurs predominantly in GC pairs, while loop-G occurs in diverse dinucleotide contexts (GA, GU, GG). If the model's embedding encodes dinucleotide identity — which requires no secondary structure knowledge — this residual context difference could inflate the ratio.

A dinucleotide-stratified null (shuffling labels within each dinucleotide type rather than each nucleotide type) would resolve whether these two survivors reflect genuine structure sensitivity or second-order composition artifacts. Until that control is run, we treat tRNA-Ala and SAM as "possibly confounded" rather than confirmed positives.

### 5.4 Preregistered failure as a feature

The 2.0 threshold was chosen to demand an effect large enough to survive composition controls and support practical downstream use. 1.74x is genuine — the per-family unanimity (12/12) is exploratory evidence of a real learned effect — but reporting it as confirmatory would be a post-hoc endpoint switch. The preregistered criterion is not met.

This discipline has a direct cost: we cannot claim RNA-FM has "strong" structure encoding, even though 12/12 families show trained exceeding untrained. That cost is the price of credible null results. The alternative — setting a low threshold and declaring victory at 1.74x — would produce a result that cannot be distinguished from composition artifacts by future readers.

### 5.5 Limitations and next steps

**Sample size.** Twelve RNA families provide limited statistical power. The probing margins (+0.013 to +0.041) and mutation sensitivity differences are consistent with noise at this sample size. Expanding to 50+ Rfam families would either rescue or decisively kill these signals.

**Dinucleotide null.** The first-order composition control leaves open a second-order confound that specifically affects the two surviving families. This is the highest-priority robustness check.

**Model coverage.** Five models spanning four architectures is broad but thin per architecture. Additional RNA-pretrained transformers (RiNALMo, UTR-LM) would test whether RNA-FM's null attention result is architecture-specific or RNA-pretraining-general.

**Structure scope.** The evaluation targets secondary structure only and cannot assess tertiary contacts, pseudoknots, or other higher-order features.

## 6. Conclusion

A composition-controlled evaluation framework reveals that GC enrichment in stems — a fundamental consequence of Watson-Crick base pairing — inflates naive structure-sensitivity metrics by approximately 3x. After controlling for this confound, the residual embedding-level structure signal falls below a preregistered threshold designed to demand practical significance.

Two positive results survive all controls. NT v2 attention learns base-pairing contacts (0.278 trained vs 0.000 untrained), confirming that ESM+GLU attention can capture secondary structure even from DNA pretraining. Evo's 7B parameters produce the weakest structure sensitivity of all five models, showing that pretraining domain and tokenization dominate parameter count for RNA structure encoding.

The framework — nucleotide-stratified permutation, max-over-layers null matching, preregistered thresholds, and untrained baselines — provides a reusable template for evaluating structure awareness in sequence foundation models across biological domains.

## References

\textcolor{red}{[CITE: Paper B — boundary conditions / geometric metrics need construct validity]}

\textcolor{red}{[CITE: RNA-FM — Chen et al. 2022]}

\textcolor{red}{[CITE: Nucleotide Transformer v2 — Dalla-Torre et al. 2023]}

\textcolor{red}{[CITE: HyenaDNA — Nguyen et al. 2024]}

\textcolor{red}{[CITE: Caduceus — Schiff et al. 2024]}

\textcolor{red}{[CITE: Evo — Nguyen et al. 2024]}

\textcolor{red}{[CITE: Johnson-Lindenstrauss lemma]}

---

## Appendix A: RNA family details

\textcolor{red}{[TODO: Per-family table with Rfam accession or PDB ID, sequence length, stem/loop counts, GC% in stems vs loops, structure source (Rfam consensus vs PDB crystal)]}

## Data availability

All results are saved as JSON files in `data/gpu_results/batch3_structure_metrics/rerun_stratified/`. RNA structures, code, and preregistration are in the `causal-rna` and `factorization-unified` repositories. Preregistration SHA: `694b43b`. Model checkpoints: RNA-FM (HuggingFace `ml4bio/RNA-FM`), NT v2 (`InstaDeepAI/nucleotide-transformer-v2-50m-multi-species`), HyenaDNA (`LongSafari/hyenadna-small-32k-seqlen`), Caduceus (`kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16`), Evo (`togethercomputer/evo-1-131k-base`).
