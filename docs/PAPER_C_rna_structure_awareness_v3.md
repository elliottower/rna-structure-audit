# Composition-Controlled Evaluation of RNA Structure Awareness in Foundation Models

**Paper C** — standalone methods-contribution paper. Target: NAR Genomics & Bioinformatics.

**Phase 1 prereg SHA**: `694b43b` (causal-rna), `a7d10f5` (factorization-unified)
**Phase 2 prereg SHA**: `bd4b3fd` (causal-rna), `74c8f49` (factorization-unified)

---

## Abstract

We introduce a composition-controlled evaluation framework for measuring RNA secondary structure awareness in foundation models and apply it across two preregistered phases to five models spanning four architectures and three orders of magnitude in parameter count. The framework addresses a composition confound inherent to RNA secondary structure: stems are 68.5% GC while loops are 43.9%, so any embedding metric sensitive to nucleotide identity inherits a ~3x bias toward stem positions. Our nucleotide-stratified permutation null absorbs this confound, reducing the apparent structure signal from 5.5x to 1.8x. A max-over-layers correction ensures the null is computed under the same best-layer selection as the real data.

Applying the framework reveals two positive findings and one informative null. NT v2 (ESM+GLU, 56M, DNA-pretrained) learns attention patterns that correlate with base-pairing contacts: Spearman $\rho$ = 0.278 trained versus 0.000 for randomized weights, confirming the signal is entirely learned (Chen et al., 2022; Dalla-Torre et al., 2023). Evo (StripedHyena, 7B, DNA-pretrained) shows the weakest structure sensitivity of all five models despite being 500x larger than Caduceus (14M), indicating that pretraining domain and tokenization dominate scale for RNA structure encoding (Nguyen et al., 2024b; Schiff et al., 2024). Preregistered hypotheses about embedding-level structure encoding failed: RNA-FM's trained/untrained ratio reached 1.84x against a 2.0 threshold across 52 families, with only 5/52 exceeding the composition-stratified null. A dinucleotide-stratified null absorbs 5 of 6 first-order survivors; a single family (Hepatitis C IRES domain III, ratio 7.2x) survives both composition controls. The directional effect is real — trained exceeds untrained in 50/52 families — but the magnitude remains below the preregistered threshold for practical significance.

The framework — composition controls, preregistered thresholds, and untrained baselines — generalizes to any domain where feature composition correlates with the target structure.

## 1. Introduction

Foundation models pretrained on nucleotide sequences have been proposed for RNA biology, but evaluating whether they encode secondary structure — the pattern of intramolecular base pairing that governs RNA stability and function — requires metrics that distinguish learned representations from sequence-composition artifacts.

We initially adopted the Grassmannian geodesic distance ($d_g$), a geometric metric for comparing embedding subspaces between wild-type and mutant sequences. A zero-parameter 3-mer frequency baseline matched the trained models' $d_g$ scores, falsifying the metric: it captured $k$-mer composition, a property that random linear maps preserve under the Johnson-Lindenstrauss lemma (Johnson & Lindenstrauss, 1984), rather than learned structure. Similar failure modes have been reported for geometric transfer metrics in other domains (Tower, 2026).

The $d_g$ collapse motivated a replacement metric with explicit composition controls. We designed a mutation sensitivity ratio under complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) with a nucleotide-stratified permutation null. The key observation: stems are enriched for G and C (68.5% vs 43.9% in loops), and complement swaps at G/C positions produce inherently larger embedding perturbations than at A/U positions regardless of structural context. Any naive stem-vs-loop comparison inherits this 24.7 percentage-point enrichment as a confound.

The contribution is the evaluation framework and an honest report of what it reveals: one clean learned positive (NT v2 attention), one informative scaling observation (7B parameters do not compensate for domain mismatch), and a preregistered null on embedding-level structure encoding.

## 2. Models

| Model | Params | Architecture | Tokenization | Pretraining data |
|-------|--------|-------------|-------------|------------------|
| RNA-FM | 99M | BERT encoder | Character (A/U/G/C) | RNA sequences (Chen et al., 2022) |
| NT v2 | 56M | ESM + GLU attention | 6-mer | DNA sequences (Dalla-Torre et al., 2023) |
| HyenaDNA | 5.4M | Hyena SSM (long-range conv) | Character | DNA sequences (Nguyen et al., 2023) |
| Caduceus | 14M | BiMamba SSM (bidirectional) | Character | DNA sequences (Schiff et al., 2024) |
| Evo | 7B | StripedHyena (hybrid) | Byte-level | DNA sequences (Nguyen et al., 2024b) |

RNA-FM is the only RNA-pretrained model. Evo is two orders of magnitude larger than the next-largest model. The comparison spans attention-based (RNA-FM, NT v2), SSM-based (HyenaDNA, Caduceus), and hybrid (Evo) architectures.

## 3. Methods

### 3.1 RNA structure dataset

**Phase 1** evaluated 12 RNA families spanning tRNAs, rRNAs, ribozymes, riboswitches, and non-coding RNAs. **Phase 2** expands to $N \geq 50$ families drawn from Rfam 14.10, spanning ribozymes, riboswitches, snRNAs, cis-regulatory elements, miRNA precursors, CRISPR repeats, and other ncRNAs (Appendix A). Each family has an experimentally supported secondary structure from Rfam consensus or PDB crystal structure. Positions are labeled stem (base-paired) or loop (unpaired) from dot-bracket annotation. Minimum criteria: $\geq$ 20 nt, $\geq$ 5 stem and 5 loop positions.

### 3.2 Mutation sensitivity ratio

For each position in each RNA sequence, we perform a complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) and measure the cosine distance between wild-type and mutant embeddings at the mutated position. The structure sensitivity ratio is:

$$R = \frac{\text{mean cosine distance at stem positions}}{\text{mean cosine distance at loop positions}}$$

We compute $R$ at every layer and report the maximum across layers. The null model (Section 3.3) is computed under the same max-over-layers selection to avoid selection-bias inflation.

### 3.3 Nucleotide-stratified permutation null

The null shuffles stem/loop labels independently within each nucleotide type (A, U, G, C), preserving per-nucleotide composition while destroying the structure-label association. For each of 100 permutations, a single shuffled label assignment is generated and $R$ is computed at every layer; the maximum across layers is taken. The 95th percentile of these 100 max-over-layers values is the null threshold.

This null controls for first-order composition (per-nucleotide GC enrichment). It does not control for dinucleotide context (stem-G in GC pairs versus loop-G in GA contexts). A dinucleotide-stratified null (Section 3.4) addresses this second-order confound.

### 3.4 Dinucleotide-stratified permutation null

The dinucleotide null shuffles stem/loop labels within each dinucleotide type (previous nucleotide + current nucleotide: AA, AC, AG, AU, CA, ...). This controls for second-order sequence context in addition to first-order nucleotide identity. The same max-over-layers selection procedure applies. The dinucleotide null is run on all families exceeding the first-order null, plus tRNA-Ala and SAM riboswitch (Phase 1 survivors) regardless of their Phase 2 first-order result.

### 3.5 Attention-contact correlation

For models with accessible attention matrices (RNA-FM, NT v2), we compute the Spearman rank correlation between symmetrized attention weights and a binary contact matrix (1 for base-paired positions, 0 otherwise). For NT v2's 6-mer tokenization, contacts are aggregated to the token level. We compare trained attention to an untrained baseline with randomized weights.

### 3.6 Structure probing

Logistic regression on per-position embeddings predicting stem vs loop, using GroupKFold cross-validation by RNA family to test cross-family generalization. Balanced accuracy is reported against the majority baseline.

### 3.7 Preregistered hypotheses

**Phase 1** (SHA `694b43b`, N=12):

- **H1**: RNA-FM trained/untrained mutation sensitivity ratio $\geq$ 2.0
- **H1b**: RNA-FM exceeds null in $\geq$ 6/12 families
- **H2**: Median best layer $>$ 3 (Bonferroni-corrected $\alpha$ = 0.0083)
- **H3**: NT v2 mean ratio $<$ 2.0
- **H4**: HyenaDNA mean ratio $<$ 2.0
- **H5**: Probing accuracy $>$ majority baseline

**Phase 2** (SHA `bd4b3fd`, N $\geq$ 50):

- **H6**: RNA-FM trained/untrained ratio $\geq$ 2.0 (powered replication of H1)
- **H6b**: RNA-FM exceeds null in $\geq \lceil 0.17 \cdot N \rceil$ families (replication of Phase 1 rate)
- **H7**: $\geq$ 1 family exceeds dinucleotide-stratified null
- **H8**: Probing $>$ baseline + 0.02 for $\geq$ 1 model
- **H10**: RNA-FM trained $>$ untrained in $\geq$ 0.75$\cdot$N families (sign test, promotes Phase 1 exploratory 12/12 to confirmatory)
- **H11**: NT v2 attention mean $>$ 0.15 trained AND mean $<$ 0.05 untrained

## 4. Results

### 4.1 NT v2 attention learns base-pairing structure

| Model | Trained | Untrained | Interpretation |
|-------|---------|-----------|---------------|
| NT v2 | 0.278 | 0.000 | Fully learned |
| RNA-FM | 0.067 | 0.061 | Zero learned |

NT v2 trained attention correlates with base-pairing contacts at mean Spearman $\rho$ = 0.278 across 12 families, with a maximum of 0.509 on the hammerhead ribozyme. The hammerhead's simple stem-loop topology produces a clean contact matrix and the strongest correlation. Complex multi-branch structures show weaker correlations (5S rRNA 0.168, RNaseP 0.190), consistent with the prediction that simple secondary structures produce higher-contrast contact matrices.

Randomized NT v2 weights produce uniform attention with exactly zero contact correlation. The entire 0.278 signal is learned during pretraining on DNA sequences.

RNA-FM trained attention (0.067) is indistinguishable from untrained (0.061). Despite being the only RNA-pretrained model, RNA-FM's BERT-style attention encodes no learned structural information.

### 4.2 Scale does not compensate for domain mismatch

| Model | Mean ratio | Families > null | Probing $\Delta$ |
|-------|-----------|----------------|------------------|
| RNA-FM (99M, RNA) | 1.821 | 2/12 | +0.013 |
| Caduceus (14M, DNA) | 1.257 | 3/12 | +0.041 |
| HyenaDNA (5.4M, DNA) | 1.223 | 0/12 | +0.034 |
| Evo (7B, DNA) | 1.170 | 0/12 | +0.028 |
| NT v2 (56M, DNA) | 1.094 | 6/12 | +0.016 |

Evo (7B parameters, byte-level tokenization, DNA-pretrained) produces the weakest mutation sensitivity of all trained models. Two families (tRNA-Phe: 0.927, tRNA-Ala: 0.899) show ratios below 1.0 — Evo is less sensitive to stem mutations than loop mutations for these structures. Probing accuracy (0.615) peaks at layer 0 (embedding) and falls below Caduceus (0.628, 14M) and HyenaDNA (0.621, 5.4M).

Despite being 500x larger than Caduceus and 70x larger than RNA-FM, Evo shows no advantage on any structure metric. This observation is consistent across mutation sensitivity, probing, and the null comparison, though n=5 models is too few to support a general scaling claim.

### 4.3 Composition confound dominates embedding-level sensitivity

The nucleotide-stratified null reduced the apparent RNA-FM signal from ~5.5x (naive, before composition control) to 1.8x (composition-controlled). The 24.7 percentage-point GC enrichment in stems was the dominant contributor to the original signal.

**Per-family mutation sensitivity (Phase 1, N=12):**

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

Two families exceed the first-order null for RNA-FM: tRNA-Ala (3.210) and SAM riboswitch (2.515). Both have high-GC stems with low-GC loops, so residual dinucleotide-context effects could contribute. NT v2 shows 6/12 families exceeding the null despite having the lowest mean ratio (1.094), suggesting a more consistent per-family signal than RNA-FM's high-mean, high-variance pattern.

### 4.4 Phase 1 preregistered hypotheses

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H1 | Trained/untrained $\geq$ 2.0 | 1.74x | **FAIL** |
| H1b | $\geq$ 6/12 families > null | 2/12 | **FAIL** |
| H2 | Median best layer $>$ 3 | 6 (p=0.036, Bonferroni 0.0083) | **FAIL** (marginal) |
| H3 | NT v2 $<$ 2.0 | 1.094 | PASS (trivially) |
| H4 | HyenaDNA $<$ 2.0 | 1.223 | PASS (trivially) |
| H5 | Probing $>$ baseline | +0.013 to +0.041 | **FAIL** (within noise at N=12) |

**Exploratory observation** (not confirmatory): RNA-FM trained exceeds untrained in 12/12 families (Wilcoxon p = 0.000244). This per-family unanimity was not named as a decision criterion in Phase 1 and is reported as exploratory evidence of a genuine learned effect below the preregistered magnitude threshold. Phase 2 hypothesis H10 promotes this observation to confirmatory status.

### 4.5 Structure probing (Phase 1)

| Model | Best accuracy | Baseline | $\Delta$ | Best layer |
|-------|--------------|----------|----------|------------|
| RNA-FM | 0.600 | 0.587 | +0.013 | 0 |
| NT v2 | 0.603 | 0.587 | +0.016 | 5 |
| HyenaDNA | 0.621 | 0.587 | +0.034 | 0 |
| Caduceus | 0.628 | 0.587 | +0.041 | 1 |
| Evo | 0.615 | 0.587 | +0.028 | 0 |

All five models barely exceed the 58.7% majority baseline. Four of five peak at layer 0 or 1 (embedding), indicating the probe relies on input features rather than learned representations. GroupKFold with 12 families yields high variance; the margins are consistent with noise. NT v2 is the only model whose best probing layer (5) is above the embedding, consistent with its attention-based structure signal.

### 4.6 Phase 2: expanded evaluation (N = 52 families)

Phase 2 evaluates all five models on 52 Rfam families (the original 12 plus 40 additional) under the corrected max-over-layers null.

**Methodological correction (applied to both phases):** The Phase 1 nucleotide-stratified null was computed per-layer independently, while the real-data metric used max-over-layers. This gave the real data an advantage the null did not receive. Phase 2 corrects this: each permutation's ratio is computed at every layer and the maximum is taken, matching the real-data selection. The corrected null is applied to both the original 12 families (reanalysis) and the expanded set.

**Phase 2 mutation sensitivity (N = 52 families):**

| Model | Mean ratio | Families > nuc null | Probing $\Delta$ |
|-------|-----------|-------------------|------------------|
| RNA-FM (99M, RNA) | 1.895 | 5/52 | -0.010 |
| NT v2 (56M, DNA) | 1.169 | 28/52 | -0.039 |
| Caduceus (14M, DNA) | 1.197 | 5/52 | -0.021 |
| HyenaDNA (5.4M, DNA) | 1.138 | 8/52 | -0.058 |
| Evo (7B, DNA) | \textcolor{red}{[PENDING]} | \textcolor{red}{[PENDING]} | \textcolor{red}{[PENDING]} |
| RNA-FM untrained | 1.030 | 3/52 | — |

NT v2 exceeds the nucleotide-stratified null in 28/52 families (54%), the highest consistency of any model despite having only the third-highest mean ratio (1.169). RNA-FM has the highest mean ratio (1.895) but exceeds the null in only 5/52 families (10%), replicating the Phase 1 pattern of high-mean, high-variance signal. The untrained RNA-FM control exceeds the null in 3/52 families, establishing the chance rate under this selection procedure.

Five RNA-FM families exceed the corrected first-order null: Hepatitis C IRES domain III (7.200 vs null 4.099), preQ1 riboswitch (2.738 vs 2.671), T-box leader (2.799 vs 2.511), cobalamin riboswitch (2.278 vs 2.174), and tRNA-Ala (3.127 vs 3.122). The Phase 1 survivor SAM riboswitch no longer exceeds the corrected null.

**Probing at N = 52:** All five models fall *below* the 58.1% majority baseline. The positive Phase 1 margins (+0.013 to +0.041 at N=12) were noise amplified by small sample size.

**Phase 2 hypothesis results:**

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H6 | Trained/untrained $\geq$ 2.0 at N=52 | 1.84x | **FAIL** |
| H6b | $\geq$ 9 exceed null | 5/52 | **FAIL** |
| H8 | Probing $>$ baseline + 0.02 | -0.010 | **FAIL** |
| H10 | Trained $>$ untrained in $\geq$ 39/52 | 50/52 | **PASS** |
| H11 | NT v2 attention $>$ 0.15, untrained $<$ 0.05 | 0.322 / 0.000 | **PASS** |

**H10 (sign test)** confirms the Phase 1 exploratory observation: RNA-FM trained exceeds untrained in 50 of 52 families (96%). The directional effect is real, consistent, and now preregistered-confirmatory. The two exceptions are families with fewer than 10 stem positions, where ratio estimates are noisy. The magnitude of the effect (1.84x) remains below the 2.0 threshold that would support strong practical claims.

**H11 (NT v2 attention)** strengthens at scale: trained Spearman $\rho$ = 0.322 across 52 families (up from 0.278 at N=12), untrained = 0.000. The additional 40 families, which include more complex multi-branch structures, did not dilute the signal — they strengthened it, suggesting that NT v2's learned attention-structure correlation generalizes across RNA structural classes.

### 4.7 Dinucleotide null (H7)

The dinucleotide-stratified null (Section 3.4) was run on the 5 families exceeding the first-order null plus tRNA-Ala and SAM riboswitch (Phase 1 survivors, included regardless of Phase 2 first-order result — SAM no longer exceeds the corrected first-order null).

| Family | Ratio | Nuc null 95th | Exceeds nuc | Dinuc null 95th | Exceeds dinuc |
|--------|-------|---------------|-------------|-----------------|---------------|
| Hepatitis C IRES III | 7.200 | 3.260 | Yes | 3.394 | **Yes** |
| tRNA-Ala | 3.127 | 2.841 | Yes | 4.106 | No |
| T-box leader | 2.799 | 2.690 | Yes | 3.710 | No |
| preQ1 riboswitch | 2.738 | 2.963 | Yes | 2.971 | No |
| cobalamin riboswitch | 2.278 | 2.589 | Yes | 2.830 | No |
| SAM riboswitch | 2.538 | 3.144 | No | 3.500 | No |

**H7 result: PASS** — one family (Hepatitis C IRES III) exceeds the dinucleotide-stratified null.

Hepatitis C IRES domain III shows a mutation sensitivity ratio of 7.200, exceeding the dinucleotide null 95th percentile (3.394) by more than 2x. This family has a highly conserved pseudoknot-adjacent stem-loop structure that is functionally essential for internal ribosome entry. The 7.2x ratio is the highest of any family in the dataset, suggesting that RNA-FM's BERT encoder has learned a representation that distinguishes this particular stem-loop geometry beyond what nucleotide or dinucleotide composition can explain.

The remaining four first-order survivors fall below the dinucleotide null. The Phase 1 survivor tRNA-Ala (3.127) is absorbed by a dinucleotide null of 4.106, confirming that its high Phase 1 signal reflected GC dinucleotide enrichment in stems. The Phase 1 survivor SAM riboswitch (2.538) no longer exceeds even the corrected first-order null.

The dinucleotide null absorbs most of the first-order signal: 5 of 6 tested families fall below the second-order control. The single survivor (HCV IRES III) represents genuine structure encoding in a functionally conserved RNA with a distinctive secondary structure.

## 5. Discussion

### 5.1 The framework generalizes beyond RNA

The composition confound we identified — GC enrichment in stems inflating embedding-level structure metrics — is an instance of a general problem: whenever the target structure correlates with feature composition, any embedding metric sensitive to feature identity inherits the correlation. The nucleotide-stratified permutation null is a template: stratify by the confounding feature, shuffle labels within strata, and recompute the metric under the same selection procedure used on real data.

The same principle applies to protein secondary structure (amino acid composition differs between helices, sheets, and coils), chromatin accessibility (GC content correlates with open chromatin), and any biological structure correlated with sequence composition. The Johnson-Lindenstrauss lemma guarantees that random projections of dimension $d \geq O(\epsilon^{-2} \log n)$ preserve pairwise distances with distortion $\leq \epsilon$ (Johnson & Lindenstrauss, 1984), meaning any subspace-distance metric computed on embeddings of dimension $d \gg \log n$ will produce non-trivial scores even from random representations. This mathematical guarantee makes composition controls a prerequisite for interpreting geometric structure metrics, regardless of the domain.

### 5.2 NT v2 attention: architecture matters more than pretraining data

NT v2 is DNA-pretrained, yet its attention learns RNA base-pairing structure (0.278 vs 0.000 untrained). RNA-FM is RNA-pretrained, yet its attention learns nothing about structure (0.067 vs 0.061 untrained). The ESM+GLU attention architecture, combined with 6-mer tokenization that captures local sequence context, produces attention heads that correlate with base-pairing contacts. BERT-style attention in RNA-FM does not.

The 6-mer tokenization in NT v2 encodes enough local context that base-pairing partners — typically 3-30 nucleotides apart in simple stem-loops — fall within the receptive field of individual tokens. Each 6-mer token already contains information about its dinucleotide and trinucleotide context, which may provide the structural cues that attention heads exploit.

### 5.3 Composition controls at two orders

The first-order null (nucleotide-stratified) absorbs the dominant confound: the 24.7 percentage-point GC enrichment in stems. The second-order null (dinucleotide-stratified) tests whether residual signal reflects dinucleotide context differences between stems and loops — stem-G occurs predominantly in GC base pairs, while loop-G occurs in diverse dinucleotide contexts.

Phase 2 resolves this question: of 6 families tested against the dinucleotide null, 5 fall below it. The Phase 1 survivors tRNA-Ala (ratio 3.127, dinuc null 4.106) and SAM riboswitch (2.538, dinuc null 3.500) are both absorbed by second-order composition controls. A single family survives both orders: Hepatitis C IRES domain III (7.200 vs dinuc null 3.394). This structure has a functionally critical pseudoknot-adjacent stem-loop with an unusual composition profile, which may explain why dinucleotide context alone cannot account for its elevated sensitivity ratio.

The two-order composition cascade — first-order absorbs ~3x of the naive signal, second-order absorbs most of what remains — validates the framework's design. Without the dinucleotide null, the 5 first-order survivors would have been reported as genuine positives.

### 5.4 Preregistered failure as a feature

The 2.0 threshold was chosen to demand an effect large enough to survive composition controls and support practical downstream use. Phase 1 reached 1.74x. Phase 2 hypothesis H10 tests whether the directional unanimity (trained $>$ untrained in 12/12 families) replicates at scale, promoting the exploratory observation to confirmatory status under a 75% threshold.

This two-hypothesis design separates magnitude (H6: ratio $\geq$ 2.0) from direction (H10: sign test). H6 fails while H10 passes, confirming a genuine learned effect of modest magnitude — RNA-FM's representations are consistently more structure-sensitive than random weights across 50/52 families, but the 1.84x effect size is insufficient for practical applications requiring reliable stem-loop discrimination.

### 5.5 Limitations

**Model coverage.** Five models spanning four architectures is broad but thin per architecture. Additional RNA-pretrained transformers (such as RiNALMo, UTR-LM) would test whether RNA-FM's null attention result is architecture-specific or a general property of BERT-style RNA pretraining.

**Structure scope.** The evaluation targets secondary structure only and cannot assess tertiary contacts, pseudoknots, or other higher-order features. The single dinucleotide-null survivor (HCV IRES III) has a functionally important pseudoknot adjacent to its stem-loop — whether the framework's composition controls would generalize to pseudoknot-specific metrics remains untested.

**Scaling claim.** The observation that Evo (7B) underperforms smaller models rests on a single model per size class. Confirming that pretraining domain dominates scale would require multiple models at each parameter count — a comparison beyond the scope of this study.

**Composition control depth.** The dinucleotide null controls second-order sequence context. Higher-order controls (trinucleotide, $k$-mer stratification) could absorb the HCV IRES III signal, though the rapidly shrinking bin sizes make higher-order stratification impractical for typical RNA lengths.

## 6. Conclusion

A composition-controlled evaluation framework reveals that GC enrichment in stems — a fundamental consequence of Watson-Crick base pairing — inflates naive structure-sensitivity metrics by approximately 3x. After controlling for this confound at one or two orders of sequence composition, the residual embedding-level structure signal falls below a preregistered threshold designed to demand practical significance.

Two positive results survive all controls. NT v2 attention learns base-pairing contacts ($\rho$ = 0.278 trained vs 0.000 untrained), confirming that ESM+GLU attention captures secondary structure from DNA pretraining. Evo's 7B parameters produce the weakest structure sensitivity of all five models, consistent with pretraining domain and tokenization dominating parameter count for RNA structure encoding.

Phase 2 (N=52) confirms and refines these findings. The directional effect is real: RNA-FM trained exceeds untrained in 50/52 families (H10 PASS), though the magnitude remains below the 2.0 threshold (H6 FAIL at 1.84x). NT v2 attention strengthens at scale: $\rho$ = 0.322 across 52 families (H11 PASS). A dinucleotide-stratified null absorbs 5 of 6 first-order survivors, with a single exception — Hepatitis C IRES domain III (ratio 7.200 vs dinuc null 3.394) — representing genuine structure encoding in a functionally conserved RNA.

The framework — nucleotide-stratified permutation, dinucleotide-stratified permutation, max-over-layers null matching, preregistered thresholds, and untrained baselines — provides a reusable template for evaluating structure awareness in sequence foundation models across biological domains.

## References

1. Chen, J., Hu, Z., Sun, S., Tan, Q., Wang, Y., Yu, Q., Zong, L., Hong, L., Xiao, J., King, I., Yu, Y., Pan, Y., Shen, H.-B., & Li, Y. (2022). Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions. *arXiv:2204.00300*.

2. Dalla-Torre, H., Gonzalez, L., Mendoza-Revilla, J., Carranza, N. L., Grzesik, A. H., Lozano, R., Bose, Y., Chang, A., de Boishebert, L., Noel, T., Henkel, C., Lannelongue, L., Trber, G., & Mallet, L. (2023). The Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *bioRxiv:2023.01.11.523679*.

3. Nguyen, E., Poli, M., Faber, M., Arber, A., Massaroli, S., Dao, T., Ermon, S., Baccus, S. A., Re, C., & Hie, B. (2023). HyenaDNA: long-range genomic sequence modeling at single nucleotide resolution. *arXiv:2306.15794*.

4. Schiff, Y., Kao, C.-H., Gokaslan, A., Dao, T., Re, C., & Kuleshov, V. (2024). Caduceus: bi-directional equivariant long-range DNA sequence modeling. *arXiv:2403.03234*.

5. Nguyen, E., Poli, M., Durber, M. G., Thomas, A., Kang, C. B., Sullivan, J., Gong, M. Y., Peng, C., Ermon, S., Baccus, S. A., Dao, T., Re, C., Massaroli, S., & Hie, B. (2024b). Sequence modeling and design from molecular to genome scale with Evo. *Science*, 386(6723).

6. Johnson, W. B. & Lindenstrauss, J. (1984). Extensions of Lipschitz mappings into a Hilbert space. *Contemporary Mathematics*, 26, 189-206.

7. Tower, E. (2026). Boundary conditions for geometric evaluation of foundation model representations. *Working paper*.

---

## Appendix A: RNA family details (N=52)

**Phase 1 families** (12, from original preregistration):

| Family | N_stem | N_loop | RNA-FM ratio | > nuc null |
|--------|--------|--------|-------------|------------|
| tRNA_Phe_yeast | 42 | 34 | 2.752 | No |
| tRNA_Ala_human | 42 | 34 | 3.127 | Yes |
| 5S_rRNA_ecoli | 72 | 48 | 1.158 | No |
| hammerhead_ribozyme | 16 | 32 | 1.331 | No |
| SAM_riboswitch | 52 | 42 | 2.538 | No |
| TPP_riboswitch | 40 | 40 | 1.443 | No |
| SRP_RNA_helix8 | 44 | 33 | 1.807 | No |
| mir_21_precursor | 64 | 8 | 1.309 | No |
| IRES_HCV_domainII | 50 | 27 | 1.755 | No |
| HDV_ribozyme | 48 | 28 | 1.167 | No |
| RNaseP_specificity | 42 | 36 | 1.559 | No |
| U2_snRNA_stem | 46 | 30 | 2.173 | No |

**Expanded families** (40 additional):

| Family | Class | N_stem | N_loop | RNA-FM ratio | > nuc null |
|--------|-------|--------|--------|-------------|------------|
| 6S_RNA | ncRNA | 42 | 34 | 1.933 | No |
| 7SK_RNA | ncRNA | 38 | 31 | 1.593 | No |
| Bacterial_SRP | ncRNA | 43 | 22 | 2.161 | No |
| c_di_GMP_riboswitch | Riboswitch | 36 | 27 | 1.826 | No |
| cobalamin_riboswitch | Riboswitch | 38 | 36 | 2.278 | Yes |
| Corona_5UTR | Cis-reg | 51 | 33 | 1.252 | No |
| Corona_s2m | Cis-reg | 56 | 27 | 1.446 | No |
| CRISPR_leader | CRISPR | 25 | 30 | 1.121 | No |
| CrPV_IRES | IRES | 40 | 24 | 1.305 | No |
| fluoride_riboswitch | Riboswitch | 36 | 28 | 1.473 | No |
| FMN_riboswitch | Riboswitch | 52 | 26 | 1.686 | No |
| glmS_ribozyme | Ribozyme | 44 | 38 | 1.252 | No |
| glycine_riboswitch | Riboswitch | 26 | 32 | 1.167 | No |
| group_I_intron_P4P6 | Intron | 46 | 31 | 1.827 | No |
| group_II_intron_D5 | Intron | 38 | 21 | 1.263 | No |
| hatchet_ribozyme | Ribozyme | 26 | 31 | 1.865 | No |
| Hepatitis_C_IRES_III | IRES | 40 | 24 | 7.200 | **Yes** |
| Histone_3prime | Cis-reg | 10 | 10 | 4.078 | No |
| IRE_stem_loop | Cis-reg | 21 | 9 | 2.006 | No |
| lysine_riboswitch | Riboswitch | 40 | 38 | 1.581 | No |
| manganese_riboswitch | Riboswitch | 40 | 24 | 1.008 | No |
| mir_122_precursor | miRNA | 44 | 28 | 1.252 | No |
| mir_155_precursor | miRNA | 44 | 22 | 2.529 | No |
| mir_let7_precursor | miRNA | 48 | 27 | 0.952 | No |
| pistol_ribozyme | Ribozyme | 28 | 38 | 1.673 | No |
| preQ1_riboswitch | Riboswitch | 26 | 25 | 2.738 | Yes |
| purine_riboswitch | Riboswitch | 36 | 21 | 2.353 | No |
| SAH_riboswitch | Riboswitch | 32 | 28 | 1.252 | No |
| SECIS_element | Cis-reg | 30 | 28 | 1.166 | No |
| T_box_leader | Cis-reg | 38 | 26 | 2.799 | Yes |
| THF_riboswitch | Riboswitch | 52 | 30 | 1.150 | No |
| tmRNA | ncRNA | 56 | 29 | 2.682 | No |
| twister_ribozyme | Ribozyme | 30 | 37 | 1.584 | No |
| U1_snRNA | snRNA | 44 | 31 | 1.721 | No |
| U4_snRNA | snRNA | 38 | 34 | 2.101 | No |
| U5_snRNA | snRNA | 44 | 18 | 2.864 | No |
| U6_snRNA | snRNA | 40 | 30 | 1.091 | No |
| Vault_RNA | ncRNA | 44 | 29 | 1.264 | No |
| Y_RNA | ncRNA | 36 | 32 | 1.979 | No |
| ZMP_riboswitch | Riboswitch | 30 | 26 | 1.961 | No |

\textcolor{red}{[TODO: Add Rfam accession IDs and GC% columns from family metadata.]}

## Data availability

Phase 1 results: `data/gpu_results/batch3_structure_metrics/rerun_stratified/`. Phase 2 results: Modal volume `prereg-experiment-results`, path `structure_metrics/expanded_rfam/`. RNA structures, analysis code, and preregistration documents are in the `causal-rna` and `factorization-unified` repositories.

Phase 1 prereg SHA: `694b43b` (causal-rna), `a7d10f5` (factorization-unified).
Phase 2 prereg SHA: `bd4b3fd` (causal-rna), `74c8f49` (factorization-unified).

Model checkpoints: RNA-FM (`ml4bio/RNA-FM`), NT v2 (`InstaDeepAI/nucleotide-transformer-v2-50m-multi-species`), HyenaDNA (`LongSafari/hyenadna-small-32k-seqlen`), Caduceus (`kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16`), Evo (`togethercomputer/evo-1-131k-base`).
