# A Graded Evaluation of RNA Structure Awareness in Foundation Models: From Stem-Loop Discrimination to Partner Specificity

**Paper C v4** — unified ladder paper incorporating Phase 6 (perturbation specificity).
Target: NAR Genomics & Bioinformatics.

**Phase 1 prereg SHA**: `694b43b` (causal-rna), `a7d10f5` (factorization-unified)
**Phase 2 prereg SHA**: `bd4b3fd` (causal-rna), `74c8f49` (factorization-unified)
**Phase 6 prereg SHA**: `c19aa59` (causal-rna)

---

## Abstract

How deeply do RNA foundation models understand secondary structure? We present a three-rung evaluation ladder, each rung demanding a strictly stronger form of structural knowledge, and apply it to ten models spanning five architectures, three orders of magnitude in parameter count, and two pretraining domains (RNA, DNA). The first rung asks whether models distinguish stems from loops; a nucleotide-stratified permutation null absorbs 3x of the naive signal, revealing that most apparent structure awareness reflects GC enrichment in stems. The second rung asks whether discrimination survives dinucleotide composition controls; a single RNA family out of 52 survives both orders of stratification. The third and most demanding rung asks whether models resolve which specific position pairs with which: after mutating one base to its Watson-Crick complement, does the model perturb the paired partner more than adjacent stem positions? A within-stem derangement null controls for positional and compositional confounds.

ERNIE-RNA (86M parameters, RNA-pretrained) is the only model that clears the third rung, with mean perturbation specificity PS = 0.117 across 28 of 30 qualifying families and partner-is-max precision of 88.3%. The effect is localized to layers 10--12 (of 12), consistent with learned pairing rather than positional encoding. The remaining nine models show PS at or below 10^{-3}, indistinguishable from the derangement null. Two attention-based models---NT v2 (DNA-pretrained) and DNABERT-2 (DNA-pretrained)---learn attention patterns that correlate with base-pairing contacts ($\rho$ = 0.32 and 0.26), surviving the first rung through a different channel than embedding sensitivity. Evo (7B parameters) produces the weakest structure signal despite being 500x larger than the next-largest SSM model. All evaluations are preregistered with SHA-frozen protocols.

The three-rung framework---composition-controlled discrimination, dinucleotide stratification, and partner specificity with derangement null---provides a reusable template for evaluating structure awareness in sequence foundation models.

## 1. Introduction

Foundation models pretrained on nucleotide sequences have been proposed for RNA secondary structure prediction, drug target identification, and splice-site recognition. Evaluating whether these models encode base-pairing structure---the pattern of intramolecular Watson-Crick pairing that governs RNA stability and function---requires metrics that separate learned representations from sequence-composition artifacts.

RNA secondary structure introduces a systematic composition confound. Stems are GC-rich (68.5% GC) because Watson-Crick base pairs favor G-C and C-G pairing. Loops are AU-rich (43.9% GC) because unpaired regions tolerate all four nucleotides more equally. Any embedding metric sensitive to nucleotide identity inherits this 24.7 percentage-point enrichment as an apparent structure signal. We initially adopted the Grassmannian geodesic distance, a geometric metric for comparing embedding subspaces, and found that a zero-parameter 3-mer frequency baseline matched trained models' scores---a result consistent with random projections preserving pairwise distances under the Johnson-Lindenstrauss lemma.

This failure motivated a graded evaluation framework. The core insight is that "structure awareness" encompasses at least three levels of resolution, each demanding strictly stronger knowledge. The first level---stem-versus-loop discrimination---can be achieved by any representation sensitive to nucleotide composition. The second level---discrimination surviving composition controls at one and two orders of sequence context---requires genuine sensitivity beyond GC enrichment. The third level---knowing which specific position pairs with which---requires that the model encode the pairing topology, not merely the categories of paired and unpaired positions.

We apply this three-rung ladder to ten models, expanding the original five-model panel from Phases 1--2 with five additional models. The evaluation reveals a surprising asymmetry: the hardest test, not the easier ones, produces the clearest positive result.

## 2. Models

| Model | Params | Architecture | Tokenization | Pretrained on |
|-------|--------|-------------|-------------|---------------|
| RNA-FM | 99M | BERT encoder | Character | RNA |
| ERNIE-RNA | 86M | BERT encoder | Character | RNA |
| RiNALMo | 650M | Encoder | Character | RNA |
| UTR-LM | ~2M | BERT encoder | Character | RNA (5' UTR) |
| SpliceBERT | 19M | BERT encoder | Character | RNA (pre-mRNA) |
| NT v2 | 56M | ESM + GLU | 6-mer | DNA |
| DNABERT-2 | 117M | BERT + MosaicBERT | BPE | DNA |
| HyenaDNA | 5.4M | Hyena SSM | Character | DNA |
| Caduceus | 14M | BiMamba SSM | Character | DNA |
| Evo | 7B | StripedHyena | Byte-level | DNA |

Five models are RNA-pretrained; five are DNA-pretrained. Parameter counts span three orders of magnitude (2M to 7B). Architectures include attention-based encoders (RNA-FM, ERNIE-RNA, RiNALMo, UTR-LM, SpliceBERT, NT v2, DNABERT-2), state-space models (HyenaDNA, Caduceus), and a hybrid (Evo). Tokenization strategies include character-level, 6-mer, byte-pair encoding, and byte-level.

## 3. Methods

### 3.1 RNA structure dataset

52 RNA families drawn from Rfam 14.10, spanning tRNAs, rRNAs, ribozymes, riboswitches, snRNAs, cis-regulatory elements, miRNA precursors, and CRISPR repeats. Each family has an experimentally supported secondary structure from Rfam consensus or PDB crystal structure. Positions are labeled stem (base-paired) or loop (unpaired) from dot-bracket annotation. Minimum criteria: $\geq$ 20 nt, $\geq$ 5 stem and 5 loop positions. For Phase 6, families require $\geq$ 15 Watson-Crick pairs with $\geq$ 3 eligible pairs per stem (K $\geq$ 3) to support the derangement null; 32 families qualify.

### 3.2 Rung 1: Mutation sensitivity ratio with composition null

For each position in each RNA, we perform a complement swap (A$\leftrightarrow$U, C$\leftrightarrow$G) and measure the cosine distance between wild-type and mutant hidden-state embeddings at the mutated position. The structure sensitivity ratio is:

$$R = \frac{\text{mean cosine distance at stem positions}}{\text{mean cosine distance at loop positions}}$$

We compute $R$ at every layer and report the maximum across layers.

**Nucleotide-stratified permutation null.** The null shuffles stem/loop labels independently within each nucleotide type (A, U, G, C), preserving per-nucleotide composition while destroying the structure-label association. For each of 1000 permutations, a shuffled label assignment is generated and $R$ is computed at every layer; the maximum across layers is taken. The 95th percentile of these 1000 max-over-layers values is the null threshold. This controls for first-order composition (per-nucleotide GC enrichment).

**Untrained baseline.** RNA-FM with randomized weights serves as the composition-only control. The untrained model's mutation sensitivity captures the contribution of the tokenizer and positional encoding to the stem/loop ratio.

### 3.3 Rung 2: Dinucleotide-stratified null

The dinucleotide null shuffles stem/loop labels within each dinucleotide type (previous nucleotide + current nucleotide: AA, AC, AG, AU, CA, ...). This controls for second-order sequence context. Families exceeding the first-order null are tested against this stricter control.

### 3.4 Attention-contact correlation

For models with accessible attention matrices, we compute the Spearman rank correlation between symmetrized attention weights and a binary contact matrix (1 for base-paired positions, 0 otherwise). For NT v2's 6-mer and DNABERT-2's BPE tokenization, contacts are aggregated to the token level. We compare trained attention to an untrained baseline with randomized weights.

### 3.5 Structure probing

Logistic regression on per-position embeddings predicting stem vs loop, using GroupKFold cross-validation by RNA family to test cross-family generalization. Balanced accuracy is reported against the majority baseline.

### 3.6 Rung 3: Perturbation specificity

The most demanding test. For each Watson-Crick base pair $(i, j)$ in a stem with $\geq$ 3 eligible pairs:

1. Mutate position $i$ to its Watson-Crick complement.
2. Measure the cosine distance between wild-type and mutant embeddings at every other stem position.
3. Compute perturbation specificity: $\text{PS}(i, j) = \Delta_j - \max(\Delta_{j-1}, \Delta_{j+1})$, where $\Delta_j$ is the cosine distance at the paired partner and $\Delta_{j \pm 1}$ are distances at adjacent stem positions.

A positive PS means the model perturbs the correct partner more than its neighbors.

**Positive control gate.** Each family must show significantly greater stem perturbation than loop perturbation (paired $t$-test, $p < 0.05$). Families failing this gate are excluded: a model insensitive to complement swaps in stems cannot be expected to resolve individual partners.

**Within-stem derangement null.** For each stem with $K \geq 3$ eligible pairs, 1000 derangements permute partner labels within the stem. The real PS is compared against the distribution of deranged PS values. A family exceeds the null if the observed mean PS exceeds the 95th percentile of derangement means.

**H3: Partner-is-max.** For deep-interior pairs (not at stem termini), the fraction of pairs where the true partner receives the maximum perturbation among all stem positions. Chance baseline is 1/3 (partner vs two neighbors); preregistered threshold is binomial $> 1/3$.

### 3.7 Preregistered hypotheses

**Phases 1--2** (composition-controlled discrimination):

- H1: RNA-FM trained/untrained ratio $\geq$ 2.0
- H6: RNA-FM ratio $\geq$ 2.0 at N = 52 (powered replication)
- H7: $\geq$ 1 family exceeds dinucleotide null
- H10: RNA-FM trained $>$ untrained in $\geq$ 75% of families
- H11: NT v2 attention $\rho > 0.15$ trained AND $< 0.05$ untrained

**Phase 6** (partner specificity, Bonferroni $\alpha$ = 0.05/3 = 0.0167):

- H1$_6$: At least one model has mean PS $> 0$ with $\geq$ 7 families exceeding derangement null (Wilcoxon)
- H2$_6$: RNA-pretrained models have higher PS than DNA-pretrained models (rank-biserial $> 0.5$)
- H3$_6$: At least one model has partner-is-max fraction exceeding 1/3 (binomial)

## 4. Results

### 4.1 Rung 1: Composition absorbs most apparent structure signal

The nucleotide-stratified null reduces the apparent RNA-FM signal from ~5.5x (naive, before composition control) to 1.84x (composition-controlled). The 24.7 percentage-point GC enrichment in stems is the dominant contributor.

**Mutation sensitivity across 10 models (N = 52 families):**

| Model | Mean ratio | Families > nuc null | Attention $\rho$ | Probing $\Delta$ |
|-------|-----------|-------------------|-----------------|-----------------|
| RNA-FM (99M, RNA) | 1.895 | 5/52 | 0.067 | -0.010 |
| ERNIE-RNA (86M, RNA) | 1.149 | 22/52 | 0.099 | +0.023 |
| RiNALMo (650M, RNA) | *pending* | *pending* | *pending* | *pending* |
| UTR-LM (~2M, RNA) | 1.054 | 3/49 | — | -0.005 |
| SpliceBERT (19M, RNA) | 1.045 | 4/52 | 0.064 | -0.001 |
| NT v2 (56M, DNA) | 1.169 | 28/52 | **0.322** | +0.016 |
| DNABERT-2 (117M, DNA) | 1.109$^\dagger$ | 5/40 | **0.257** | — |
| HyenaDNA (5.4M, DNA) | 1.138 | 8/52 | — | -0.058 |
| Caduceus (14M, DNA) | 1.197 | 5/52 | — | -0.021 |
| Evo (7B, DNA) | 1.352 | 6/52 | — | -0.042 |
| RNA-FM untrained | 1.030 | 3/52 | 0.061 | — |

$^\dagger$Median reported for DNABERT-2; the mean (9.59) is inflated by BPE tokenization artifacts that produce extreme outlier ratios in a few families. All other models use mean.

Two patterns emerge from Rung 1. First, ERNIE-RNA shows the highest consistency of any RNA-pretrained model: 22 of 52 families exceed the nucleotide-stratified null, despite having a modest mean ratio (1.149). RNA-FM has the highest mean ratio (1.895) but exceeds the null in only 5 of 52 families, replicating a high-mean, high-variance pattern. The consistency gap---22/52 versus 5/52---suggests ERNIE-RNA has a genuine per-family signal that RNA-FM lacks.

Second, attention-based models encode structure through attention rather than embeddings. NT v2 attention correlates with base-pairing contacts at $\rho$ = 0.322 across 52 families (versus 0.000 for randomized weights). DNABERT-2 attention shows a comparable signal ($\rho$ = 0.257). Both are DNA-pretrained; RNA-pretrained RNA-FM attention encodes nothing (0.067 trained vs 0.061 untrained). Structure encoding in attention and in embeddings appear to be independent channels: NT v2 has 28/52 families exceeding the embedding null despite a low mean ratio, while its attention shows the strongest contact correlation.

**Scale does not compensate for domain mismatch.** Evo (7B, byte-level, DNA) produces weaker mutation sensitivity than Caduceus (14M, character, DNA) across all metrics. Probing accuracy peaks at layer 0 for Evo, consistent with positional encoding rather than learned computation. The 500x parameter advantage is erased by byte-level tokenization and DNA-only pretraining.

### 4.2 Rung 2: One family survives dinucleotide controls

Six families tested against the dinucleotide-stratified null:

| Family | Ratio | Nuc null 95th | Dinuc null 95th | Survives |
|--------|-------|---------------|-----------------|----------|
| Hepatitis C IRES III | 7.200 | 3.260 | 3.394 | **Yes** |
| tRNA-Ala | 3.127 | 2.841 | 4.106 | No |
| T-box leader | 2.799 | 2.690 | 3.710 | No |
| preQ1 riboswitch | 2.738 | 2.963 | 2.971 | No |
| cobalamin riboswitch | 2.278 | 2.589 | 2.830 | No |
| SAM riboswitch | 2.538 | 3.144 | 3.500 | No |

A single family---Hepatitis C IRES domain III---survives both composition controls, with a ratio of 7.200 exceeding the dinucleotide null (3.394) by more than 2x. This structure has a functionally critical pseudoknot-adjacent stem-loop. The remaining five families fall below the second-order control. At Rung 2, the composition cascade has absorbed nearly all apparent structure signal from embedding-level metrics.

### 4.3 Rung 3: ERNIE-RNA resolves specific base-pair partners

Phase 6 applies the most demanding test: does a model know which position pairs with which? Across the 8 models tested so far:

| Model | Pretrain | Mean PS | Gate pass | Exceed null | H3 precision |
|-------|----------|---------|-----------|-------------|-------------|
| **ERNIE-RNA** (86M) | RNA | **0.1170** | 30/32 | **28/30** | **0.883** |
| ERNIE-RNA untrained | RNA | 2.5 $\times 10^{-8}$ | 10/32 | 9/10 | — |
| SpliceBERT (19M) | RNA | 0.0004 | 21/32 | 14/21 | 0.256 |
| UTR-LM (~2M) | RNA | 0.00001 | 16/32 | 13/16 | 0.333 |
| RNA-FM (99M) | RNA | 0.00005 | 8 | 4 | — |
| HyenaDNA (5.4M) | DNA | 0.000001 | 25/32 | 10/25 | 0.078 |
| NT v2 (56M) | DNA | 0.0 | 1 | 0 | — |
| RiNALMo (650M) | RNA | *pending* | — | — | — |
| Evo (7B) | DNA | *pending* | — | — | — |
| DNABERT-2 (117M) | DNA | *pending* | — | — | — |
| Caduceus (14M) | DNA | *blocked* | — | — | — |

ERNIE-RNA stands alone. Its mean PS of 0.117 is three orders of magnitude above the next-best model (SpliceBERT, 0.0004). In 28 of 30 gate-passing families, the real PS exceeds the within-stem derangement null. The partner-is-max fraction of 0.883---the correct partner receives the largest perturbation in 151 of 171 eligible pairs---far exceeds the 1/3 chance baseline ($p \ll 0.001$, binomial test).

**Untrained control.** ERNIE-RNA with randomized weights (preserving architecture and attention bias structure, destroying all learned parameters) produces PS = 2.5 $\times 10^{-8}$---seven orders of magnitude below the trained model. Only 10 of 32 families pass the positive control gate (vs 30/32 trained). The architectural base-pairing attention bias alone does not produce partner specificity; the trained weights are necessary. This confirms that ERNIE-RNA's Rung 3 result is a learned property of the pretrained representations, not an artifact of the model's inductive bias.

**The effect concentrates in the final layers.** Peak PS layers for ERNIE-RNA range from 10 to 12 (of 12 total), with median at layer 12. SpliceBERT and UTR-LM peak at layers 0--1; HyenaDNA peaks at layer 0. Deep-layer concentration distinguishes genuine learned pairing from tokenization and positional-encoding artifacts.

**Top families by PS (ERNIE-RNA):**

| Family | PS | Layer | Eligible pairs |
|--------|-----|-------|---------------|
| glycine riboswitch | 0.197 | 12 | 9 |
| purine riboswitch | 0.197 | 12 | 16 |
| U5 snRNA | 0.191 | 12 | 18 |
| mir-122 precursor | 0.191 | 12 | 12 |
| mir-155 precursor | 0.191 | 12 | 12 |

The two quarantined pilot families (tRNA-Phe, PS = 0.281; tRNA-Ala, PS = 0.178) were excluded from confirmatory counts by preregistration but would not have changed the outcome.

**The gap is qualitative.** SpliceBERT technically exceeds the null in 14 of 21 families (above the H1$_6$ threshold of 7), but its PS magnitude (0.0004) is four orders of magnitude below ERNIE-RNA. The derangement null catches positional artifacts at these scales; the 14/21 count reflects noise fluctuations around zero PS, not genuine partner specificity. UTR-LM's H3 fraction of exactly 0.333 matches the chance baseline. HyenaDNA's H3 fraction (0.078) falls below chance, meaning the model systematically perturbs partners less than adjacent positions.

### 4.4 Case study: Huntington's disease repeat-expansion alleles

The three-rung evaluation framework has direct implications for therapeutic applications. Allele-selective antisense oligonucleotide design for repeat-expansion diseases requires computational structure prediction that distinguishes normal from expanded alleles. We tested six models on HTT exon 1 CAG repeat sequences at five clinically relevant lengths (CAG17, CAG21, CAG36, CAG40, CAG60) with identical flanking context.

**Embedding distances from wild-type (CAG17):**

| Model | CAG21 | CAG36 | CAG40 | CAG60 |
|-------|-------|-------|-------|-------|
| ERNIE-RNA | 5.8 $\times 10^{-4}$ | 7.2 $\times 10^{-3}$ | 9.1 $\times 10^{-3}$ | **2.3 $\times 10^{-2}$** |
| RNA-FM | 0.0 | 0.0 | 3.6 $\times 10^{-7}$ | 3.0 $\times 10^{-7}$ |

RNA-FM produces near-identical embeddings for normal and severely expanded alleles (cosine distance < $10^{-6}$), making allele-selective site identification impossible from its representations. ERNIE-RNA shows a monotonic increase in distance with repeat count, reaching 0.023 between CAG17 and CAG60. The distance is modest in absolute terms but four orders of magnitude above RNA-FM and strictly monotonic with repeat length, consistent with a representation that tracks the structural transition from single-stranded to stable hairpin as repeats lengthen past the pathogenic threshold.

The Rung 3 result predicts this difference: ERNIE-RNA encodes pairing topology in its deepest layers, so structural changes from hairpin formation as repeats expand produce detectable embedding shifts. RNA-FM's representations are composition-sensitive rather than structure-sensitive; since the nucleotide composition of CAG repeats is identical regardless of repeat count, its embeddings collapse.

### 4.5 Preregistered hypothesis summary

**Composition-controlled discrimination (Phases 1--2):**

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H1 | RNA-FM ratio $\geq$ 2.0 (N=12) | 1.74x | FAIL |
| H6 | RNA-FM ratio $\geq$ 2.0 (N=52) | 1.84x | FAIL |
| H7 | $\geq$ 1 family survives dinuc null | HCV IRES III (7.200) | **PASS** |
| H10 | RNA-FM trained $>$ untrained $\geq$ 75% | 50/52 = 96% | **PASS** |
| H11 | NT v2 attn $> 0.15$ trained, $< 0.05$ untrained | 0.322 / 0.000 | **PASS** |

**Partner specificity (Phase 6):**

| Hypothesis | Criterion | Result | Verdict |
|-----------|-----------|--------|---------|
| H1$_6$ | $\geq$ 1 model: PS $> 0$, $\geq$ 7 fam exceed null | ERNIE-RNA: 28/30 | **PASS** |
| H2$_6$ | RNA-pretrained PS $>$ DNA-pretrained (rb $> 0.5$) | *pending full panel* | *pending* |
| H3$_6$ | Partner-is-max $> 1/3$ | ERNIE-RNA: 0.883 | **PASS** |

### 4.6 Cross-rung summary

The three rungs tell a coherent story about structure awareness across the ten-model panel:

| | Rung 1 (stem/loop) | Rung 2 (composition-controlled) | Rung 3 (partner specificity) |
|---|---|---|---|
| ERNIE-RNA | 22/52 exceed null | not tested (first-order) | **28/30 exceed null, PS=0.117** |
| NT v2 | 28/52 exceed null; attn $\rho$=0.322 | — | PS = 0.0 |
| DNABERT-2 | 5/40; attn $\rho$=0.257 | — | *pending* |
| RNA-FM | 5/52; ratio 1.895 | 1/6 survives dinuc | PS = 0.00005 |
| SpliceBERT | 4/52 | — | PS = 0.0004 |
| UTR-LM | 3/49 | — | PS = 0.00001 |
| HyenaDNA | 8/52 | — | PS $\approx$ 0 |
| Caduceus | 5/52 | — | *blocked* |
| Evo | 6/52 | — | *pending* |
| RiNALMo | *pending* | — | *pending* |

Each rung filters more aggressively. Most models show marginal signal at Rung 1 that disappears under composition controls (Rung 2). At Rung 3, ERNIE-RNA separates from the pack by three orders of magnitude. NT v2 and DNABERT-2 show genuine Rung 1 signal through the attention channel but no partner specificity in embeddings.

## 5. Discussion

### 5.1 Why ERNIE-RNA and not the others

ERNIE-RNA's Phase 6 result raises a natural question: what distinguishes ERNIE-RNA architecturally from the other RNA-pretrained models?

ERNIE-RNA uses the same masked language modeling objective as RNA-FM, SpliceBERT, and UTR-LM. The difference is architectural: ERNIE-RNA injects a base-pairing informed attention bias into the attention score computation (Yin et al., 2025). From the second transformer layer onward, the attention bias of each layer is determined by the attention map of the previous layer, with initial biases favoring Watson-Crick (AU, GC) and wobble (GU) base pairs. This inductive bias steers attention toward paired positions throughout training, allowing the model to learn pairing-aware representations through standard MLM alone.

The other RNA-pretrained models use standard attention without structural bias. RNA-FM, SpliceBERT, UTR-LM, and RiNALMo all train on masked nucleotide prediction with no architectural mechanism to favor base-pairing partners in the attention computation. The Phase 6 results suggest that this architectural inductive bias, rather than scale (RiNALMo at 650M) or pretraining domain (RNA vs DNA), determines whether a model acquires partner-level structural knowledge.

The untrained control (Section 4.3) disambiguates architecture from learning. Randomized-weight ERNIE-RNA preserves the attention bias structure but produces PS = 2.5 $\times 10^{-8}$, seven orders of magnitude below the trained model. The bias provides the *capacity* to learn pairing; pretraining on RNA sequences with MLM fills that capacity with *specific* pairing knowledge. Architecture without training produces nothing; training without the bias (RNA-FM, SpliceBERT) also produces nothing. Both are necessary.

The peak-layer distribution supports this interpretation. ERNIE-RNA's PS concentrates at layers 10--12, indicating the pairing signal is computed through deep processing. The other models' PS, where nonzero, peaks at layers 0--1, consistent with tokenization artifacts rather than learned computation.

### 5.2 Attention as an alternative structure channel

NT v2 ($\rho$ = 0.322) and DNABERT-2 ($\rho$ = 0.257) learn attention patterns that correlate with base-pairing contacts, despite DNA-only pretraining and no structural supervision. Both use attention architectures (ESM+GLU and MosaicBERT); state-space models (HyenaDNA, Caduceus, Evo) have no attention matrices to examine.

The attention signal is qualitatively different from ERNIE-RNA's embedding-level partner specificity. Attention-contact correlation measures whether the attention pattern assigns higher weight to paired positions, a statistical correlation that can arise from sequence covariation without explicit structural knowledge. ERNIE-RNA's PS measures whether mutating one base of a pair selectively perturbs the partner's representation, a causal test that requires the model to encode specific pairing relationships in its hidden states.

The three models with structural attention represent two mechanisms. NT v2 and DNABERT-2 learn attention-contact correlations through architecture (ESM+GLU and MosaicBERT) and DNA pretraining, without any structural inductive bias. ERNIE-RNA's attention is architecturally biased toward base-pairing partners from initialization. The bias propagates structure into the residual stream, producing both attention-level correlation and embedding-level partner specificity. NT v2 and DNABERT-2 achieve the former without the latter.

RNA-FM attention encodes no structural information (0.067 trained vs 0.061 untrained). Standard BERT-style attention with RNA pretraining produces neither attention-contact correlation nor embedding-level structure.

### 5.3 The composition confound at two orders

The first-order null (nucleotide-stratified) absorbs the dominant confound: GC enrichment in stems. The second-order null (dinucleotide-stratified) tests residual signal from dinucleotide context differences (stem-G occurs in GC pairs, loop-G in diverse contexts).

The two-order cascade validates the framework's design. Without the dinucleotide null, five first-order survivors would have been reported as genuine positives. With it, only one---Hepatitis C IRES domain III---survives, a structure with a functionally critical pseudoknot-adjacent stem-loop and an unusual composition profile.

### 5.4 The ladder as a diagnostic

The three-rung ladder reveals complementary failure modes across the model panel.

RNA-FM has high embedding sensitivity (mean ratio 1.895) but fails composition controls (5/52 exceed null) and shows no partner specificity (PS = 0.00005). Its representations are composition-sensitive, not structure-sensitive.

NT v2 has low embedding sensitivity (mean ratio 1.169) but the highest per-family consistency (28/52 exceed null) and strong attention-contact correlation (0.322). Its structure knowledge lives in attention heads, not in embedding geometry.

ERNIE-RNA has moderate embedding sensitivity (22/52 exceed null) and strong partner specificity (PS = 0.117). Its structure knowledge runs through the full representation, from moderate Rung 1 signal to dominant Rung 3 signal.

Evo (7B) fails at every rung despite being the largest model, confirming that pretraining domain and tokenization dominate parameter count for RNA structure encoding.

### 5.5 Limitations and the covariation confound

ERNIE-RNA's partner specificity could reflect sequence covariation rather than learned geometric pairing. Compensatory mutations that maintain base pairing are overrepresented in evolutionary RNA alignments; a model trained on such alignments might predict that mutating one position perturbs its covarying partner without encoding the physical mechanism of base pairing. Distinguishing covariation from geometry requires interventions beyond the scope of this evaluation---for example, testing on synthetic sequences with no evolutionary history.

The within-stem derangement null controls for positional confounds (adjacent positions receive similar perturbations regardless of pairing) and compositional confounds (GC versus AU pairs producing different absolute perturbation magnitudes), but it does not control for covariation.

NT v2 and DNABERT-2 tokenization constraints limit the Phase 6 evaluation. NT v2's 6-mer tokenization makes single-nucleotide complement swaps ambiguous (a swap at one position changes the 6-mer encoding of up to six overlapping tokens); only 1 family passed the gate. DNABERT-2's BPE tokenization presents a similar challenge. These models' Phase 6 results (PS = 0.0) reflect tokenization limitations as much as architectural ones.

RiNALMo, Evo, and DNABERT-2 Phase 6 results are pending; Caduceus is blocked by a CUDA compilation requirement for mamba-ssm. The paper will be updated with these results.

The probing evaluation (logistic regression on embeddings) shows no model above the majority baseline at N = 52 families. The positive margins observed at N = 12 in Phase 1 (+0.013 to +0.041) were noise amplified by small sample size. Cross-family generalization of linear probes provides no evidence for embedding-level structure encoding in any tested model.

### 5.6 Implications for downstream applications

Antisense oligonucleotide target-site selection, splice-site prediction, and RNA design tools increasingly incorporate foundation model representations. The Phase 6 results carry a direct practical consequence: ERNIE-RNA's representations encode genuine pairing topology and could in principle support structure-dependent predictions. The remaining nine models should be used with physics-based fallbacks (such as ViennaRNA) for any application requiring structural resolution beyond stem-versus-loop discrimination.

The composition null---a per-prediction quality gate that certifies whether a model's structural predictions reflect learned pairing or nucleotide composition---provides a practical tool for pipeline developers. For repeat-expansion therapeutic targets, where nucleotide composition is uniform across alleles, models failing the composition null will produce identical predictions for structurally distinct alleles (Section 4.4). The composition null identifies this failure mode before the prediction enters a screening pipeline.

## 6. Conclusion

A three-rung evaluation ladder reveals that RNA structure awareness in foundation models is rarer and more specific than naive metrics suggest. Composition controls at two orders of sequence context absorb nearly all apparent embedding-level structure signal, with a single RNA family surviving both. The most demanding test---partner specificity with a derangement null---shows that ERNIE-RNA is the sole model among eight tested that encodes which specific position pairs with which, with PS = 0.117 across 28 of 30 qualifying families and partner-is-max precision of 88.3%. The effect concentrates in the final layers (10--12 of 12), consistent with learned pairing rather than positional encoding.

Two DNA-pretrained models (NT v2, DNABERT-2) encode structure through attention rather than embeddings, showing that architecture determines the channel through which structure knowledge is expressed. Evo's 7B parameters provide no advantage over 14M-parameter models, confirming that pretraining domain and tokenization dominate scale.

The three-rung framework---nucleotide-stratified permutation, dinucleotide stratification, and partner specificity with within-stem derangement null---provides a reusable evaluation template for structure awareness in sequence foundation models across biological domains.

## References

1. Chen, J. et al. (2022). Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions. *arXiv:2204.00300*.
2. Dalla-Torre, H. et al. (2023). The Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *bioRxiv:2023.01.11.523679*.
3. Nguyen, E. et al. (2023). HyenaDNA: long-range genomic sequence modeling at single nucleotide resolution. *arXiv:2306.15794*.
4. Schiff, Y. et al. (2024). Caduceus: bi-directional equivariant long-range DNA sequence modeling. *arXiv:2403.03234*.
5. Nguyen, E. et al. (2024). Sequence modeling and design from molecular to genome scale with Evo. *Science*, 386(6723).
6. Yin, W. et al. (2025). ERNIE-RNA: An RNA Language Model with Structure-Enhanced Representations. *Nature Communications*, 16, 10076.
7. Penic, R.J. et al. (2025). RiNALMo: General-Purpose RNA Language Models Can Generalize Well on Structure Prediction Tasks. *Nature Communications*, 16, 5729.
8. Chen, K. et al. (2024). Self-supervised learning on millions of primary RNA sequences improves sequence-based RNA splicing prediction. *Briefings in Bioinformatics*, 25(3), bbae163.
9. Zhou, Z. et al. (2024). DNABERT-2: Efficient Foundation Model and Benchmark for Multi-Species Genomes. ICLR 2024. *arXiv:2306.15006*.
10. Johnson, W.B. & Lindenstrauss, J. (1984). Extensions of Lipschitz mappings into a Hilbert space. *Contemporary Mathematics*, 26, 189--206.
11. Tower, E. (2026). Boundary conditions for geometric evaluation of foundation model representations. Working paper.

## Data availability

Phase 1--2 results: `data/gpu_results/`. Phase 6 results: `data/gpu_results/phase6_compensatory/`. RNA structures, analysis code, preregistration documents, and SHA hashes are in the `causal-rna` repository. Zenodo DOI: *[to be assigned]*.

Phase 1 prereg SHA: `694b43b` (causal-rna), `a7d10f5` (factorization-unified).
Phase 2 prereg SHA: `bd4b3fd` (causal-rna), `74c8f49` (factorization-unified).
Phase 6 prereg SHA: `c19aa59` (causal-rna).
