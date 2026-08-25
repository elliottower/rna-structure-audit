# ICBINB-BIO submission plan

8 pages excluding references and appendices. Double-blind. Non-archival, so the
RNA journal version is unaffected. Deadline 29 August 2026, 11:59pm AoE.

## Title

Unchanged from the journal version:

**A Graded Evaluation of RNA Structure Awareness in Foundation Models:
From Stem--Loop Discrimination to Partner Specificity**

"RNA structure awareness" is awareness of RNA structure, which is what every
model in the panel is tested for whether it was pretrained on RNA or DNA, and
"in Foundation Models" is already domain-agnostic. Nothing needs adding.

## Strategy: one paper, three venues

ICBINB-BIO is non-archival and welcomes concurrent submissions, so this is not a
second paper. It is the same work in three forms.

| venue | date | form |
|---|---|---|
| ICBINB-BIO | 29 Aug 2026 | 8 pages, double-blind, non-archival |
| ICLR 2027 | abstract 18 Sep, paper 25 Sep 2026 | main track, full |
| RNA journal | open | full, no clock -- the submission was never approved |

The workshop lands four weeks before the ICLR deadline and does not count as
prior publication, so it is a dry run on the strongest version of the work with
outside reviews in hand before the main-track submission.

## The thesis

A three-rung ladder, each rung harder to pass by accident than the last, is
climbed at every rung by a model with no learned weights. ERNIE-RNA with random
weights exceeds its Rung 1 null in 10 of 47 families against a chance expectation
of 2.4, retains 6 at Rung 2, and reaches 0.490 per-pair precision at Rung 3.
The cause is a buffer holding a Watson-Crick table that weight randomization
cannot reach, and an ablation prices it: the random-init arm falls from +0.249
excess over its own chance rate to +0.012, while the trained arm keeps +0.725.

Two models resolve base-pair partners once that is controlled. Pretraining domain
does not predict which.

## Abstract, content order

1. The confound: stems are GC-rich, loops are AU-rich, so any
   composition-sensitive embedding inherits apparent structure awareness.
2. The design: a three-rung ladder, ten models, 47 Rfam families.
3. Rungs 1 and 2: a nucleotide-stratified null absorbs most apparent signal; a
   dinucleotide null separates two models from eight.
4. Rung 3: two models resolve which position pairs with which.
5. The finding that qualifies all of it: a model with no learned weights climbs
   every rung, because a buffer holding a Watson-Crick table survives weight
   randomization. An ablation prices it at +0.249 of excess, falling to +0.012.
6. What survives: two models resolve partners; pretraining domain does not
   predict which, and a 14M DNA model outscores three of the five RNA models.

## Sections and page budget

| section | pages | content |
|---|---|---|
| 1 Introduction | 1.0 | The composition confound. Why a graded ladder. What the paper finds. |
| 2 The ladder | 1.0 | Three rungs, the metric at each, the null at each. |
| 3 Does the probe work | 0.75 | Oracle 1.000, local 0.155, each against its own chance rate. The instrument before the claims. |
| 4 Results across the ladder | 1.75 | Rungs 1 and 2 tables, Rung 3 table, the two models that pass. |
| 5 A control that climbs every rung | 1.5 | Random-init ERNIE-RNA at 10/47, 6, and 0.490. The buffer. The ablation, four arms. |
| 6 What does not predict partner specificity | 1.0 | Domain, parameter count, architecture family. Caduceus above three RNA models. |
| 7 What the measurement can resolve | 0.5 | The no-op floor per model; one model excluded because its floor exceeds its signal. |
| 8 Limitations | 0.5 | Derangement-null circularity; the unexplained sub-1/3 chance rate; two unattributed random-init excesses; no rRNA in the panel. |

## Displays

1. **Baseline hierarchy** (table). Trained / random-init / ablated / derangement,
   with ERNIE-RNA's four numbers. The paper's centrepiece.
2. **Probe validation** (table). Oracle, local, and the ten models, precision
   against own chance.
3. **Registered vs corrected** (table). Per model, both verdicts, three moving
   and three of the four movements running against the paper's own conclusions.
4. **Resolution floor** (figure or small table). Per-model no-op floor on a log
   axis, with mean PS overlaid; DNABERT-2 above its own signal.

Four displays over roughly 5,000 words is about 0.8 per 1000, under the mechval
benchmark, so a fifth is affordable if Section 5 needs one.

## Decisions applied

- **D3** DNABERT-2 excluded from Rung 3, floor as the stated reason, and the
  exclusion reported as a finding.
- **D4** the control arm is "random-init (architecture preserved)" throughout.
- **D5** both registered and corrected verdicts reported, in one table.
- **D6** H17 is FAIL. The count and effect size sit adjacent; the arbitrariness
  of a 75% family-fraction threshold goes in Limitations.
- **D7** the RiNALMo and RNA-FM random-init excesses are reported as open, with
  the buffer inventory.
- **D9** the line: a defect is a result where it changes what a number means, and
  repository material otherwise. The buffer, the floor, the test size and the
  metric are results. The twenty-entry defect register is not.
- **D10** no ribosomal RNA in the panel, stated as a limitation.

## Cut from the journal version

- The synthetic covariation comparison. Those cells were never re-run on the
  repaired panel and are not part of this thesis.
- Rungs 1 and 2 in detail; one paragraph, since partner specificity carries the
  argument.
- The per-family tables. Appendix, if at all.

## Open before drafting

- Anonymized repository link, via an anonymizing mirror.
- Whether Section 5 earns 1.5 pages or should compress to 1.0 to give Section 4
  more room.
