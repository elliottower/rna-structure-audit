> SUPERSEDED 2026-08-25. The RNA paper goes to the RNA journal only. XAI4Science
> was the closest fit of the NeurIPS workshops -- its scope names spurious
> correlation and probe faithfulness -- but the NeurIPS effort this year goes to
> Interpretability as a Science and TAI-Eval, both Sydney, both due 29 August, and
> both carrying different papers. The journal submission has no clock: it was never
> approved out of the proof queue. Kept for the venue survey and the framing notes,
> which still apply to the journal version.

# XAI4Science submission plan

Supersedes `ICBINB_BIO_PLAN.md`. The sections, numbers and displays carry over;
what changes is the framing, because this venue names the paper's subject in its
own scope statement rather than in its remit.

**XAI4Science: Knowledge Discovery and Trust through Interpretable Foundation
Models.** Sydney, 11--12 December 2026. Deadline 29 August 2026, 11:59pm AoE.
8 pages excluding references and appendices, NeurIPS 2026 template, OpenReview,
non-archival with concurrent submission explicitly allowed.

## Why this venue, in their words

> "XAI methods can help verify that FMs rely on physically meaningful patterns
> rather than spurious correlations, detect and diagnose physical
> hallucinations."

Stems are GC-rich and loops are AU-rich, so a composition-sensitive embedding
inherits apparent structure awareness. Separating the two is the paper.

Their second topic area is "attribution, probing, and mechanistic
interpretability methods, including rigorous evaluation of their faithfulness".
Section 3 is a faithfulness evaluation of the probe: an oracle whose partner
specificity is present by construction reads 1.000, and a local model with none
reads 0.155 against its own chance rate of 0.173.

Both are stated in their scope. Nothing has to be argued into fit.

## Title

Unchanged: **A Graded Evaluation of RNA Structure Awareness in Foundation
Models: From Stem--Loop Discrimination to Partner Specificity**

## Framing shift from the ICBINB plan

| | ICBINB framing | XAI4Science framing |
|---|---|---|
| lead | a control climbs every rung, and three more silent failures | a confound, a probe validated against it, and what the probe finds |
| the buffer | the headline | a major result, arriving after the instrument is established |
| the resolution floor | one of four failure modes | why one model is excluded, stated once |
| register | negative results | verification |

The negative findings do not go away. They stop being the organizing principle.

## Abstract, content order

1. Stems are GC-rich and loops are AU-rich, so apparent structure awareness may
   be composition. Separating them needs a probe whose faithfulness is
   established rather than assumed.
2. A three-rung ladder over ten RNA- and DNA-pretrained models on 47 Rfam
   families, each rung harder to pass by composition alone.
3. The probe reads 1.000 on a system with partner specificity by construction
   and its own chance rate on one without.
4. Two models resolve which position pairs with which; eight sit at their own
   chance rates. Pretraining domain does not predict which, and a 14M DNA model
   outscores three of the five RNA models.
5. A model with no learned weights climbs every rung, because a buffer holding a
   Watson-Crick table survives weight randomization. An ablation prices it: the
   random-init arm falls from +0.249 excess to +0.012 while the trained arm
   keeps +0.725.
6. What that implies for probes generally: a chance rate is a property of a
   model, not a constant, and a random-init control is not a null.

## Sections and page budget

| section | pages | content |
|---|---|---|
| 1 Introduction | 1.0 | The composition confound. Why partner specificity is the strong test. |
| 2 The probe | 1.0 | PS, the derangement null, the ladder. |
| 3 Is the probe faithful | 1.0 | Oracle 1.000, local 0.155, each against its own chance rate. Their topic 2. |
| 4 What the ladder finds | 1.75 | Rungs 1--3, the two models that pass, domain not predicting. |
| 5 A control that climbs every rung | 1.75 | 10/47, 6, 0.490. The buffer. The four-arm ablation. |
| 6 What this implies for probes | 0.75 | Chance is model-specific; random-init is not a null; measure the floor. |
| 7 Limitations | 0.5 | Derangement circularity; the unexplained sub-1/3 chance rate; two unattributed excesses; no rRNA. |

## Displays

1. Baseline hierarchy: trained / random-init / ablated / derangement, ERNIE-RNA.
2. Probe validation: oracle, local, and the ten models against their own chance.
3. Registered versus corrected verdicts, both columns.
4. Resolution floor per model, with mean PS overlaid.

## Open

- Anonymized repository link; double-blind, and linked material must preserve it.
- NeurIPS 2026 template port from the journal class.
- Whether Section 6 earns 0.75 pages or folds into the discussion.
