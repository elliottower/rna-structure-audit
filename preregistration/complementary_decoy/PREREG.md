# Does a mutation move its pairing partner more than a complementary decoy?

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

## Research questions or hypotheses

Rung 3 compares the partner $j$ of a mutated position $i$ with $j$'s two stem neighbors. A
representation that responds wherever the nucleotide of $j$ occurs wins that comparison
whenever the neighbors carry other nucleotides, without encoding which occurrence is the
partner. The earlier registration (`PREREGISTRATION_PHASE6_V2.md`) describes its derangement
as destroying only the partner assignment; the derangement does not keep the reassigned
partner complementary to the mutated base. This registration repairs that limitation.

A *matched decoy* for the pair $(i, j)$ is the partner-side position $j'$ of another
deep-interior pair of the same family with the nucleotide of $j$, excluding $j \pm 1$
(`scripts/decoy_comparison.py`, `matched_decoys`). A pair's score is the mean over its decoys
of 1 where $\Delta_j > \Delta_{j'}$, one half for a tie and 0 otherwise; the family *decoy
rate* is the mean of its pair scores. A response shared by every position carrying the
partner's nucleotide gives 0.5.

**H1.** For RiNALMo, across families whose decoy rate is not exactly 0.5, rates occur above
0.5 more often than below it.

**H2.** For ERNIE-RNA, the same.

**G1, G2.** The distance guard, for RiNALMo and ERNIE-RNA: the same holds when decoys are
restricted to positions no farther from $i$ than $j$ is,
$|i - j'| \le |i - j|$. G1 is read as a confirmation only if H1 holds, and G2 only if H2
holds.

The validation gate carries the design: without it H1 and H2 would hold for any
representation that passes Rung 3. H1 and H2 holding strengthen the paper's central result.
H1 or H2 failing means that this analysis does not establish that the model's Rung 3 result
exceeds the equal-response nucleotide-identity explanation; with 16 or 17 families a failure
can also be a lack of power. No hypothesis is registered for Caduceus or the remaining
models. On the stored results each of them other than RNA-FM keeps too few families for a
decision (Sample size); RNA-FM's population under the corrected adapter is not yet known.
All are analyzed descriptively.

## Foreknowledge of data or evidence

Seen: every Rung 3 result of the paper for all ten models, including per-pair precision, the
layer-selection audit and which families pass each model's stem-versus-loop positive control;
H1 and H2 were chosen after seeing them. The coverage of the comparison
(`results/audit/matched_decoy_counts.json`), from sequences and structures only. The
validation gate's output (`results/audit/decoy_control_validation_v2.json`): passed. One
post hoc diagnostic on stored results (`results/audit/precision_vs_complementarity.json`):
at each model's panel layer, RiNALMo's Rung 3 per-pair precision is 0.881 and ERNIE-RNA's
0.887 against 0.467 for the equal-response complementarity control on the same 35 families,
above it in 28 of 30 and 26 of 28 families that differ; the randomly initialized ERNIE-RNA is
at 0.479. (0.467 is the unweighted mean over families; pooled over the 190 deep-interior
pairs the control is at 87 of 190, 0.458.) That control moves every position of the partner's
nucleotide equally and is one nucleotide-only baseline, not a bound on all of them.

Not seen: the decoy rate of any model on any family. The delta profiles it is computed from
were never stored, so it cannot have been computed.

## Explanation of foreknowledge and managing unintended influences

The stored comparison makes H1 and H2 likely to hold under the identity criterion; it says
nothing about the distance guard. The decoy definition, the estimator, the layer rule and
the thresholds were fixed before any forward pass of the rerun, and the code that computes
them is sealed (`decoy_comparison.py`, hash `07dbec71…`). This is a targeted follow-up to the
Phase 6 confirmation and not part of it.

## Study type

Observational: fixed pretrained models, fixed sequences, a deterministic perturbation.

## Intention for causal interpretation

N/A — the perturbation is an intervention on the input and no causal claim about training is made.

## Blinding of experimental treatments

N/A — no treatment is administered by a person.

## Additional blinding during research or analysis

N/A — the statistic and its code are fixed before the run.

## Study design

Each eligible pair's first position is complement-swapped and cosine distance $\Delta_k$ between
wild-type and mutant embeddings is read at every other position and every layer
(`scripts/phase6_compensatory_mutation_v3.py`). Decoy rates, contributing pairs, comparisons
and ties are stored per family and layer for both scopes and both decoy sets.

## Randomization

N/A — the decoy comparison enumerates every decoy; nothing is sampled.

## Data collection procedures

The panel of `data/rfam_families` (sealed, tree hash `6a01623f…`), 36 confirmatory families;
the move from the 32 families of the earlier registration is the `DEVIATIONS.md` entry
"Unattributed annotations repaired; confirmatory set moves to N = 36". Models: RiNALMo,
ERNIE-RNA, RNA-FM (through `scripts/rnafm_adapter_v2.py`), UTR-LM, SpliceBERT, Caduceus,
HyenaDNA and Evo, and two control arms: the randomly initialized ERNIE-RNA with its
pairing-bias buffer and with the buffer ablated. NT v2 and DNABERT-2 are outside this analysis: their tokens span several
nucleotides, and DNABERT-2 does not rise above its no-op floor. Its result therefore does not
bear on a ranking of ten models.

The run, launched by `scripts/modal_decoy_panel.py` to a new volume:

| arm | models | stages |
|---|---|---|
| Phase 6 | RiNALMo, ERNIE-RNA, UTR-LM, SpliceBERT, Caduceus, HyenaDNA, Evo, NT v2, DNABERT-2 | Phase 6 only |
| Phase 6, randomly initialized | ERNIE-RNA, RiNALMo, SpliceBERT, UTR-LM, NT v2, DNABERT-2 | Phase 6 only |
| Phase 6, randomly initialized, pairing-bias buffer ablated | ERNIE-RNA (`--ablate-bias`) | Phase 6 only |
| RNA-FM | RNA-FM, randomly initialized RNA-FM | every stage |
| RNA-FM, other arms | transversion, synthetic sequences, multi-sequence replication, HTT fragments; randomly initialized with buffers reseeded | as in the deposited run |

The registered analysis reads the ten folders named above under Models. NT v2, DNABERT-2
and the other randomly initialized arms are rerun so that every row of the Rung 3 table
carries a chance rate measured on the pairs its precision is read on; they enter no
hypothesis here. Stages other than Phase 6 are rerun for RNA-FM alone, because no other
model's adapter or image changes.

## Data collection procedures - File upload

N/A — no file beyond the repository.

## Sample size

Counted from structure alone (`matched_decoy_counts.json`). Under family scope 17 families
have three or more pairs with a decoy: 116 pairs, 276 comparisons. 69% of those decoys are
farther from the mutation than the true partner; across all 276 comparisons, decoy distance
minus partner distance averages +34 nucleotides, which is why the distance guard is
registered. The no-farther set keeps 9 families, 38 pairs and 72 comparisons. Stem scope
keeps 5 families and is reported only as a description.

The population is narrower than the structure allows (Data inclusion and exclusion). On the
stored positive-control outcomes (`results/audit/decoy_population_counts.json`), RiNALMo keeps 16 of the 17 families and 9 of the 9;
ERNIE-RNA keeps 17 and 9. Caduceus keeps 6 and 4, and UTR-LM, SpliceBERT, HyenaDNA and Evo
between 5 and 7, so none of them can reach a decision.

With no family at exactly 0.5, an exact one-sided sign test at 0.025 needs 13 of 17 families
above 0.5 (13 of 16 as well), or 8 of 9; the requirement is recomputed after dropping families at
exactly 0.5. The identity criterion can resolve a clear majority of families; the distance
guard can resolve only a near-unanimous result and will miss a moderate one.

## Sample size rationale

Fixed by the panel; no family is added or removed for this analysis.

## Starting and stopping rules

One run per model; a run is complete when all 36 families are recorded.

## Manipulated variables

The nucleotide at the first position of each eligible pair, swapped to its Watson-Crick complement.

## Measured variables

$\Delta_k$ per position and layer; the decoy rate per family and layer; mean PS per family and layer.

## Measured variables - File upload

N/A — defined in code.

## Indices

The family decoy rate at a layer chosen without that family: the layer maximizing mean PS
over the model's other population families (leave one family out). Every population family
votes for the layer, whether or not it has a matched decoy, because having a decoy depends
on structure and not on outcome. PS must be finite at every layer for every population
family. The layer is chosen on PS and never on the decoy rate. A family enters a test when
three or more of its pairs contribute at that layer.

## Indices - File upload

N/A — defined in code.

## Statistical models

Exact one-sided sign test on the family decoy rates against 0.5, per model
(`scripts/decoy_audit_v2.py`). The test asks whether rates that are not exactly 0.5 fall above
0.5 more often than below it; it is not a test that the mean rate exceeds 0.5. Families at exactly 0.5 are counted, reported and dropped from
the denominator, $n_{\mathrm{sign}} = n_{>0.5} + n_{<0.5}$.

## Statistical models - File upload

N/A — defined in code.

## Transformations

None.

## Inference criteria

| decision | holds when |
|---|---|
| H1 | RiNALMo, decoy set "all": not void, and sign test $p < 0.025$ |
| H2 | ERNIE-RNA, decoy set "all": not void, and sign test $p < 0.025$ |
| G1 | H1 holds, and RiNALMo, decoy set "close": not void, and sign test $p < 0.025$ |
| G2 | H2 holds, and ERNIE-RNA, decoy set "close": not void, and sign test $p < 0.025$ |

H1 and H2 are each tested at 0.025, which is 0.05 divided over the two models. Each guard is
tested at 0.025 only behind its own hypothesis, so the two chains H1 then G1 and H2 then G2
keep the family-wise rate at 0.05. A guard is computed and reported whatever its hypothesis
does, and supports no confirmatory statement when that hypothesis fails.

**Validation gate.** `scripts/validate_decoy_controls_v2.py` must pass on the panel before
the analysis is read, and `decoy_audit_v2.py` refuses to run without its passing record. For every family, scope and decoy set, the number of contributing pairs
must equal the structure-only count. Where a matched pair exists the partner oracle must
score exactly 1 and the complementarity control exactly 0.5; a condition with no matched pair
is recorded as unavailable and not evaluated. **All decisions are void if the gate fails.**

**A decision is void with fewer than 14 families under "all" or 8 under "close"**, four
fifths of the 17 and 9 families the structure makes available.

Reported with every decision: families above, below and at 0.5; mean and median family rate;
pairs, comparisons and tie fraction; the layers chosen; the families excluded and why; and
the rate at the single panel layer, as a description.

## Data inclusion and exclusion

The population is that of the earlier registration: the 36 confirmatory families (the two
pilot families stay quarantined) that pass the model's stem-versus-loop positive control.
The same population chooses the layer. Which families pass is already known from the stored
results for every model except RNA-FM; the decoy rates are not. A comparison is dropped only
when a position's delta is absent from the profile at that layer.

## Missing data

A model with a failed run is rerun; no partial run is analyzed.

## Other planned analysis

Reported with no prediction and no decision: both decoy sets for Caduceus, RNA-FM, UTR-LM,
SpliceBERT, HyenaDNA and Evo, and for the randomly initialized ERNIE-RNA with and without its
pairing-bias buffer; stem scope for every model. For each of these the number of families at
exactly 0.5 is reported. The Table 3 columns restricted to the families every model passes
the gate on are reported as a description.

## Context and additional information

A registered prospective test of whether a partner-localized perturbation response, already
observed, exceeds one explicitly defined nucleotide-identity explanation. The same run
recomputes RNA-FM through a corrected adapter and stores chance rates on matched pairs at
every layer; those are corrections recorded in `DEVIATIONS.md`, not decisions of this
registration.

Sealed before this registration, by hash in `.results/ledger.jsonl`:

| file | hash |
|---|---|
| `scripts/stem_pairs.py` | `6abd3b7f…` |
| `scripts/decoy_comparison.py` | `07dbec71…` |
| `scripts/decoy_audit_v2.py` | `d6011ae9…` |
| `scripts/validate_decoy_controls_v2.py` | `2f5d6636…` |
| `scripts/count_matched_decoys.py` | `82feba72…` |
| `tests/test_decoy_comparison.py` | `b399b4d7…` |
| `tests/test_decoy_audit_v2.py` | `d952cfb5…` |
| `scripts/count_decoy_population.py` | `c0a8845a…` |
| `data/rfam_families` (tree) | `6a01623f…` |
| `results/audit/matched_decoy_counts.json` | `e1a53e04…` |
| `results/audit/decoy_control_validation_v2.json` | `14b45397…` |
| `results/audit/decoy_population_counts.json` | `90c45fde…` |
| `results/audit/precision_vs_complementarity.json` | `d17d0fca…` |
| `results/repaired_panel_v3/rinalmo/rinalmo_phase6_ps.json` | `36d9a36f…` |
| `results/repaired_panel_v3/ernierna/ernierna_phase6_ps.json` | `ebb01560…` |
| `results/repaired_panel_v3/ernierna_untrained/ernierna_untrained_phase6_ps.json` | `e6d20006…` |

`scripts/decoy_audit.py` and `tests/test_decoy_audit.py` are superseded by the v2 files above
and are not used. `scripts/phase6_compensatory_mutation_v3.py` and
`scripts/rnafm_adapter_v2.py` produce the delta profiles. Their tests need the model
libraries and run in the launch environment on a synthetic eight-nucleotide sequence. After
those tests pass and before the first forward pass on any study family or other analysis
input, the exact tested source files are sealed and their hashes written to this
registration's log with `prereg log`. The launcher records the commit and the full SHA-256
of every tested file beside the test report and refuses to run any model unless the files
it is about to run have those hashes. The statistic computed from the profiles is fixed by
the hashes above whatever those two files contain.

The corrected RNA-FM run changes the adapter implementation and the software environment
together (torch 2.4.1 and transformers 4.44.2 to torch 2.6.0, transformers 5.14.1 and
multimolecule 0.2.0). Differences from the deposited RNA-FM results therefore cannot be
attributed to tokenization, architecture or library versions singly. The deposited RNA-FM
results are superseded because the former adapter did not implement the pinned checkpoint;
the corrected run is read on its own and not as an estimate of any one component's effect.

**Maximum claim under this registration.** If the gate passes and H1 or H2 holds: at a layer
chosen without the family, the named model's response to a mutation is preferentially
localized to the annotated partner relative to other deep-interior paired positions carrying
the same nucleotide, so the response is not one shared by every such position. If G1 or G2
also holds: the same is true against decoys no farther from the mutation in linear sequence,
so nucleotide identity combined with ordinary distance attenuation does not explain it
either. If a guard is void or fails, the claim stops at the first sentence and the guard's
outcome is stated beside it.

Outside what the paper can say confirmatorily under this registration: anything about
Caduceus or the other models beyond a described rate; that any model lacks pairing
information recoverable by another readout; that the response encodes secondary-structure
topology and not some other feature correlated with the partner; anything about
non-canonical pairs; and any ranking of the ten models.
