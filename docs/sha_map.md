# Commit id map

`26` commits, rewritten from `rna-structure-awareness` into this repository with
git-filter-repo. Every id in this repository's history is the new one; documents
frozen before the rewrite quote the old one. Both are given at seven characters,
the length used throughout the manuscripts.

Documents quoting old ids are left unedited: a preregistration is a record of what
was frozen, and correcting its self-reference after the fact would make it a record
of something else. Resolve the id here instead.

| old | new | subject |
|---|---|---|
| `5427db9` | `c6a7e41` | Freeze batch-1 scorer and preregistration for SHA freeze |
| `5184151` | `66e6b25` | Fill SHA c6a7e41 into frozen preregistration |
| `3b61069` | `96c1aa4` | Batch 2 experiments B4-B7: EAP-IG, DAS, phase separation, splice-site DI |
| `6ec1c4e` | `f7d79a1` | Dense DAS sweep (n=71), behavioral DAS, CAA comparison, cross-disease replication |
| `9e6a2e4` | `9c0e16a` | Cross-gene direction cosine similarity: two-cluster structure, mean |cos|=0.73 |
| `0f1b478` | `2c0220d` | Add main_v4.tex: incorporate batch-2 causal probing results |
| `694b43b` | `a207535` | Preregistration: cross-architecture RNA structure metrics |
| `ae38128` | `69b9fb1` | Fill preregistration commit SHAs — freeze complete |
| `bd4b3fd` | `ae6712e` | Preregistration: expanded Rfam evaluation (H6-H11), Paper C drafts v1+v2 |
| `15fab48` | `1086452` | Annotate prereg with freeze SHAs: ae6712e (causal-rna), 74c8f49 (factorization-unified) |
| `01c5c84` | `322fdb6` | Paper C v3: Phase 2 results filled (H6 FAIL, H7 PASS, H10 PASS, H11 PASS, H8 FAIL; Evo pending) |
| `a956c2b` | `c6d611c` | Paper C v3: fill Evo results (1.352 mean ratio, 6/52 > null) — all 6 models complete |
| `672f269` | `4929ae8` | Phase 4: add ERNIE-RNA, SpliceBERT, DNABERT-2 results + paper v11 |
| `e9f2149` | `ef74dff` | Phase 6 preregistration: compensatory double mutation experiment |
| `9bdb2bb` | `a6d2f70` | Phase 6 script validated locally: coupling ratio metric works |
| `c19aa59` | `891d6af` | Phase 6 perturbation specificity: prereg V2, PS script, Modal wrapper, analysis |
| `6c46aba` | `d99e4df` | Freeze prereg Phase 6 V2 with SHA 891d6af |
| `e0dfd66` | `b904640` | Add rna-structure-audit pip-installable benchmark package |
| `b700a51` | `6f9b855` | Add experiment results, scripts, and data |
| `9c70d23` | `9503114` | Preregister Phase 6 untrained RiNALMo control |
| `5f82efc` | `01d4624` | Add Modal script for untrained RiNALMo Phase 6 control |
| `ec6f7fe` | `15a7bb2` | Add PyPI publish workflow via trusted publisher |
| `2772244` | `6c83904` | Add 10 pre-built adapters, notebooks, and Colab badges |
| `8d939d3` | `38b43ef` | Fix publish workflow: add contents:read permission for private repo checkout |
| `ce1b92b` | `fd43dde` | Update data availability: pip install rna-structure-audit now live on PyPI |
| `a7af0e6` | `c87c85b` | Rename submission files: rna-structure-audit.tex + cover-letter.tex |

## Ids this map does not cover

`0aff46c`, `36da8a1`, `74c8f49`, `a7d10f5`, `bbae163`

These are quoted in superseded drafts and preregistrations under `docs/` and are not
in the rewrite's commit map. Two of them (`0aff46c`, `36da8a1`) resolve in the
pre-rewrite bundle `rna-structure-audit-BACKUP.git`, which is kept outside this
repository; the rest predate that bundle or are not commit ids at all -- `bbae163`
is a *Briefings in Bioinformatics* article number. No manuscript in `paper/` cites
any of them: `paper_v11.tex` pins `a207535`, `ae6712e` and `891d6af`, all of which
resolve here.
