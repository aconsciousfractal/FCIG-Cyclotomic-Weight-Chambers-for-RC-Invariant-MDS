# Thematic artifact map

The complete article is in [paper/](../paper/main.tex). General proofs remain in the main text and technical appendices. Stable labels and evidence scope are indexed in [THEOREM_MAP.json](THEOREM_MAP.json).

| Mathematical role | Delivered evidence | Supported entry point |
| --- | --- | --- |
| Small chamber regression | [finite_assurance](../evidence/finite_assurance/README.md): 30 rows at the declared orders/landmarks | `scripts/verify.py --profile core` |
| Single-odd-prime formulas and domain checks | [theorem](../evidence/theorem/README.md) | `scripts/verify.py --profile core` |
| Arithmetic realization and selector extensions | [extended](../evidence/extended/README.md) | `scripts/verify_extended.py` |
| Complete fixed-W atlas and the 21-component union | [advances](../evidence/advances/README.md), the component panel plus `evidence/extended/atlas` and `evidence/extended/retiming` | `scripts/verify_advances.py` |
| The 8230-word M107 witness and component optimum | `evidence/advances/r5/witnesses/M107.words.txt`, `r5/certificates/M107.json` | `scripts/verify_advances.py` |
| Stable completion, charts, heights and scoped exclusions | `evidence/advances/r8` through `r13`, with dependency capsules named by their local manifests | `scripts/verify_advances.py` |
| Balanced finite-pool closure and class strictness | [balanced_loads](../evidence/balanced_loads/README.md): literal core, signed lift and five complete return families | `scripts/verify_balanced_loads.py` |

Paths abbreviated as `rN/...` in the paper are relative to `evidence/advances/`. Historic round names identify certificate provenance; they are not the article's reading order. Repeated payloads are retained where standalone capsules and their exact checks depend on them. Three unused historical report builders requiring undistributed documents are omitted; supported scientific replay is retained.

## Identities and boundaries

- [results/evidence_identity.json](../results/evidence_identity.json) identifies the delivered evidence tree and is quoted in the PDF.
- [results/paper_build.json](../results/paper_build.json) records two identical clean builds, source identity and exact statement pages.
- [MANIFEST_SHA256.txt](../MANIFEST_SHA256.txt) and [RELEASE_ATTESTATION.json](../RELEASE_ATTESTATION.json) bind the full distribution.
- [PUBLIC_CLAIM_BOUNDARY.md](PUBLIC_CLAIM_BOUNDARY.md) distinguishes W, one ordinary component, fixed-exterior fibres, real-gradient classes and the global problem.
- [OPEN_PROBLEMS.md](OPEN_PROBLEMS.md) contains scientific questions only.

The old k=30 negative-domain test remains valid: the later two-prime formula is a different theorem, not an extension of the old formula outside its hypotheses.
