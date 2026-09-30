# Reading guide

The article follows one chain: fixed weights → selector and chambers → complete-move components → the certified 141 barrier → an RC MDS at 107 outside those components → exact completion tools and their limitations.

## Main text

| Section | Question | What to check |
| --- | --- | --- |
| 1. Problem and main results | What is proved, on which family? | W=(-rho,-1,1,rho), rho>1, fixed labels and ties; component versus global bounds |
| 2. Construction and chambers | How does the weight selector vary? | Imported construction, two-family decision rule, equality cases and finite chambers |
| 3. Retiming and RC sectors | How is an entire component optimized? | Integer gradients, cycle ratios, rounding and exactly two paired-move sectors |
| 4. Certified separation | Why are 141 and 107 exact within their stated domains? | Complete W-atlas, cycles proving lower bounds, endpoint witnesses proving upper bounds |
| 5. Endpoint completion | What remains after fixing nonstable choices? | Conditional 2-SAT, bounded integral charts and heights; no global polynomial-time claim |
| 6. Balanced walks | What can all combinations of a finite pool prove? | Aggregate budget b>=0, signed real lifts, and the scoped strictness example |
| 7. Arithmetic structure | Which resonances can the chambers have? | Main arithmetic interfaces and the technical proofs referenced in the supplement |
| 8. Limits and evidence | What is still open and where are the certificates? | Global cardinality gap, other labelled weight families, companion identity and replay scope |

The technical appendices retain the arithmetic proofs, historical explicit 131 construction, finite exclusion domains and assurance details. They are part of the complete PDF. The prime-sine independence argument appears once, with cross-references elsewhere.

## From a statement to its evidence

Use the stable LaTeX label in [THEOREM_MAP.json](THEOREM_MAP.json); numerical theorem labels may change between versions. The build receipt [paper_build.json](../results/paper_build.json) records the number and page in this exact PDF. The [artifact map](ARTIFACT_MAP.md) groups the companion by mathematical role rather than historical round number.

For the finite separation, inspect `prop:component-optima` and `thm:construction-107`. For the new obstruction mechanism, inspect `thm:finite-pool-closure` and `prop:core-pool-strictness`. Replay confirms their finite inputs; it does not replace reading the general proofs.
