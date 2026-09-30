# Guide to verification

1. Read the paper's units, [claim boundary](docs/PUBLIC_CLAIM_BOUNDARY.md) and [open problems](docs/OPEN_PROBLEMS.md).
2. Run `python -X utf8 -B scripts/verify.py --profile core`, then the optimized commands in [REPRODUCE.md](REPRODUCE.md).
3. Inspect the positional necessity/sufficiency and active-midpoint proofs; the retiming Boolean step uses nonnegative integers summing to one on every PCR.
4. Check RC mate sums and tight-SCC rounding, then the parity paths forcing an odd value of twice an integer.
5. Check the current endpoint: 8230 words, rank and attaining path 107, one-hit 108-cycle. Distinguish component optimum 107 from the open global optimum.
6. In stable completion include empty reachability paths, test D before adding arcs, and use restoration rather than deletion literals. The excess-one extension may force both literals equal.
7. Inspect the bounded-chart necessity proof, its integral clock, the height rounding argument, and both types of surplus correction.
8. In the finite-pool theorem retain multiplicities, signed real lift coordinates, the aggregate budget b >= 0 and the length condition T >= k. Follow both directions of normalization and LP duality. Check that the 8235 lower applies only to the named nonstable real-gradient class, and distinguishes return lengths from closed-walk lengths.
9. Check each exclusion against its exact base, domain and operation order. In particular, radius before retiming differs from radius after changing the base by retiming.

The base-four endpoint checker and excess-one regression import no submitted mathematical checker. The literal-word modules are also separate from discovery code, but cumulative rounds share reductions and data. The old two-family regression reuses a sign oracle, so its independence is narrower. All work is AI assisted. The revised manuscript and companion have distinct AI-assisted readings in addition to preparation self-review. The balanced-load capsule reuses a source graph kernel. External specialist review is still needed.
