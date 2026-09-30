# Proof and computation boundary

The manuscript is the authority for the all-order results T01--T11. Their
proofs use no finite atlas row.

`verify_single_odd_prime_embedding.py` is a standard-library,
proof-independent finite assurance layer for T09. It checks the cyclotomic
composition identity for 18 `(p,a,b)` triples; reconstructs exact binary word
counts at `k=6,10,12,18,20,24`; verifies the `k=12` PCR count; and uses
`k=30` as a negative-domain control. The `k=30` row is deliberately outside
T09 because 30 has two distinct odd prime divisors. Passing these checks
guards proof typing and arithmetic but does not replace the all-order proof.

The executable capsule is the authority for the bounded proposition F01:

- all 3 point/interval cells at `k=4`;
- all 23 point/interval cells at `k=6`; and
- 4 named algebraic landmarks at `k=8`.

This is an exact finite theorem, not empirical sampling. It does not extend
its quantifiers beyond those thirty rows.
