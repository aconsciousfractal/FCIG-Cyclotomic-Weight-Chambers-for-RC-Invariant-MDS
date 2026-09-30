# Validation protocol

For each of the thirty declared rows, both implementations independently:

1. construct the complete occurrence-level decision-form catalog;
2. evaluate the exact sign covector;
3. build the seed, midpoint ideal, and direct terminal;
4. replay the legal firings and compare with the direct formula;
5. check Booleanity, literal reverse complement, and exactly one word per PCR;
6. construct every residual arc;
7. prove acyclicity by Kahn elimination and depth-first search;
8. compute a topological order, longest-path histogram, sources, sinks, and
   degree extrema; and
9. compare every terminal and midpoint ideal within the same order.

Before constructing any atlas row, the validator independently rebuilds the
complete `k=4,6` occurrence catalogs, canonical forms, critical roots,
supporting-form and occurrence multiplicities, ordered point/interval cells,
adjacency, representatives, landmark locations, and all associated hashes.
It therefore does not accept the critical-ratio order merely because the
input receipt is self-consistent.

The validator also checks both JSON schemas, exact canonical encoding, the
producer payload hash, complete equality with its independent reconstruction,
twelve native branch witnesses, and twenty-nine typed semantic mutations. The
last mutation replaces the first `k=6` critical ratio while coherently updating
the adjacent cells and ordered-root digest; the independently reconstructed
geometry rejects it.

The complete declared range is `k=4,6`. The `k=8` rows are four landmarks,
not a complete diagram.
