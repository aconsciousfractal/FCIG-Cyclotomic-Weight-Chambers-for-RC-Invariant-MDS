# Comparison with the published DNA order-eight path-length bound

## Quantity and antecedent

Let D(4,8) be the directed de Bruijn graph whose vertices are the DNA words of
length eight, with one-symbol-shift edges and loops. A minimum decycling set M
has 8,230 vertices. Let ell(M) count the vertices of a longest directed path in
the residual graph D(4,8) minus M. The optimization problem is to minimize
ell(M) over all such minimum decycling sets.

Marçais, DeBlasio and Kingsford study this optimization in Section 5.2 of
*Sketching Methods with Small Window Guarantee Using Minimum Decycling Sets*,
[Journal of Computational Biology 31(7), 597–615 (2024),
DOI 10.1089/cmb.2024.0544](https://journals.sagepub.com/doi/10.1089/cmb.2024.0544).
The [initial preprint](https://arxiv.org/html/2311.03592v1) was posted on
6 November 2023. Table 2 gives 145 for both the Mykkeltveit construction and
the SA Min search at alphabet size four and order eight. For this parameter
pair, SA Min is an achieved search value, not a proved minimum.

## The certified improvement

The construction in Section 4 of the present manuscript provides an MDS of
8,230 words with ell(M)=107. It also satisfies reverse-complement invariance.
Thus the ordinary minimum-cardinality optimization has upper bound 107,
improving the cited achieved value 145. The RC property is additional: it is
not an extra premise silently inserted into the earlier ordinary problem.

The certificate proves 107 is optimal in the witness's own ordinary
complete-move component. It does not prove a global optimum over all MDSs.
The 141 bound elsewhere in the manuscript concerns only the union of the
21 components met by the fixed labelled weight family W.

## Counting convention

The manuscript counts residual arcs in the lower graph B(4,7). Its arcs are
the eight-letter words, so a residual lower-graph path of 107 arcs corresponds
to 107 vertices and 106 edges in the upper residual graph D(4,8). The vertex
count is the convention used in the cited Table 2, consistently with its
Proposition 7.

For a small exact calibration, all 30 binary order-four MDSs have longest
residual paths between five and seven vertices, or four and six edges. The
cited table reports five and seven. In DNA order eight the M107 witness has a
longest avoiding sequence of 114 letters. Every sequence of 115 letters, or
every window of 108 consecutive eight-mers, contains a selected word.

## Relation to symmetric cardinality and scope

The earlier parity paper addresses the minimum number of selected words under
reverse-complement symmetry, following Marçais, Elder and Kingsford,
[*k-nonical space: sketching with reverse complements*, Bioinformatics 40(11),
btae629 (2024)](https://doi.org/10.1093/bioinformatics/btae629), Section 6.
The present comparison concerns path length at already minimum cardinality.
In Table 3 of that source, the DNA order-eight value 89 belongs to an RC union
of 16,440 words and therefore does not improve a bound at cardinality 8,230.

This note makes a comparison with specific published values. It does not
certify an exhaustive literature-wide record, establish biological performance,
or resolve the general connectivity or asymptotic conjectures in the sources.
The global cardinality at lower-graph path bound 106 remains between 8,230
and 8,232. See [PUBLIC_CLAIM_BOUNDARY.md](PUBLIC_CLAIM_BOUNDARY.md).

## Reproduction

The literal witness is
[`M107.words.txt`](../evidence/advances/r5/witnesses/M107.words.txt), SHA-256
`4e0c44d87ea9782f261cf3ed2c0e918784fa404bbf0f57041b494509ff7d5edc`.
The [endpoint verifier](../evidence/advances/verify_endpoints.py) checks
cardinality, RC closure, necklace coverage, acyclicity, the exact path length
and the component certificate. The core profile includes the advances checks
and the fixed-family atlas. From the repository root, run:

```text
python -X utf8 -B scripts/verify.py --profile core
```

For a targeted invocation of the endpoint verifier, its required `--out`
argument can name a file in an existing temporary directory outside the
distribution. These are exact finite certificates, distinct from the general
proofs in the manuscript and from independent specialist review.
