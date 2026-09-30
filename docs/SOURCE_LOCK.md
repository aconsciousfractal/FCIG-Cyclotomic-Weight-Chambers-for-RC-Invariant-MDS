# Source lock and role map

The bibliography in `paper/refs.bib` is the authoritative citation record.
This document records why each source class is used; it does not redistribute
third-party full text.

| Role | Primary sources and locators | Imported content |
| --- | --- | --- |
| MDS foundations | Mykkeltveit (1972), DOI `10.1016/0095-8956(72)90006-8`; Champarnaud--Hansel--Perrin (2004), DOI `10.1142/S0218196704001700` | classical seed, unavoidable sets, cardinality context |
| weighted RC selector | Babanskyy (2026), *A Parity Theorem for Reverse-Complement-Invariant Minimum Decycling Sets in q-ary de Bruijn Graphs*, release `v2.0.0`, commit `89a29a90962d02add6690d6411db927c4788371a`, tree `7421eb0006306f1350a20208776cb8a9b3e1a337`, PDF SHA-256 `468b627fe9d3e4c05e723f4fdd9c8caf10abc88466cd0bee9d8a64beac7e0b6`; Theorems 4.3, 6.3, 7.2 and Corollary 7.3 | imported selector interface and parity theorem |
| genomic/sketching context | Zheng--Kingsford--Marçais (2021), DOI `10.1089/cmb.2020.0432`; Pellow et al. (2023), DOI `10.1101/gr.277644.123`; Marçais--DeBlasio--Kingsford (2024), DOI `10.1089/cmb.2024.0544`; Marçais--Elder--Kingsford (2024), DOI `10.1093/bioinformatics/btae629` | motivation; the Marçais--DeBlasio--Kingsford component/cycle interface is also used as detailed below |
| cyclotomic fields | Washington (1997), 2nd ed., Proposition 2.7 and standard conductor/ramification consequences, DOI `10.1007/978-1-4612-1934-7` | field inclusions and conductor facts |
| root-of-unity relations | Mann (1965); Conway--Jones (1976); Lam--Leung (2000); Steinberger (2008); Koley--Arikatla (2019); Christie--Dykema--Klep (2020); Kiss et al. (2026); De Carli--Laporta (2026) | abstract relation and divisibility context |
| exact algebraic sign | Emiris--Tsigaridas (2005), arXiv `cs/0512072`, Theorem 37 | sign of represented real algebraic expressions after dense variable-degree reduction, degree/height accounting, and sign-preserving padding stated in the proof |
| covector language | Björner et al. (1999), *Oriented Matroids*, 2nd ed., Section 1.2 | terminology only |

Every citation key used by the manuscript appears exactly once in
`paper/refs.bib`. Context sources are not represented as proving the new
weight-to-terminal or certificate statements.

The uniform radical family, twice-prime embedding counts, rational
finite-energy corollary, and typed equality-level Galois transport are proved
in the manuscript from the imported standard algebraic facts. No new
bibliographic source is represented as subsuming those selector-specific
statements. Only the declared additive instance and component optimization are classified; the general terminal-arrangement and global optimization programmes remain open.

## Version 2.0.0 additions

Martin--Reiner, Israel J. Math. 150 (2005), 229–240, arXiv math/0402206: squarefree tensor kernel context; the bounded positional realization step is explicit here.

Marçais--DeBlasio--Kingsford, arXiv 2311.03592v1 §§4–5 / DOI 10.1089/cmb.2024.0544: complete-move components and cycle invariants, now mathematical antecedents beyond the original motivation role.

Leiserson--Saxe, Algorithmica 6 (1991), 5–35, DOI 10.1007/BF01759032: retiming context, with a self-contained application proof. Lahiri--Musuvathi, MSR-TR-2005-67 §4.2: tight-SCC/UTVPI integer feasibility. RC/PCR specialization is explicit.

The old exclusion of all additive arrangements is superseded for the proved Freiman criterion and three-magnitude order-two instance. General terminal arrangements remain open.

## Version 3.0.0 additions

Aspvall--Plass--Tarjan (1979), Information Processing Letters 8(3), 121–123, DOI 10.1016/0020-0190(79)90002-4, with the 1982 erratum: classical implication-graph 2-SAT mechanism. Publisher metadata and scope were consulted; no claim of a new 2-SAT algorithm is made. The RC reachability reduction and length formula are proved in the manuscript.

Kenkre--Pandit--Purohit--Saket, arXiv:1507.00662v1, abstract and formulation: classical Max-k-Ordering antecedent. The bounded antisymmetric RC height equivalence is proved in the manuscript; the approximation theorems are not imported.

The finite inputs are frozen in evidence/advances/PROVENANCE.json. The two new consolidation corollaries combine stable completion with surplus classification and transport two-arm cuts correctly across cardinality regimes. This is not an exhaustive novelty search.

## Version 3.2.0 consolidation

The mechanism comparison is sharpened: Marçais–DeBlasio–Kingsford, arXiv 2311.03592v1, §§4.1–4.2 and Proposition 4 supplies the cycle-signature component antecedent. Leiserson–Saxe, Theorem 11 and Corollary 12, supplies the unit-delay retiming comparison. Lahiri–Musuvathi §4.2, Algorithm 1 supplies the tight-SCC integer-feasibility antecedent. Martin–Reiner §2 supplies the tensor/direct-sum kernel context; De Carli–Laporta, Theorem 1.1 and §3 concerns coefficient divisibility, rather than the bounded DNA positional alphabet. Kenkre et al. §2, Definition 1 and Lemma 1 supplies the ordering/bounded-DAG equivalence used as context.

Boyd–Vandenberghe, *Convex Optimization* (2004), Chapter 5, is the source for finite LP duality. The authors’ [book and lecture material](https://web.stanford.edu/~boyd/cvxbook/) state the feasibility/finite-optimum duality used here. The signed-lift theorem writes the primal and dual explicitly; it does not claim a new LP duality theorem.

The balanced-load capsule preserves seven imported payloads byte for byte, listed with hashes in evidence/balanced_loads/PROVENANCE.json. Its verifier uses a shared exact graph kernel. Independent implementation is not asserted for that capsule. The paper uses tail-minus-head gradients; imported head-minus-tail data are translated with the opposite sign.

## Fixed-family scope of the finite atlas

The imported selector permits general injective antisymmetric labelled weights. The delivered complete DNA8 atlas specializes to W=(-rho,-1,1,rho), rho>1, in labelled order (A,C,G,T), with the fixed A<C<G<T tie rules. The 21-component classification and 141 union optimum refer only to this specialization. Relabelling weights while leaving the tie rules unchanged is not an equivalence used by the article.
