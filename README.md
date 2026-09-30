# Cyclotomic Weight Chambers for RC-Invariant Minimum Decycling Sets

**Oleksiy Babanskyy — version 3.2.0**

An exact comparison between a specified weight construction and the wider space of reverse-complement-invariant minimum decycling sets (RC MDS).

For DNA order eight, fix the labelled weights

\[
\mathcal W=\{(s(A),s(C),s(G),s(T))=(-\rho,-1,1,\rho):\rho>1\}
\]

and the declared lexicographic tie rules with `A<C<G<T`. This family meets **21 ordinary complete-move components**. The minimum residual length over their entire union is **141**, including endpoints that are not RC invariant. A certified **8230-word RC MDS has length 107** in another component, whose optimum is exactly 107.

The 21-component classification and the 141 union optimum concern **this fixed labelled family**. They do not classify all injective antisymmetric assignments of DNA weights. Changing the signs/order of the labelled weights while keeping ties fixed can change the component; other weight assignments are a separate open domain.

| Result | Meaning | Limit |
| --- | --- | --- |
| 21 components; union optimum 141 | Complete atlas and component certificates for the fixed family W | Not all possible labelled weight maps |
| 8230 words, residual length 107 | RC MDS outside that union; optimum of its own component | Not a global optimality theorem |
| 8232 words, residual length 105 | RC feedback-set construction | Minimum cardinality at length 106 remains unknown |
| Finite-pool closure | Exact signed-real characterization of a necessary walk test | Does not produce a Boolean endpoint |

Lengths count residual arcs in the lower de Bruijn graph. The global frontier remains **8230 ≤ tau_RC(106) ≤ 8232**. The paper does not settle the global optimum at MDS cardinality.

## Three reading paths

1. **Understand the result.** Start with the [article](paper/Cyclotomic_Weight_Chambers_for_RC_Invariant_MDS.pdf), main results and the explicit family W, then the certified separation section. The [reading guide](docs/READING_GUIDE.md) maps the eight main sections and technical appendices.
2. **Check 141 and 107.** Follow the [thematic certificate map](docs/ARTIFACT_MAP.md), then run `python -X utf8 -B scripts/verify.py --profile core`. The extended suite checks the complete fixed-family atlas and component lower bounds; the advances suite checks the literal 107 witness and its component certificate.
3. **Study completion and its obstructions.** Read the sections on stable completion and balanced walks, then the [balanced-load capsule](evidence/balanced_loads/README.md). One core passes an entire specified 1804-walk test but forty additional returns exclude its real-gradient class at cost below 8235. That lower bound is scoped, not global or asserted attained.

The arithmetic strand explains why only two form families determine the selector, the rational Farey spectrum and eventual cyclotomic realization. General proofs remain in the article and its technical appendices; [THEOREM_MAP.json](docs/THEOREM_MAP.json) maps formal labels to hypotheses and evidence.

## Verify

    python -m pip install -r requirements.txt
    python -X utf8 -B scripts/verify.py --profile core
    python -X utf8 -B scripts/verify_extended.py --optimized
    python -X utf8 -B scripts/verify_advances.py --optimized

Scientific replay uses exact arithmetic and local data, without a discovery solver or a private repository. The core profile includes all 22 advances groups and balanced-load replay in both Python modes. The balanced-load capsule also runs with the standard library alone. See [REPRODUCE.md](REPRODUCE.md) and [INDEPENDENT_VERIFICATION.md](INDEPENDENT_VERIFICATION.md) for the supported entry points and their scope.

The PDF identifies this companion by repository, version and an immutable digest of the delivered evidence tree. [results/evidence_identity.json](results/evidence_identity.json), [MANIFEST_SHA256.txt](MANIFEST_SHA256.txt) and [RELEASE_ATTESTATION.json](RELEASE_ATTESTATION.json) identify the evidence and complete distribution; these records do not imply an online release.

See the [claim boundaries](docs/PUBLIC_CLAIM_BOUNDARY.md), [scientific open problems](docs/OPEN_PROBLEMS.md), [contribution map](docs/ADVANCES.md) and [source lock](docs/SOURCE_LOCK.md). The imported parity interface is fixed at v2.0.0, commit `89a29a90962d02add6690d6411db927c4788371a`. Classical cyclotomic, retiming, SCC, 2-SAT, ordering and LP mechanisms are attributed. AI assistance and its limits are disclosed in the paper and [authorship statement](docs/AUTHORSHIP.md). Computational checks do not establish specialist acceptance or exhaustive novelty.

The MIT licence applies as stated in [LICENSE_SCOPE.md](LICENSE_SCOPE.md). Cite using [CITATION.cff](CITATION.cff).
