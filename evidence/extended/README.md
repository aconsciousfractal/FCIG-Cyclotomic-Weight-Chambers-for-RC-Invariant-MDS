# Extended exact evidence

All scientific inputs are included. No inherited archive, network, solver service or private repository is required.

| Directory | Coverage |
|---|---|
| algebra | Rational spectra, witnesses, two-prime counts, midpoint classification, additive arrangement, two-family and max-plus regression |
| atlas | Complete k=4,6,8 atlas: 199 strata; 404 regenerated products; 11,432,704 membership decisions; 18 mutations |
| walls | Mixed-wall identities; 3857 paired moves; 213 separators; 15 mutations; bounded probes |
| retiming | 72 MDS certificate occurrences; 42 sector solutions; 30 lower cycles; 31 separators; four patches to 131; 24 mutations and small exhaustive controls |
| checks | Separate physical-geometry and polynomial-remainder implementations |

Run from the repository root:

    python -X utf8 -B scripts/verify_extended.py

Add --optimized for optimized Python subprocesses. Products are regenerated in a temporary tree and compared to the mathematical payloads. Only explicitly named elapsed-time fields are excluded. Uniform proof claims are justified in the manuscript, not by finite replay.

The retiming checker joins all 25 base bitsets directly to the atlas. Attachment identities and adaptations are in docs/EXTENSION_PROVENANCE.json. Serialization is UTF-8/LF with POSIX manifest keys.

The 130-path probe is restricted to its sampled paired paths. Its failure to improve is not a global negative result. Time-limited discovery logs are not mathematical nonexistence evidence and are omitted.
