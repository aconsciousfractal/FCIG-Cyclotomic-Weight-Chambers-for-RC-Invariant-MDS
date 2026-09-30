#!/usr/bin/env python3
"""Build the acyclic public release attestation."""

from __future__ import annotations

import json

from package_common import (
    ATTESTATION, MANIFEST, PDF, PDF_RELATIVE, ROOT, canonical, load_json,
    parse_manifest, sha256_bytes, sha256_file,
)


ATLAS = ROOT / "evidence" / "finite_assurance" / "results" / "certified_terminal_atlas.json"
VALIDATION = ROOT / "evidence" / "finite_assurance" / "results" / "certified_validation.json"
SINGLE_ODD_PRIME = ROOT / "evidence" / "theorem" / "results" / "single_odd_prime_embedding_checks.json"


def build() -> dict[str, object]:
    rows = parse_manifest()
    atlas = load_json(ATLAS, canonical_lf=True)
    validation = load_json(VALIDATION, canonical_lf=True)
    single_odd_prime = load_json(SINGLE_ODD_PRIME, canonical_lf=True)
    assert isinstance(atlas, dict) and isinstance(validation, dict)
    assert isinstance(single_odd_prime, dict)
    return {
        "finite_assurance": {
            "atlas_payload_sha256": atlas["scientific_payload_sha256"],
            "single_odd_prime_embedding_payload_sha256": single_odd_prime["scientific_payload_sha256"],
            "validation_payload_sha256": validation["scientific_payload_sha256"],
            "advances_provenance_sha256": sha256_file(ROOT / "evidence/advances/PROVENANCE.json"),
            "balanced_loads_provenance_sha256": sha256_file(ROOT / "evidence/balanced_loads/PROVENANCE.json"),
        },
        "manifest": {
            "entries": len(rows),
            "sha256": sha256_file(MANIFEST),
            "tree_sha256": sha256_bytes(canonical(rows)),
        },
        "paper": {
            "bytes": PDF.stat().st_size,
            "path": PDF_RELATIVE,
            "sha256": sha256_file(PDF),
        },
        "publication_status": "OUT_OF_SCOPE_NOT_ASSERTED",
        "schema": "fcig-cwc-release-attestation-v1",
        "status": "PASS_STANDALONE_RELEASE_ATTESTATION",
        "version": "3.2.0",
    }


def main() -> int:
    value = build()
    ATTESTATION.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(f"wrote {ATTESTATION.name} manifest_entries={value['manifest']['entries']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
