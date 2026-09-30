#!/usr/bin/env python3
"""Verify the release attestation against current public bytes."""

from __future__ import annotations

from build_release_attestation import build
from package_common import ATTESTATION, load_json


def main() -> int:
    try:
        actual = load_json(ATTESTATION)
        expected = build()
    except (OSError, ValueError, KeyError) as exc:
        print(f"FAIL release_attestation:{type(exc).__name__}:{exc}")
        return 1
    if actual != expected:
        print("FAIL release_attestation_mismatch")
        return 1
    print(
        "PASS release attestation "
        f"entries={expected['manifest']['entries']} pdf={expected['paper']['sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
