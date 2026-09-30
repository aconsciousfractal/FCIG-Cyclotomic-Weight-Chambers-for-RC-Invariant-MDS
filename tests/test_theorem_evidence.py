from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "evidence" / "theorem" / "results" / "single_odd_prime_embedding_checks.json"


def canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def test_single_odd_prime_embedding_result() -> None:
    raw = RESULT.read_bytes()
    value = json.loads(raw)
    assert raw == canonical(value) + b"\n"
    payload = value["scientific_payload"]
    assert hashlib.sha256(canonical(payload)).hexdigest() == value["scientific_payload_sha256"]
    assert value["status"] == "PASS_SINGLE_ODD_PRIME_EMBEDDING_FINITE_ASSURANCE"
    assert len(payload["cyclotomic_composition_cases"]) == 18
    assert [row["k"] for row in payload["positive_binary_word_cases"]] == [6, 10, 12, 18, 20, 24]
    assert all(row["observed"] == row["expected"] for row in payload["positive_binary_word_cases"])
    assert all(row["rejected"] is True for row in payload["negative_controls"])
    assert payload["orbit_controls"] == [
        {"k": 12, "observed_pcrs": 19, "predicted_pcrs": 19},
        {"k": 30, "observed_zero_words": 146854, "predicted_pcrs": 6069},
    ]
