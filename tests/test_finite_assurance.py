from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_HASH = "aca5d5fc81c839c218b6eefaaf00acdc733b6e86c342e3d13546c82b0ce07f1f"
VALIDATION_HASH = "b79051b4190d6f13b9bb183fe1299fb67b34a1eef21bbed5bf09c516f243ce46"


def canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def load(relative: str) -> dict[str, object]:
    path = ROOT.joinpath(*relative.split("/"))
    raw = path.read_bytes()
    value = json.loads(raw)
    assert raw == canonical(value) + b"\n"
    return value


def test_atlas_identity_and_scope() -> None:
    atlas = load("evidence/finite_assurance/results/certified_terminal_atlas.json")
    payload = atlas["scientific_payload"]
    assert hashlib.sha256(canonical(payload)).hexdigest() == ATLAS_HASH
    assert atlas["scientific_payload_sha256"] == ATLAS_HASH
    assert len(payload["rows"]) == 30
    assert [len(order["ordered_row_ids"]) for order in payload["orders"]] == [3, 23, 4]
    assert [order["unique_terminal_geometry_count"] for order in payload["orders"]] == [2, 17, 4]
    assert payload["assurance_scope"]["complete_k8"] is False


def test_validation_identity_and_mutations() -> None:
    validation = load("evidence/finite_assurance/results/certified_validation.json")
    payload = validation["scientific_payload"]
    assert hashlib.sha256(canonical(payload)).hexdigest() == VALIDATION_HASH
    assert validation["scientific_payload_sha256"] == VALIDATION_HASH
    assert payload["validated_row_count"] == 30
    mutations = payload["mutation_panel"]["mutations"]
    assert len(mutations) == 29
    assert all(row["rejected"] is True for row in mutations)
    by_id = {row["mutation_id"]: row for row in mutations}
    assert by_id["critical_ratio_skeleton_root_substitution"] == {
        "error_code": "ATLAS_SKELETON_ORDER_K6",
        "mutation_id": "critical_ratio_skeleton_root_substitution",
        "rejected": True,
    }


def test_every_row_is_exact_and_acyclic() -> None:
    atlas = load("evidence/finite_assurance/results/certified_terminal_atlas.json")
    for row in atlas["scientific_payload"]["rows"]:
        terminal = row["terminal"]
        residual = row["residual_dag"]
        assert terminal["selected_word_count"] == terminal["minimum_cardinality"]
        assert terminal["rc_mismatch_count"] == 0
        assert terminal["pcr_underselection_count"] == 0
        assert terminal["pcr_overselection_count"] == 0
        assert terminal["direct_firing_mismatch_count"] == 0
        assert residual["kahn_acyclic"] is True
        assert residual["dfs_acyclic"] is True
        assert residual["longest_path_vertices"] == residual["longest_path_edges"] + 1
