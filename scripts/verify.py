#!/usr/bin/env python3
"""Orchestrate manifest, exact-replay, and isolated verification profiles."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from package_common import ROOT, load_json


def run(command: list[str], cwd: Path, timeout: int = 900) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, cwd=cwd, text=True, encoding="utf-8", errors="replace",
        capture_output=True, timeout=timeout, check=False,
    )


def require_process(name: str, command: list[str], cwd: Path, rows: list[dict[str, object]], timeout: int = 900) -> None:
    completed = run(command, cwd, timeout)
    rows.append({
        "check": name,
        "passed": completed.returncode == 0,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip()[-2000:],
        "stderr": completed.stderr.strip()[-2000:],
    })


def manifest_checks(rows: list[dict[str, object]]) -> None:
    for script, args in (
        ("build_evidence_identity.py", ["--check"]),
        ("check_manifest.py", ["--closed-tree"]),
        ("check_release_attestation.py", []),
        ("validate_repository.py", []),
    ):
        require_process(script, [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / script), *args], ROOT, rows)


def core_checks(rows: list[dict[str, object]]) -> None:
    shipped_atlas = ROOT / "evidence" / "finite_assurance" / "results" / "certified_terminal_atlas.json"
    shipped_validation = ROOT / "evidence" / "finite_assurance" / "results" / "certified_validation.json"
    with tempfile.TemporaryDirectory(prefix="fcig-cwc-core-") as temporary:
        capsule = Path(temporary) / "finite_assurance"
        shutil.copytree(ROOT / "evidence" / "finite_assurance", capsule)
        first = Path(temporary) / "atlas_1.json"
        second = Path(temporary) / "atlas_2.json"
        normal = Path(temporary) / "validation.json"
        optimized = Path(temporary) / "validation_O.json"
        producer = capsule / "scripts" / "produce.py"
        validator = capsule / "scripts" / "validate.py"
        require_process("producer_first", [sys.executable, "-X", "utf8", "-B", str(producer), "--output", str(first)], capsule, rows)
        require_process("producer_optimized", [sys.executable, "-O", "-X", "utf8", "-B", str(producer), "--output", str(second)], capsule, rows)
        if first.is_file() and second.is_file():
            rows.append({"check": "producer_byte_identity", "passed": first.read_bytes() == second.read_bytes() == shipped_atlas.read_bytes()})
        require_process("validator_normal", [sys.executable, "-X", "utf8", "-B", str(validator), "--receipt", str(first), "--output", str(normal)], capsule, rows)
        require_process("validator_optimized", [sys.executable, "-O", "-X", "utf8", "-B", str(validator), "--receipt", str(first), "--output", str(optimized)], capsule, rows)
        if normal.is_file() and optimized.is_file():
            rows.append({"check": "validator_byte_identity", "passed": normal.read_bytes() == optimized.read_bytes() == shipped_validation.read_bytes()})
            value = load_json(normal, canonical_lf=True)
            assert isinstance(value, dict)
            payload = value.get("scientific_payload", {})
            rows.append({
                "check": "mutation_panel_29_of_29",
                "passed": payload.get("mutation_panel", {}).get("mutation_count") == 29
                and all(item.get("rejected") is True for item in payload.get("mutation_panel", {}).get("mutations", [])),
            })
        theorem_normal = Path(temporary) / "single_odd_prime.json"
        theorem_optimized = Path(temporary) / "single_odd_prime_O.json"
        theorem_script = ROOT / "evidence" / "theorem" / "verify_single_odd_prime_embedding.py"
        theorem_shipped = ROOT / "evidence" / "theorem" / "results" / "single_odd_prime_embedding_checks.json"
        require_process(
            "single_odd_prime_normal",
            [sys.executable, "-X", "utf8", "-B", str(theorem_script), "--output", str(theorem_normal)],
            ROOT, rows,
        )
        require_process(
            "single_odd_prime_optimized",
            [sys.executable, "-O", "-X", "utf8", "-B", str(theorem_script), "--output", str(theorem_optimized)],
            ROOT, rows,
        )
        if theorem_normal.is_file() and theorem_optimized.is_file():
            rows.append({
                "check": "single_odd_prime_byte_identity",
                "passed": theorem_normal.read_bytes() == theorem_optimized.read_bytes() == theorem_shipped.read_bytes(),
            })
    require_process(
        "pytest",
        [sys.executable, "-X", "utf8", "-B", "-m", "pytest", "-q", "-p", "no:cacheprovider"],
        ROOT, rows,
    )

    require_process(
        "extended_exact_suite",
        [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / "verify_extended.py")],
        ROOT, rows,
    )
    require_process(
        "advances_exact_suite",
        [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / "verify_advances.py")],
        ROOT, rows, timeout=2400,
    )
    for optimized in (False, True):
        require_process(
            "balanced_loads_optimized" if optimized else "balanced_loads_normal",
            [sys.executable, *(["-O"] if optimized else []), "-X", "utf8", "-B",
             str(ROOT / "scripts" / "verify_balanced_loads.py")], ROOT, rows,
        )


def full_checks(rows: list[dict[str, object]]) -> None:
    for script in ("verify_manifest_only.py", "verify_git_checkout.py"):
        require_process(script, [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / script)], ROOT, rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("manifest", "core", "full"), default="core")
    args = parser.parse_args()
    rows: list[dict[str, object]] = []
    manifest_checks(rows)
    if args.profile in {"core", "full"} and all(bool(row.get("passed")) for row in rows):
        core_checks(rows)
    if args.profile == "full" and all(bool(row.get("passed")) for row in rows):
        full_checks(rows)
    passed = bool(rows) and all(bool(row.get("passed")) for row in rows)
    result = {
        "checks": rows,
        "checks_passed": sum(bool(row.get("passed")) for row in rows),
        "checks_total": len(rows),
        "pass": passed,
        "profile": args.profile,
        "schema": "fcig-cwc-verification-v1",
        "status": "PASS_STANDALONE_PACKAGE_VERIFICATION" if passed else "FAIL_STANDALONE_PACKAGE_VERIFICATION",
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
