from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(script: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / script), *arguments],
        cwd=ROOT, text=True, capture_output=True, check=False, timeout=120,
    )


def test_manifest_is_closed() -> None:
    result = run("check_manifest.py", "--closed-tree")
    assert result.returncode == 0, result.stdout + result.stderr


def test_release_attestation() -> None:
    result = run("check_release_attestation.py")
    assert result.returncode == 0, result.stdout + result.stderr


def test_public_surface() -> None:
    result = run("validate_repository.py")
    assert result.returncode == 0, result.stdout + result.stderr
