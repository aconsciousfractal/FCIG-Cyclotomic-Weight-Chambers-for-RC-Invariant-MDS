#!/usr/bin/env python3
"""Re-seal and reject hostile public-package mutations."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parent.parent
Mutation = Callable[[Path], None]
BS = chr(92)


def write_text(root: Path, relative: str, value: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def append(root: Path, relative: str, data: bytes) -> None:
    path = root / relative
    path.write_bytes(path.read_bytes() + data)


def alter_atlas(root: Path) -> None:
    path = root / "evidence" / "finite_assurance" / "results" / "certified_terminal_atlas.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    value["scientific_payload_sha256"] = "0" * 64
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")


def mutations() -> list[tuple[str, Mutation]]:
    return [
        (
            "yaml_private_path",
            lambda root: write_text(
                root, "README.md", (root / "README.md").read_text(encoding="utf-8") + "\nsource: " + "C:" + BS + "Temp" + BS + "private.txt\n"
            ),
        ),
        (
            "extensionless_private_path",
            lambda root: write_text(root, "README.md", (root / "README.md").read_text(encoding="utf-8") + "\nsource: /" + "home/reviewer/private.txt\n"),
        ),
        (
            "unix_workspace_path",
            lambda root: write_text(root, "evidence/rogue.json", json.dumps({"base": "/mnt/" + "data/private.json"}) + "\n"),
        ),
        (
            "environment_secret",
            lambda root: write_text(root, "README.md", (root / "README.md").read_text(encoding="utf-8") + "\nAPI_TOKEN=" + "sk-" + "A" * 32 + "\n"),
        ),
        (
            "unc_private_path",
            lambda root: write_text(
                root,
                "README.md",
                (root / "README.md").read_text(encoding="utf-8")
                + "\nsource=" + BS + BS + "private-host" + BS + "hidden-share" + BS + "input.txt\n",
            ),
        ),
        (
            "raw_transcript",
            lambda root: write_text(root, "README.md", (root / "README.md").read_text(encoding="utf-8") + "\nSYSTEM" + ": hidden instructions\n"),
        ),
        (
            "tar_archive", lambda root: (root / "rogue.tar").write_bytes(b"archive")
        ),
        (
            "invalid_utf8_text",
            lambda root: append(root, "README.md", b"\n\xff\xfe"),
        ),
        (
            "additional_pdf",
            lambda root: (root / "private.pdf").write_bytes(
                (root / "paper" / "Cyclotomic_Weight_Chambers_for_RC_Invariant_MDS.pdf").read_bytes()
            ),
        ),
        (
            "active_pdf_marker",
            lambda root: append(root, "paper/Cyclotomic_Weight_Chambers_for_RC_Invariant_MDS.pdf", b"\n/JavaScript\n"),
        ),
        ("atlas_payload_tamper", alter_atlas),
        ("missing_claim_boundary", lambda root: (root / "docs" / "PUBLIC_CLAIM_BOUNDARY.md").unlink()),
        ("unexpected_binary", lambda root: (root / "rogue.bin").write_bytes(b"binary")),
    ]


def run(script: Path, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script)],
        cwd=cwd,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )


def main() -> int:
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="fcig_cwc_hygiene_") as temporary:
        scratch = Path(temporary)
        for index, (name, mutate) in enumerate(mutations()):
            work = scratch / f"case_{index:02d}"
            shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            mutate(work)
            manifest_build = run(work / "scripts" / "build_manifest.py", work)
            release_build = run(work / "scripts" / "build_release_attestation.py", work)
            manifest_check = run(work / "scripts" / "check_manifest.py", work)
            repository_check = run(work / "scripts" / "validate_repository.py", work)
            sealed = manifest_build.returncode == 0 and release_build.returncode == 0 and manifest_check.returncode == 0
            rejected = sealed and repository_check.returncode != 0
            rows.append(
                {
                    "name": name,
                    "rejected": rejected,
                    "manifest_exit": manifest_check.returncode,
                    "repository_exit": repository_check.returncode,
                }
            )

    passed = len(rows) == len(mutations()) and all(bool(row["rejected"]) for row in rows)
    result = {
        "status": "PASS_RELEASE_HYGIENE_MUTATIONS" if passed else "FAIL_RELEASE_HYGIENE_MUTATIONS",
        "pass": passed,
        "mutations_rejected": sum(bool(row["rejected"]) for row in rows),
        "mutations_total": len(rows),
        "rows": rows,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
