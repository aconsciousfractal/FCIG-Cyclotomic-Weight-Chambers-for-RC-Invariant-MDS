#!/usr/bin/env python3
"""Reconstruct and validate a temporary tree solely from manifest entries."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from package_common import ATTESTATION, MANIFEST, ROOT, parse_manifest


def main() -> int:
    rows = parse_manifest()
    with tempfile.TemporaryDirectory(prefix="fcig-cwc-manifest-") as temporary:
        target = Path(temporary) / "package"
        target.mkdir()
        for row in rows:
            source = ROOT.joinpath(*row["path"].split("/"))
            destination = target.joinpath(*row["path"].split("/"))
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        shutil.copy2(MANIFEST, target / MANIFEST.name)
        shutil.copy2(ATTESTATION, target / ATTESTATION.name)
        for script, arguments in (
            ("check_manifest.py", ["--closed-tree"]),
            ("check_release_attestation.py", []),
            ("validate_repository.py", []),
        ):
            completed = subprocess.run(
                [sys.executable, "-X", "utf8", "-B", str(target / "scripts" / script), *arguments],
                cwd=target, text=True, capture_output=True, check=False, timeout=120,
            )
            if completed.returncode:
                print(f"FAIL manifest_only:{script}:{completed.stdout}{completed.stderr}")
                return 1
    print(f"PASS manifest-only reconstruction entries={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
