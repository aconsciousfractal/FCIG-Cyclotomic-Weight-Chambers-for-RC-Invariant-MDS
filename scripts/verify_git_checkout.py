#!/usr/bin/env python3
"""Verify clean committed bytes through a safe Git archive materialization."""

from __future__ import annotations

import io
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

from package_common import ROOT


def git(*args: str, text: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ROOT, text=text, capture_output=True, check=False)


def main() -> int:
    status = git("status", "--porcelain")
    if status.returncode or status.stdout.strip():
        print("FAIL git_checkout_dirty_or_unavailable")
        return 1
    commits = git("rev-list", "--all")
    roots = git("rev-list", "--max-parents=0", "--all")
    commit_rows = [row for row in commits.stdout.splitlines() if row]
    root_rows = [row for row in roots.stdout.splitlines() if row]
    if not commit_rows or len(root_rows) != 1:
        print(f"FAIL git_history commits={len(commit_rows)} roots={len(root_rows)}")
        return 1
    archive = git("archive", "--format=tar", "HEAD", text=False)
    if archive.returncode:
        print("FAIL git_archive")
        return 1
    with tempfile.TemporaryDirectory(prefix="fcig-cwc-git-") as temporary:
        target = Path(temporary) / "checkout"
        target.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive.stdout), mode="r:") as bundle:
            for member in bundle.getmembers():
                path = PurePosixPath(member.name)
                if path.is_absolute() or ".." in path.parts or member.issym() or member.islnk():
                    print(f"FAIL unsafe_git_archive_member:{member.name}")
                    return 1
            bundle.extractall(target, filter="data")
        completed = subprocess.run(
            [sys.executable, "-X", "utf8", "-B", str(target / "scripts" / "verify.py"), "--profile", "manifest"],
            cwd=target, text=True, capture_output=True, check=False, timeout=180,
        )
        if completed.returncode:
            print(f"FAIL git_materialized_manifest:{completed.stdout}{completed.stderr}")
            return 1
    history = subprocess.run(
        [sys.executable, "-X", "utf8", "-B", str(ROOT / "scripts" / "check_git_history.py")],
        cwd=ROOT, text=True, capture_output=True, check=False, timeout=600,
    )
    if history.returncode:
        print(f"FAIL git_history_privacy:{history.stdout}{history.stderr}")
        return 1
    print("PASS clean single-root Git checkout and archive materialization")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
