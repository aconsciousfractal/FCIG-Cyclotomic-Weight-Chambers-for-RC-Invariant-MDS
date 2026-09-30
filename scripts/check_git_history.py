#!/usr/bin/env python3
"""Fail closed unless the distributable Git history is a clean history with one root."""

from __future__ import annotations

import re
import gzip
import subprocess
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parent.parent
FORBIDDEN_BASENAMES = {
    ".env", "id_dsa", "id_ecdsa", "id_ed25519", "id_rsa",
    "secrets.json", "credentials.json", "private.key",
}
WINDOWS_PATH = re.compile(
    r"(?i)\b[A-Z]:[\\/](?:[A-Za-z0-9_. -]{2,}[\\/]"
    r"|[A-Za-z0-9_. -]+\.[A-Za-z0-9]{1,8}(?:\s|$))"
)
UNC_PATH = re.compile(r"\\\\[A-Za-z0-9._$-]{2,}[\\/][A-Za-z0-9._$-]{2,}")
UNIX_PRIVATE_PATH = re.compile(r"/(?:Users|home)/[^/\s]+/")
UNIX_WORKSPACE_PATH = re.compile(re.escape("/mnt/" + "data/"))
PRIVATE_TEXT = {
    "private_key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    "github_token": re.compile(r"ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}"),
    "api_token": re.compile(r"sk-[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}"),
    "raw_transcript": re.compile(r"(?m)^\s*(?:SYSTEM|USER|ASSISTANT)\s*:"),
}
TEXT_SUFFIXES = {
    ".md", ".txt", ".tex", ".bib", ".py", ".ps1", ".cff", ".json",
    ".yaml", ".yml", ".toml", ".ini", ".cfg",
}
TEXT_NAMES = {".gitignore", ".gitattributes", "LICENSE"}


def git(*args: str, text: bool = True) -> str | bytes:
    completed = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        text=text,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() if text else completed.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git {' '.join(args)} failed: {detail}")
    return completed.stdout


def main() -> None:
    failures: list[str] = []
    try:
        commits = [line for line in str(git("rev-list", "--all")).splitlines() if line]
        roots = [
            line
            for line in str(git("rev-list", "--max-parents=0", "--all")).splitlines()
            if line
        ]
        if not commits:
            failures.append(f"reachable commit count is {len(commits)}, expected at least 1")
        if len(roots) != 1:
            failures.append(f"root commit count is {len(roots)}, expected 1")

        object_rows = [
            line for line in str(git("rev-list", "--objects", "--all")).splitlines() if line
        ]
        seen_objects: set[str] = set()
        for row in object_rows:
            parts = row.split(" ", 1)
            object_id = parts[0]
            path = parts[1] if len(parts) == 2 else None
            seen_objects.add(object_id)
            if path is not None and PurePosixPath(path).name.lower() in FORBIDDEN_BASENAMES:
                failures.append(f"private historical path {path}")

            object_type = str(git("cat-file", "-t", object_id)).strip()
            if object_type != "blob":
                continue
            if path is None:
                continue
            posix = PurePosixPath(path)
            compressed = path.endswith(".json.gz")
            if not compressed and posix.suffix.lower() not in TEXT_SUFFIXES and posix.name not in TEXT_NAMES:
                continue
            data = bytes(git("cat-file", "blob", object_id, text=False))
            if compressed:
                data = gzip.decompress(data)
            decoded = data.decode("utf-8", errors="ignore")
            for name, pattern in PRIVATE_TEXT.items():
                if pattern.search(decoded):
                    failures.append(f"{name} in historical blob {path or object_id}")
            if WINDOWS_PATH.search(decoded):
                failures.append(f"machine-local path in historical blob {path or object_id}")
            if UNC_PATH.search(decoded):
                failures.append(f"UNC path in historical blob {path or object_id}")
            if UNIX_PRIVATE_PATH.search(decoded):
                failures.append(f"user-home path in historical blob {path or object_id}")
            if UNIX_WORKSPACE_PATH.search(decoded):
                failures.append(f"workspace path in historical blob {path or object_id}")
    except Exception as exc:
        failures.append(f"history inspection failed: {type(exc).__name__}: {exc}")
        object_rows = []
        seen_objects = set()

    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        raise SystemExit(1)
    print(
        "PASS Git history single-root "
        f"commits={len(commits)} reachable_objects={len(seen_objects)} listed_rows={len(object_rows)}"
    )


if __name__ == "__main__":
    main()
