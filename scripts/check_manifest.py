#!/usr/bin/env python3
"""Independently verify manifest hashes and exact closed-tree coverage."""

from __future__ import annotations

import argparse

from package_common import ROOT, distributable_paths, parse_manifest, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--closed-tree", action="store_true")
    parser.parse_args()
    failures: list[str] = []
    try:
        rows = parse_manifest()
    except (OSError, ValueError) as exc:
        print(f"FAIL {exc}")
        return 1
    seen = {row["path"] for row in rows}
    for row in rows:
        relative = row["path"]
        path = ROOT.joinpath(*relative.split("/"))
        if not path.is_file():
            failures.append(f"missing {relative}")
            continue
        if sha256_file(path) != row["sha256"]:
            failures.append(f"mismatch {relative}")
    expected_set = {path.relative_to(ROOT).as_posix() for path in distributable_paths()}
    for relative in sorted(expected_set - seen):
        failures.append(f"unmanifested {relative}")
    for relative in sorted(seen - expected_set):
        failures.append(f"unexpected manifest entry {relative}")
    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        return 1
    print(f"PASS manifest hashes and coverage {len(seen)}/{len(expected_set)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
