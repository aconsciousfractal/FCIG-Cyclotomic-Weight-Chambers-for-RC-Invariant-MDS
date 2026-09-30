#!/usr/bin/env python3
"""Build the closed-tree SHA-256 manifest."""

from __future__ import annotations

from package_common import MANIFEST, ROOT, distributable_paths, sha256_file


def main() -> int:
    paths = distributable_paths()
    lines = [
        "# SHA-256 manifest for FCIG-Cyclotomic-Weight-Chambers-for-RC-Invariant-MDS",
        "# RELEASE_ATTESTATION.json, this manifest, Git metadata, caches, and TeX auxiliaries are excluded.",
    ]
    for path in paths:
        lines.append(f"{sha256_file(path).upper()}  {path.relative_to(ROOT).as_posix()}")
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {MANIFEST.name} entries={len(paths)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
