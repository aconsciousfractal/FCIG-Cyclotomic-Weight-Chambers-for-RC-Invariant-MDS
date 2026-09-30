#!/usr/bin/env python3
"""Shared deterministic helpers for the public package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST_SHA256.txt"
ATTESTATION = ROOT / "RELEASE_ATTESTATION.json"
PDF_RELATIVE = "paper/Cyclotomic_Weight_Chambers_for_RC_Invariant_MDS.pdf"
PDF = ROOT.joinpath(*PurePosixPath(PDF_RELATIVE).parts)
BUILD_SUFFIXES = (
    ".aux", ".bbl", ".blg", ".log", ".out", ".toc", ".fls",
    ".fdb_latexmk", ".synctex.gz",
)
EXCLUDED_DIRS = {".git", "__pycache__", ".pytest_cache", "dist"}


def canonical(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def distributable(path: Path) -> bool:
    if not path.is_file() or path in {MANIFEST, ATTESTATION}:
        return False
    relative = path.relative_to(ROOT)
    if any(part in EXCLUDED_DIRS for part in relative.parts):
        return False
    return not any(relative.name.endswith(suffix) for suffix in BUILD_SUFFIXES)


def distributable_paths() -> tuple[Path, ...]:
    return tuple(sorted(
        (path for path in ROOT.rglob("*") if distributable(path)),
        key=lambda path: path.relative_to(ROOT).as_posix(),
    ))


def no_duplicate_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate_json_key:{key}")
        result[key] = value
    return result


def load_json(path: Path, *, canonical_lf: bool = False) -> object:
    raw = path.read_bytes()
    value = json.loads(raw, object_pairs_hook=no_duplicate_object)
    if canonical_lf and raw != canonical(value) + b"\n":
        raise ValueError(f"noncanonical_json:{path.relative_to(ROOT).as_posix()}")
    return value


def parse_manifest() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for number, raw in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            digest, relative = line.split(None, 1)
        except ValueError as exc:
            raise ValueError(f"manifest_malformed_line:{number}") from exc
        posix = PurePosixPath(relative)
        if posix.is_absolute() or ".." in posix.parts or relative in seen:
            raise ValueError(f"manifest_unsafe_or_duplicate:{relative}")
        if len(digest) != 64 or any(char not in "0123456789abcdefABCDEF" for char in digest):
            raise ValueError(f"manifest_bad_digest:{relative}")
        seen.add(relative)
        rows.append({"path": relative, "sha256": digest.lower()})
    return rows
