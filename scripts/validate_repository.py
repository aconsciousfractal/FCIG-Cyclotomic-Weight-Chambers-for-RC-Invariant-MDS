#!/usr/bin/env python3
"""Fail-closed validation of the public repository surface."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

from pypdf import PdfReader
from pypdf.generic import ArrayObject, DictionaryObject, IndirectObject

from package_common import PDF, PDF_RELATIVE, ROOT, load_json, sha256_file


REQUIRED = {
    ".gitattributes", ".gitignore", ".github/workflows/verify.yml",
    "CITATION.cff", "INDEPENDENT_VERIFICATION.md", "LICENSE",
    "LICENSE_SCOPE.md", "MANIFEST_SHA256.txt", "README.md",
    "RELEASE_ATTESTATION.json", "REPRODUCE.md", "THIRD_PARTY_NOTICES.md",
    "docs/ACCESSIBILITY.md", "docs/ARTIFACT_MAP.md", "docs/AUTHORSHIP.md",
    "docs/PRIOR_ART_BOUNDARY.md", "docs/PUBLIC_CLAIM_BOUNDARY.md",
    "docs/SOURCE_LOCK.md", "docs/THEOREM_MAP.json",
    "docs/VERIFICATION_ARCHITECTURE.md",
    "evidence/finite_assurance/README.md",
    "evidence/finite_assurance/VALIDATION_PROTOCOL.md",
    "evidence/finite_assurance/inputs/critical_ratio_skeleton.json",
    "evidence/finite_assurance/results/certified_terminal_atlas.json",
    "evidence/finite_assurance/results/certified_validation.json",
    "evidence/finite_assurance/schemas/atlas_row.schema.json",
    "evidence/finite_assurance/schemas/terminal_atlas.schema.json",
    "evidence/finite_assurance/scripts/produce.py",
    "evidence/finite_assurance/scripts/validate.py",
    "evidence/theorem/README.md",
    "evidence/theorem/results/single_odd_prime_embedding_checks.json",
    "evidence/theorem/verify_single_odd_prime_embedding.py",
    "paper/BUILD.md", "paper/abstract.tex",
    "paper/macros.tex", "paper/main.tex", "paper/refs.bib", PDF_RELATIVE,
    "requirements.txt", "scripts/build_manifest.py",
    "scripts/build_paper.py", "scripts/build_release_attestation.py",
    "scripts/check_git_history.py", "scripts/check_manifest.py",
    "scripts/build_evidence_identity.py", "scripts/manuscript_sources.py",
    "paper/evidence_identity.tex", "results/evidence_identity.json",
    "scripts/check_release_attestation.py", "scripts/package_common.py",
    "scripts/test_release_hygiene_mutations.py", "scripts/validate_repository.py",
    "scripts/verify.py", "scripts/verify_git_checkout.py",
    "scripts/verify_manifest_only.py", "tests/test_finite_assurance.py",
    "tests/test_repository.py", "tests/test_theorem_evidence.py",
    "docs/ADVANCES.md", "docs/OPEN_PROBLEMS.md",
    "evidence/advances/README.md", "evidence/advances/PROVENANCE.json",
    "scripts/verify_advances.py", "pytest.ini",
    "evidence/balanced_loads/README.md", "evidence/balanced_loads/PROVENANCE.json",
    "evidence/balanced_loads/verify_case.py", "evidence/balanced_loads/MANIFEST_SHA256.txt",
    "scripts/verify_balanced_loads.py", "tests/test_balanced_loads.py",
}
TEXT_SUFFIXES = {
    ".bib", ".cff", ".json", ".md", ".py", ".tex", ".txt", ".yaml", ".yml",
}
TEXT_NAMES = {".gitattributes", ".gitignore", "LICENSE", "pytest.ini"}
FORBIDDEN_ENDINGS = (
    ".aux", ".bbl", ".blg", ".log", ".out", ".toc", ".pyc", ".pyo",
    ".fls", ".fdb_latexmk", ".synctex.gz", ".zip", ".7z", ".tar",
    ".tar.gz", ".tgz", ".rar",
)
PRIVATE_PATTERNS = {
    "windows_absolute_path": re.compile(
        r"(?i)\b[A-Z]:[\\/](?:[A-Za-z0-9_. -]{2,}[\\/]"
        r"|[A-Za-z0-9_. -]+\.[A-Za-z0-9]{1,8}(?:\s|$))"
    ),
    "windows_unc_path": re.compile(
        re.escape("\\" + "\\")
        + r"[A-Za-z0-9._$-]{2,}[\\/][A-Za-z0-9._$-]{2,}"
    ),
    "unix_private_path": re.compile(r"/(?:Users|home)/[^/\s]+/"),
    "unix_workspace_path": re.compile(re.escape("/mnt/" + "data/")),
    "private_key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    "github_token": re.compile(r"ghp_[A-Za-z0-9]{20,}|github_" + r"pat_[A-Za-z0-9_]{20,}"),
    "api_token": re.compile(r"sk-[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}"),
    "raw_transcript": re.compile(r"(?m)^\s*(?:SYSTEM|USER|ASSISTANT)\s*:"),
}
PDF_TITLE = "Cyclotomic Weight Chambers for Reverse-Complement-Invariant Minimum Decycling Sets"
ATLAS_HASH = "aca5d5fc81c839c218b6eefaaf00acdc733b6e86c342e3d13546c82b0ce07f1f"
VALIDATION_HASH = "b79051b4190d6f13b9bb183fe1299fb67b34a1eef21bbed5bf09c516f243ce46"


def check(name: str, passed: bool, detail: object, rows: list[dict[str, object]]) -> None:
    rows.append({"check": name, "detail": detail, "passed": bool(passed)})


def tex_files() -> tuple[Path, ...]:
    from manuscript_sources import active_tex_files
    return active_tex_files(ROOT / "paper")


def manuscript_surface(rows: list[dict[str, object]]) -> None:
    source = "\n".join(path.read_text(encoding="utf-8") for path in tex_files())
    labels = set(re.findall(r"\\label\{([^}]+)\}", source))
    theorem_map = load_json(ROOT / "docs" / "THEOREM_MAP.json")
    assert isinstance(theorem_map, dict)
    claims = theorem_map.get("claims", [])
    wanted = {row["label"] for row in claims if isinstance(row, dict)}
    imported = theorem_map.get("imported_interface", {})
    if isinstance(imported, dict):
        wanted.add(imported["label"])
    check("theorem_label_map", wanted <= labels, {"mapped": len(wanted), "labels": len(labels)}, rows)

    citation_keys: set[str] = set()
    for group in re.findall(r"\\cite\w*(?:\[[^]]*\])?\{([^}]+)\}", source):
        citation_keys.update(item.strip() for item in group.split(",") if item.strip())
    bibliography = (ROOT / "paper" / "refs.bib").read_text(encoding="utf-8")
    bibliography_keys = set(re.findall(r"@[A-Za-z]+\{([^,]+),", bibliography))
    check(
        "bibliography_equals_citations",
        citation_keys == bibliography_keys and bool(citation_keys),
        {"citations": len(citation_keys), "bibliography": len(bibliography_keys)}, rows,
    )
    formal_pattern = re.compile(
        r"\\begin\{(theorem|proposition|lemma|corollary|definition|remark|example)\}"
        r"(?:\[[^]]*\])?(.*?)\\end\{\1\}",
        flags=re.DOTALL,
    )
    actual: list[dict[str, str | None]] = []
    for match in formal_pattern.finditer(source):
        environment_labels = [
            label for label in re.findall(r"\\label\{([^}]+)\}", match.group(2))
            if not label.startswith("eq:")
        ]
        actual.append({
            "kind": match.group(1),
            "label": environment_labels[0] if environment_labels else None,
        })
    inventory = theorem_map.get("formal_inventory", [])
    expected = [
        {"kind": row.get("kind"), "label": row.get("label")}
        for row in inventory if isinstance(row, dict)
    ]
    check(
        "formal_environment_inventory",
        actual == expected and len(actual) == theorem_map["formal_count"],
        {"actual": actual, "expected": expected},
        rows,
    )


DANGEROUS_PDF_KEYS = {
    "/AA", "/EmbeddedFile", "/JavaScript", "/JS", "/Launch",
    "/RichMedia", "/RichMediaExecute",
}
DANGEROUS_PDF_ACTIONS = {
    "/GoToR", "/ImportData", "/JavaScript", "/Launch", "/Movie",
    "/Rendition", "/RichMediaExecute", "/Sound", "/SubmitForm",
}


def pdf_object_findings(reader: PdfReader) -> tuple[list[str], list[str]]:
    findings: list[str] = []
    uris: list[str] = []
    visited_indirect: set[tuple[int, int]] = set()
    visited_direct: set[int] = set()

    def resolve(value: object) -> object:
        return value.get_object() if isinstance(value, IndirectObject) else value

    def walk(value: object, location: str) -> None:
        if isinstance(value, IndirectObject):
            identity = (value.idnum, value.generation)
            if identity in visited_indirect:
                return
            visited_indirect.add(identity)
            value = value.get_object()
        if isinstance(value, (DictionaryObject, ArrayObject)):
            identity = id(value)
            if identity in visited_direct:
                return
            visited_direct.add(identity)
        if isinstance(value, DictionaryObject):
            for raw_key, raw_child in value.items():
                key = str(raw_key)
                child_location = f"{location}/{key.lstrip('/')}"
                if key in DANGEROUS_PDF_KEYS:
                    findings.append(f"dangerous_key:{child_location}")
                if key == "/OpenAction":
                    open_action = resolve(raw_child)
                    if isinstance(open_action, DictionaryObject):
                        action_type = str(resolve(open_action.get("/S", "")))
                        if action_type != "/GoTo":
                            findings.append(f"dangerous_open_action:{child_location}:{action_type}")
                    elif not isinstance(open_action, ArrayObject):
                        findings.append(f"invalid_open_action:{child_location}")
                if key == "/S":
                    action = str(resolve(raw_child))
                    if action in DANGEROUS_PDF_ACTIONS:
                        findings.append(f"dangerous_action:{child_location}:{action}")
                if key == "/URI":
                    uri = str(resolve(raw_child))
                    uris.append(uri)
                    if not re.fullmatch(r"https?://[^\s]+", uri):
                        findings.append(f"non_http_uri:{child_location}:{uri}")
                walk(raw_child, child_location)
        elif isinstance(value, ArrayObject):
            for index, child in enumerate(value):
                walk(child, f"{location}[{index}]")

    walk(reader.trailer, "trailer")
    return sorted(set(findings)), sorted(set(uris))


def pdf_surface(rows: list[dict[str, object]]) -> None:
    reader = PdfReader(PDF)
    metadata = reader.metadata or {}
    title = str(metadata.get("/Title", ""))
    author = str(metadata.get("/Author", ""))
    subject = str(metadata.get("/Subject", ""))
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    build_receipt = load_json(ROOT / "results" / "paper_build.json")
    assert isinstance(build_receipt, dict)
    expected_pages = build_receipt.get("pages")
    check("pdf_page_count", len(reader.pages) == expected_pages, {"pages": len(reader.pages), "expected": expected_pages}, rows)
    check("pdf_title", title == PDF_TITLE, {"title": title}, rows)
    check("pdf_author", author == "Oleksiy Babanskyy", {"author": author}, rows)
    check(
        "pdf_subject_public",
        bool(subject) and not re.search(r"(?i)draft|internal|pending|review manuscript", subject),
        {"subject": subject}, rows,
    )
    check(
        "pdf_text_markers",
        "certified 141/107 component separation" in text
        and "Main arithmetic results" in text
        and "evidence availability" in text,
        {"characters": len(text)}, rows,
    )
    text_privacy = [name for name, pattern in PRIVATE_PATTERNS.items() if pattern.search(text)]
    check("pdf_text_privacy", not text_privacy, {"findings": text_privacy}, rows)
    check(
        "pdf_reference_links_visible",
        "doi.org/" in re.sub(r"\s+", "", text)
        and "github.com/aconsciousfractal/FCIG-RC-Invariant" in re.sub(r"\s+", "", text),
        {}, rows,
    )
    check("pdf_no_attachments", not list(reader.attachments), {"attachments": len(reader.attachments)}, rows)
    object_findings, uris = pdf_object_findings(reader)
    check("pdf_recursive_action_audit", not object_findings, {"findings": object_findings, "uris": len(uris)}, rows)
    raw = PDF.read_bytes()
    active = [
        token.decode("ascii")
        for token in (b"/JavaScript", b"/JS", b"/Launch", b"/EmbeddedFile", b"/AA")
        if token in raw
    ]
    check("pdf_no_active_content", not active, {"markers": active}, rows)
    raw_text = raw.decode("latin-1", errors="ignore")
    raw_privacy = [
        name for name in ("windows_absolute_path", "windows_unc_path", "unix_private_path", "unix_workspace_path", "private_key", "github_token", "api_token")
        if PRIVATE_PATTERNS[name].search(raw_text)
    ]
    check("pdf_raw_privacy", not raw_privacy, {"findings": raw_privacy}, rows)


def evidence_surface(rows: list[dict[str, object]]) -> None:
    atlas = load_json(
        ROOT / "evidence" / "finite_assurance" / "results" / "certified_terminal_atlas.json",
        canonical_lf=True,
    )
    validation = load_json(
        ROOT / "evidence" / "finite_assurance" / "results" / "certified_validation.json",
        canonical_lf=True,
    )
    single_odd_prime = load_json(
        ROOT / "evidence" / "theorem" / "results" / "single_odd_prime_embedding_checks.json",
        canonical_lf=True,
    )
    assert isinstance(atlas, dict) and isinstance(validation, dict)
    assert isinstance(single_odd_prime, dict)
    check(
        "atlas_payload",
        atlas.get("scientific_payload_sha256") == ATLAS_HASH
        and len(atlas.get("scientific_payload", {}).get("rows", [])) == 30,
        {"sha256": atlas.get("scientific_payload_sha256")}, rows,
    )
    payload = validation.get("scientific_payload", {})
    check(
        "validation_payload",
        validation.get("scientific_payload_sha256") == VALIDATION_HASH
        and payload.get("mutation_panel", {}).get("mutation_count") == 29
        and payload.get("validated_row_count") == 30,
        {"sha256": validation.get("scientific_payload_sha256")}, rows,
    )
    theorem_payload = single_odd_prime.get("scientific_payload", {})
    check(
        "single_odd_prime_embedding_finite_assurance",
        single_odd_prime.get("status") == "PASS_SINGLE_ODD_PRIME_EMBEDDING_FINITE_ASSURANCE"
        and len(theorem_payload.get("cyclotomic_composition_cases", [])) == 18
        and [row.get("k") for row in theorem_payload.get("positive_binary_word_cases", [])]
        == [6, 10, 12, 18, 20, 24]
        and all(row.get("rejected") is True for row in theorem_payload.get("negative_controls", []))
        and theorem_payload.get("orbit_controls") == [
            {"k": 12, "observed_pcrs": 19, "predicted_pcrs": 19},
            {"k": 30, "observed_zero_words": 146854, "predicted_pcrs": 6069},
        ],
        {"sha256": single_odd_prime.get("scientific_payload_sha256")}, rows,
    )
    provenance = load_json(ROOT / "evidence/advances/PROVENANCE.json")
    mismatches = []
    seen = set()
    for item in provenance["inputs"]:
        rel = item["path"]
        path = ROOT / "evidence/advances" / rel
        if rel in seen or not path.resolve().is_relative_to((ROOT / "evidence/advances").resolve()):
            mismatches.append(rel)
        elif not path.is_file() or sha256_file(path) != item["packaged_sha256"]:
            mismatches.append(rel)
        seen.add(rel)
    actual = {p.relative_to(ROOT / "evidence/advances").as_posix()
              for p in (ROOT / "evidence/advances").rglob("*")
              if p.is_file() and p != ROOT / "evidence/advances/PROVENANCE.json" and "__pycache__" not in p.parts}
    # The root capsule ledger lists imported inputs; public capsule helpers are additional.
    check("advances_provenance", not mismatches and seen <= actual and len(seen) == 361,
          {"files": len(seen), "mismatches": mismatches}, rows)


    balanced = ROOT / "evidence/balanced_loads"
    provenance = load_json(balanced / "PROVENANCE.json")
    mismatches = []
    seen = set()
    for item in provenance["sources"]:
        rel = item["path"]
        path = balanced / rel
        if rel in seen or not path.resolve().is_relative_to(balanced.resolve()):
            mismatches.append(rel)
        elif (not path.is_file() or sha256_file(path) != item["sha256"]
              or path.stat().st_size != item["bytes"]):
            mismatches.append(rel)
        seen.add(rel)
    check("balanced_loads_provenance", not mismatches and len(seen) == 7,
          {"files": len(seen), "mismatches": mismatches}, rows)


def file_surface(rows: list[dict[str, object]]) -> None:
    failures: list[str] = []
    files: set[str] = set()
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        rel = relative.as_posix()
        if relative.parts[:1] == (".git",) or any(part in {"__pycache__", ".pytest_cache", "dist"} for part in relative.parts):
            continue
        if path.is_symlink():
            failures.append(f"symlink:{rel}")
            continue
        if not path.is_file():
            continue
        files.add(rel)
        if path.stat().st_size > 5 * 1024 * 1024:
            failures.append(f"oversize:{rel}")
        lower = relative.name.lower()
        if any(lower.endswith(suffix) for suffix in FORBIDDEN_ENDINGS):
            failures.append(f"generated_or_archive:{rel}")
        data = path.read_bytes()
        if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
            failures.append(f"lfs_pointer:{rel}")
        if relative.suffix.lower() == ".pdf":
            if rel != PDF_RELATIVE:
                failures.append(f"additional_pdf:{rel}")
            continue
        if relative.suffix.lower() in TEXT_SUFFIXES or relative.name in TEXT_NAMES or not relative.suffix:
            try:
                text = data.decode("utf-8", errors="strict")
            except UnicodeDecodeError:
                failures.append(f"invalid_utf8:{rel}")
                continue
            for name, pattern in PRIVATE_PATTERNS.items():
                if pattern.search(text):
                    failures.append(f"{name}:{rel}")
        elif relative.suffix.lower() in {".bits", ".levels"} and relative.parts[:2] == ("evidence", "extended"):
            if not data:
                failures.append(f"empty_scientific_binary:{rel}")
        elif relative.name.endswith(".json.gz") and relative.parts[:2] in {("evidence", "advances"), ("evidence", "balanced_loads")}:
            try:
                decoded = gzip.decompress(data).decode("utf-8", errors="strict")
                json.loads(decoded)
                for name, pattern in PRIVATE_PATTERNS.items():
                    if pattern.search(decoded):
                        failures.append(f"{name}:compressed:{rel}")
            except (ValueError, OSError, UnicodeError) as exc:
                failures.append(f"invalid_compressed_json:{rel}:{type(exc).__name__}")
        else:
            failures.append(f"unsupported_binary:{rel}")
    missing = sorted(REQUIRED - files)
    failures.extend(f"required_missing:{path}" for path in missing)
    check("file_hygiene", not failures, {"files": len(files), "failures": failures}, rows)


def validate() -> dict[str, object]:
    rows: list[dict[str, object]] = []
    file_surface(rows)
    manuscript_surface(rows)
    pdf_surface(rows)
    evidence_surface(rows)
    from build_evidence_identity import verify as verify_evidence_identity
    check("evidence_identity", verify_evidence_identity(), {}, rows)
    licence = (ROOT / "LICENSE").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    check(
        "licence",
        licence.startswith("MIT License\n")
        and "Permission is hereby granted, free of charge" in licence
        and "Copyright (c) 2026 Oleksiy Babanskyy" in licence,
        {}, rows,
    )
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    check("citation_release_metadata", "license: MIT" in cff and "date-released:" not in cff, {}, rows)
    dependency = "89a29a90962d02add6690d6411db927c4788371a"
    source_lock = (ROOT / "docs" / "SOURCE_LOCK.md").read_text(encoding="utf-8")
    check("imported_dependency_identity", dependency in readme and dependency in source_lock, {}, rows)
    check("author_identity", "Oleksiy Babanskyy" in readme, {}, rows)
    passed = all(bool(row["passed"]) for row in rows)
    return {
        "checks": rows,
        "checks_passed": sum(bool(row["passed"]) for row in rows),
        "checks_total": len(rows),
        "paper_sha256": sha256_file(PDF) if PDF.is_file() else None,
        "pass": passed,
        "schema": "fcig-cwc-package-validation-v1",
        "status": "PASS_PUBLIC_PACKAGE_VALIDATION" if passed else "FAIL_PUBLIC_PACKAGE_VALIDATION",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = validate()
    except Exception as exc:
        result = {
            "checks": [], "checks_passed": 0, "checks_total": 0, "pass": False,
            "schema": "fcig-cwc-package-validation-v1",
            "status": "FAIL_PUBLIC_PACKAGE_VALIDATION",
            "error": f"{type(exc).__name__}:{exc}",
        }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
            encoding="utf-8", newline="\n",
        )
    print(json.dumps({key: result[key] for key in ("status", "pass", "checks_passed", "checks_total")}, sort_keys=True))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
