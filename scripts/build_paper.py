#!/usr/bin/env python3
"""Build the article twice with Tectonic and verify the committed PDF."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from pypdf import PdfReader

from package_common import PDF, ROOT, sha256_file
from manuscript_sources import active_tex_files
from build_evidence_identity import verify as verify_evidence_identity


PAPER = ROOT / "paper"
OUTPUT_NAME = "Cyclotomic_Weight_Chambers_for_RC_Invariant_MDS.pdf"
RECEIPT = ROOT / "results" / "paper_build.json"
SOURCE_DATE_EPOCH = "1786924800"
EXPECTED_TOOLCHAIN = "Tectonic 0.16.9"
CLAIM_LABELS = (
    "lem:stable-pcr-parity", "lem:signed-reconstruction", "lem:one-fixed-arc",
    "thm:headline-separation",
    "thm:component-field", "cor:dna-field-gate", "cor:uniform-radical",
    "thm:twice-prime-structural", "thm:twice-prime-embedding",
    "cor:rational-lower-energy",
    "thm:terminal-covector", "thm:dna-certificate", "thm:rational-twice-prime",
    "thm:galois-certificate", "thm:membership-complexity", "prop:finite-atlas",
    'ext:thm:two',
    'ext:thm:farey',
    'ext:thm:freiman',
    'ext:thm:count',
    'ext:thm:mid',
    'ext:thm:algebraic',
    'ext:thm:paths',
    'thm:regular-wall',
    'thm:mixed-wall',
    'thm:complete-atlas',
    'prop:quadratic-embedding',
    'prop:family-optima',
    'thm:unit-retiming',
    'thm:rc-retiming',
    'thm:rc-rounding',
    'thm:two-rc-sectors',
    'prop:family-components',
    'prop:component-optima',
    'prop:construction-131',
    'cor:global-upper-131',
    'prop:sharp-parity', 'thm:construction-107', 'thm:bounded-chart',
    'thm:stable-completion', 'prop:stable-fibre', 'cor:excess-one-oracle',
    'prop:boolean-ranks', 'cor:surplus-two-arm', 'thm:height-formulation',
    'prop:8232', 'prop:affine-pruning', 'prop:local-exclusions',
    'lem:walk-budget', 'prop:balanced-normal-form', 'thm:finite-pool-closure',
    'cor:mds-pool', 'prop:core-pool-strictness',
)
BAD_LOG = (
    r"undefined references", r"Citation .* undefined", r"Overfull \\[hv]box",
)


def source_files() -> tuple[Path, ...]:
    return tuple(sorted(
        (*active_tex_files(PAPER), PAPER / "refs.bib"),
        key=lambda path: path.relative_to(PAPER).as_posix(),
    ))


def source_digest() -> str:
    digest = hashlib.sha256()
    for path in source_files():
        digest.update(path.relative_to(PAPER).as_posix().encode("utf-8") + b"\0")
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest()


def parse_labels(auxiliary: Path) -> dict[str, dict[str, object]]:
    pattern = re.compile(r"^\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{([^}]*)\}")
    labels: dict[str, dict[str, object]] = {}
    for line in auxiliary.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match and match.group(1) in CLAIM_LABELS:
            label, number, page = match.groups()
            labels[label] = {"number": number, "page": int(page)}
    if set(labels) != set(CLAIM_LABELS):
        raise RuntimeError(f"claim_label_set_drift:{sorted(set(CLAIM_LABELS)-set(labels))}")
    return dict(sorted(labels.items()))


def one_build(parent: Path, label: str, tectonic: str) -> tuple[Path, dict[str, dict[str, object]], str]:
    work = parent / label
    work.mkdir(parents=True)
    for source in source_files():
        target = work / source.relative_to(PAPER)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    out = work / "out"
    out.mkdir()
    env = os.environ.copy()
    env["SOURCE_DATE_EPOCH"] = SOURCE_DATE_EPOCH
    command = [tectonic, "-X", "compile", "main.tex", "--outdir", str(out), "--keep-logs", "--keep-intermediates"]
    completed = subprocess.run(command, cwd=work, env=env, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"tectonic_failed:{completed.stdout}\n{completed.stderr}")
    log = (out / "main.log").read_text(encoding="utf-8", errors="replace")
    problems = [pattern for pattern in BAD_LOG if re.search(pattern, log, flags=re.IGNORECASE)]
    if problems:
        raise RuntimeError(f"latex_warning_gate:{problems}\\n"+log)
    pdf = out / "main.pdf"
    if not pdf.is_file():
        raise RuntimeError("paper_pdf_missing")
    labels = parse_labels(out / "main.aux")
    version = subprocess.run([tectonic, "--version"], text=True, capture_output=True, check=True).stdout.splitlines()[0]
    return pdf, labels, version.strip()


def build_pair() -> tuple[bytes, dict[str, dict[str, object]], str]:
    tectonic = os.environ.get("TECTONIC") or shutil.which("tectonic")
    if not tectonic:
        raise RuntimeError("tectonic_0_16_9_required")
    with tempfile.TemporaryDirectory(prefix="fcig-cwc-paper-") as temporary:
        parent = Path(temporary)
        first, labels1, version = one_build(parent, "run1", tectonic)
        second, labels2, version2 = one_build(parent, "run2", tectonic)
        first_bytes = first.read_bytes()
        if first_bytes != second.read_bytes() or labels1 != labels2:
            raise RuntimeError("paper_build_not_deterministic")
        if version != version2 or version != EXPECTED_TOOLCHAIN:
            raise RuntimeError(
                f"paper_toolchain_identity:{version!r}:{version2!r}:expected={EXPECTED_TOOLCHAIN!r}"
            )
    return first_bytes, labels1, version


def receipt(pdf_bytes: bytes, labels: dict[str, dict[str, object]], version: str) -> dict[str, object]:
    reader = PdfReader(io.BytesIO(pdf_bytes))
    return {
        "byte_identical_runs": 2,
        "bytes": len(pdf_bytes),
        "claim_labels": labels,
        "output": f"paper/{OUTPUT_NAME}",
        "pages": len(reader.pages),
        "schema": "fcig-cwc-paper-build-v1",
        "sha256": hashlib.sha256(pdf_bytes).hexdigest(),
        "source_date_epoch": SOURCE_DATE_EPOCH,
        "source_tree_sha256": source_digest(),
        "status": "PASS_DETERMINISTIC_PAPER_BUILD",
        "toolchain": version,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-rebuild", action="store_true")
    args = parser.parse_args()
    if not verify_evidence_identity():
        raise RuntimeError("evidence_identity_drift; run scripts/build_evidence_identity.py")
    fresh, labels, version = build_pair()
    value = receipt(fresh, labels, version)
    if args.verify_rebuild:
        if not PDF.is_file() or fresh != PDF.read_bytes():
            raise RuntimeError("committed_pdf_byte_drift")
        if not RECEIPT.is_file():
            raise RuntimeError("paper_build_receipt_missing")
        receipt_bytes = RECEIPT.read_bytes()
        try:
            recorded = json.loads(receipt_bytes)
        except json.JSONDecodeError as exc:
            raise RuntimeError("paper_build_receipt_json") from exc
        expected_receipt_bytes = (
            json.dumps(recorded, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        if receipt_bytes != expected_receipt_bytes:
            raise RuntimeError("paper_build_receipt_serialization")
        if recorded != value:
            raise RuntimeError("paper_build_receipt_drift")
        print(json.dumps({
            "pages": value["pages"],
            "sha256": value["sha256"],
            "status": "PASS_PAPER_BYTE_EXACT_REBUILD_COMPARISON",
            "toolchain": version,
        }, sort_keys=True))
        return 0
    PDF.write_bytes(fresh)
    RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": value["status"], "sha256": value["sha256"], "pages": value["pages"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
