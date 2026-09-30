"""Regression checks for byte identity, active sources, and obvious path hygiene."""
from pathlib import Path
import gzip
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build_evidence_identity import outputs, verify
from manuscript_sources import active_tex_files
import validate_repository
import check_git_history


def forbidden_path() -> str:
    # Assemble the fixture so the distributed test source itself is clean.
    return "/mnt/" + "data/hidden/work.json"


def test_evidence_identity_detects_mutation_and_has_no_self_reference(tmp_path):
    (tmp_path / "evidence").mkdir()
    (tmp_path / "paper").mkdir()
    (tmp_path / "results").mkdir()
    payload = tmp_path / "evidence/input.json"
    payload.write_bytes(b'{"witness":107}\n')
    receipt, tex = outputs(tmp_path)
    (tmp_path / "results/evidence_identity.json").write_bytes(receipt)
    (tmp_path / "paper/evidence_identity.tex").write_bytes(tex)
    assert verify(tmp_path)
    (tmp_path / "paper/main.tex").write_text("Changed manuscript")
    assert outputs(tmp_path) == (receipt, tex)
    payload.write_bytes(b'{"witness":106}\n')
    assert not verify(tmp_path)
    assert outputs(tmp_path) != (receipt, tex)
    value = json.loads(receipt)
    assert value["files"][0]["path"] == "evidence/input.json"
    assert len(value["sha256"]) == 64


def test_active_source_tree_excludes_unused_files(tmp_path):
    (tmp_path / "main.tex").write_text("% \\input{unused}\n\\input{active}\n")
    (tmp_path / "active.tex").write_text("A theorem")
    (tmp_path / "unused.tex").write_text("Old theorem")
    assert [p.name for p in active_tex_files(tmp_path)] == ["main.tex", "active.tex"]


def test_resealed_private_path_is_still_rejected(tmp_path, monkeypatch):
    (tmp_path / "evidence").mkdir()
    path = tmp_path / "evidence/rogue.json"
    path.write_text(json.dumps({"base": forbidden_path()}))
    # A freshly valid digest/manifest cannot make the content public.
    from package_common import sha256_file
    (tmp_path / "MANIFEST_SHA256.txt").write_text(sha256_file(path) + "  evidence/rogue.json\n")
    import package_common
    import check_manifest
    monkeypatch.setattr(package_common, "ROOT", tmp_path)
    monkeypatch.setattr(package_common, "MANIFEST", tmp_path / "MANIFEST_SHA256.txt")
    monkeypatch.setattr(package_common, "ATTESTATION", tmp_path / "RELEASE_ATTESTATION.json")
    monkeypatch.setattr(check_manifest, "ROOT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["check_manifest.py", "--closed-tree"])
    assert check_manifest.main() == 0
    monkeypatch.setattr(validate_repository, "ROOT", tmp_path)
    monkeypatch.setattr(validate_repository, "REQUIRED", set())
    rows = []
    validate_repository.file_surface(rows)
    assert "unix_workspace_path:evidence/rogue.json" in rows[0]["detail"]["failures"]


def test_compressed_private_path_is_rejected(tmp_path, monkeypatch):
    target = tmp_path / "evidence/advances/rogue.json.gz"
    target.parent.mkdir(parents=True)
    target.write_bytes(gzip.compress(json.dumps({"base": forbidden_path()}).encode()))
    monkeypatch.setattr(validate_repository, "ROOT", tmp_path)
    monkeypatch.setattr(validate_repository, "REQUIRED", set())
    rows = []
    validate_repository.file_surface(rows)
    assert "unix_workspace_path:compressed:evidence/advances/rogue.json.gz" in rows[0]["detail"]["failures"]


def test_reachable_history_rejects_removed_path(tmp_path, monkeypatch, capsys):
    def git(*args):
        result = subprocess.run(["git", *args], cwd=tmp_path, capture_output=True, check=True)
        return result.stdout
    git("init", "-q")
    (tmp_path / "payload.json").write_text(json.dumps({"base": forbidden_path()}))
    git("add", "payload.json")
    git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "historical input")
    (tmp_path / "payload.json").write_text('{}')
    git("add", "payload.json")
    git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "remove path")
    monkeypatch.setattr(check_git_history, "ROOT", tmp_path)
    import pytest
    with pytest.raises(SystemExit) as caught:
        check_git_history.main()
    assert caught.value.code == 1
    assert "workspace path in historical blob payload.json" in capsys.readouterr().out


def test_scanner_sources_and_historical_hashes_are_clean():
    pattern = validate_repository.PRIVATE_PATTERNS["unix_workspace_path"]
    assert not pattern.search(Path(__file__).read_text())
    assert not pattern.search((ROOT / "scripts/validate_repository.py").read_text())
    assert not pattern.search((ROOT / "scripts/check_git_history.py").read_text())
    historical = json.dumps({"original_hash": "b72f90514a267f8b5e79743ce84a95672b872837af1c663941ce5a79ebf9226f"})
    assert not pattern.search(historical)
