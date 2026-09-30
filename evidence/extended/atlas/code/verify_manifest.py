"""Check SHA-256 byte identities. This does not validate mathematical claims."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    manifest = ROOT / 'MANIFEST_SHA256.txt'
    seen = set()
    for number, line in enumerate(manifest.read_text(encoding='utf8').splitlines(), 1):
        if not line.strip():
            continue
        digest, rel = line.split('  ', 1)
        path = ROOT / rel
        if (len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest)
                or rel in seen or not path.resolve().is_relative_to(ROOT.resolve())):
            raise SystemExit(f'MANIFEST_FORMAT_OR_PATH line {number}')
        seen.add(rel)
        if not path.is_file():
            raise SystemExit('MISSING_FILE ' + rel)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise SystemExit('HASH_MISMATCH ' + rel)
    if not seen:
        raise SystemExit('EMPTY_MANIFEST')
    print(f'PASS: {len(seen)} listed files have matching SHA-256 identities.')


if __name__ == '__main__':
    main()
