"""Verify SHA-256 identities and the complete distributed payload inventory."""
from __future__ import annotations
import hashlib
from pathlib import Path
import sys


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = root / 'MANIFEST_SHA256.txt'
    expected: set[str] = set()
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, relative = line.split('  ', 1)
        path = Path(relative)
        if path.is_absolute() or '..' in path.parts or relative in expected:
            raise ValueError('UNSAFE_OR_DUPLICATE_MANIFEST_PATH: ' + relative)
        expected.add(relative)
        target = root / path
        if not target.is_file():
            raise ValueError('MISSING_FILE: ' + relative)
        actual = hashlib.sha256(target.read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError('HASH_MISMATCH: ' + relative)
    actual_files = {p.relative_to(root).as_posix() for p in root.rglob('*')
                    if p.is_file() and '__pycache__' not in p.parts
                    and p != manifest}
    if expected != actual_files:
        raise ValueError('INVENTORY_MISMATCH: ' + repr(sorted(expected ^ actual_files)))
    print('MANIFEST_PASS', len(expected), 'files; identity is not mathematical validity')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
