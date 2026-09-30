"""Rebuild every Round 3 JSON product in a temporary directory and compare bytes.

Run from any working directory using Python 3.10+. No network or third-party
packages are needed.  --optimized runs each scientific child under python -O.
The output distributed in results/ is never overwritten.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
from core import check, locate_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--parent', type=Path, help='Optional already-extracted frozen Round 2 root.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    parent, parent_temp = locate_parent(args.parent)
    try:
        with tempfile.TemporaryDirectory(prefix='fcig_round3_replay_') as name:
            out = Path(name)
            flags = ['-O'] if args.optimized or not __debug__ else []
            for script, extra in [('produce.py', []),
                                  ('verify.py', ['--parent-replay']),
                                  ('mutations.py', []), ('probe_paths.py', [])]:
                command = [sys.executable, *flags, '-B', str(root/'code'/script),
                           '--parent', str(parent), '--out', str(out), *extra]
                print('RUN', script, 'optimized=' + str(bool(flags)), flush=True)
                subprocess.run(command, check=True)
            expected = {p.name for p in (root/'results').glob('*.json')}
            actual = {p.name for p in out.glob('*.json')}
            check(expected == actual, 'REPLAY_FILE_SET_MISMATCH: ' + repr(expected ^ actual))
            for name in sorted(expected):
                check((out/name).read_bytes() == (root/'results'/name).read_bytes(),
                      'REPLAY_BYTE_MISMATCH: ' + name)
            print('ROUND3_REPLAY_PASS', len(expected), 'JSON products byte-identical', flush=True)
    finally:
        if parent_temp:
            parent_temp.cleanup()


if __name__ == '__main__':
    main()
