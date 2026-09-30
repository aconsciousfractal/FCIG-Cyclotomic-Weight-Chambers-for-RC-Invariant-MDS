"""Replay the balanced-load certificate in isolation, using standard Python."""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--optimized', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='cwc-balanced-') as temporary:
        capsule = Path(temporary) / 'balanced_loads'
        shutil.copytree(ROOT / 'evidence/balanced_loads', capsule,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        environment = {**os.environ, 'PYTHONUTF8': '1', 'PYTHONDONTWRITEBYTECODE': '1'}
        environment.pop('PYTHONOPTIMIZE', None)
        command = [sys.executable, '-X', 'utf8', '-B']
        if args.optimized or sys.flags.optimize:
            command.append('-O')
        command += [str(capsule / 'verify_case.py'), '--capsule', str(capsule)]
        completed = subprocess.run(command, cwd=capsule, env=environment,
            capture_output=True, text=True, encoding='utf-8', timeout=900)
        if completed.returncode:
            print(completed.stdout, end='')
            print(completed.stderr, end='', file=sys.stderr)
            return completed.returncode
        result = json.loads(completed.stdout)
        print(json.dumps(result, sort_keys=True))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
