#!/usr/bin/env python3
"""Byte-identity check, not a mathematical proof."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
expected={}
for line in (ROOT/'MANIFEST_SHA256.txt').read_text(encoding='utf8').splitlines():
 h,p=line.split('  ',1)
 if p in expected or Path(p).is_absolute() or '..' in Path(p).parts:raise ValueError('MANIFEST_PATH')
 expected[p]=h
actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
if actual!=set(expected):raise ValueError('MANIFEST_FILE_SET')
for p,h in expected.items():
 if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise ValueError('MANIFEST_HASH:'+p)
print(json.dumps({'status':'PASS','files':len(expected),'purpose':'BYTE_IDENTITY_NOT_TRUTH'},indent=2))
