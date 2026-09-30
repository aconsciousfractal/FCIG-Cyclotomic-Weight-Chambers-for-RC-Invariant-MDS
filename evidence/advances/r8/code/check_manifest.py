#!/usr/bin/env python3
"""Check byte identity and complete file coverage. Hashes do not prove mathematics."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();lines=(ROOT/'MANIFEST_SHA256.txt').read_text().splitlines();seen=set()
 for line in lines:
  sha,rel=line.split('  ',1)
  if rel in seen or not (ROOT/rel).is_file() or '..' in Path(rel).parts:raise ValueError('MANIFEST_PATH')
  seen.add(rel)
  if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha:raise ValueError('MANIFEST_HASH: '+rel)
 # Include the parent's manifest as data; exclude only the manifest being checked.
 actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=ROOT/'MANIFEST_SHA256.txt'}
 if seen!=actual:raise ValueError('MANIFEST_COVERAGE')
 r={'status':'PASS','files':len(seen),'purpose':'BYTE_IDENTITY_NOT_TRUTH'}
 if a.out:a.out.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
