#!/usr/bin/env python3
"""Extract an archive into a NEW directory and rerun all checks.
Outputs a receipt outside the archive to avoid a circular archive digest.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,time,zipfile

def main():
 ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('--extract-to',type=Path,required=True);ap.add_argument('--receipt',type=Path,required=True);a=ap.parse_args()
 if a.extract_to.exists():raise SystemExit('Extraction directory must not already exist')
 a.extract_to.mkdir(parents=True);t=time.time()
 with zipfile.ZipFile(a.archive) as z:
  if z.testzip() is not None:raise ValueError('ZIP_CRC')
  for name in z.namelist():
   if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError('UNSAFE_ARCHIVE_PATH')
  z.extractall(a.extract_to)
 roots=[p for p in a.extract_to.iterdir() if p.is_dir()]
 if len(roots)!=1:raise ValueError('ARCHIVE_ROOT_COUNT')
 root=roots[0];checks=[]
 files=['parent_round7/code/verify_round7.py','parent_round7/code/verify_model.py','parent_round7/code/audit_search.py','parent_round7/code/test_mutations.py','parent_round7/code/check_manifest.py','code/verify_round8.py','code/verify_model8.py','code/audit_round8.py','code/test_mutations.py','code/check_manifest.py']
 for relative in files:
  start=time.time();r=subprocess.run([sys.executable,'-B',str(root/relative)],cwd=root,capture_output=True,text=True,timeout=120)
  outpath=a.extract_to/(relative.replace('/','_')+'.log');outpath.write_text(r.stdout+r.stderr)
  parsed=None
  try:parsed=json.loads(r.stdout)
  except json.JSONDecodeError:pass
  checks.append({'script':relative,'returncode':r.returncode,'seconds':time.time()-start,'stdout_sha256':hashlib.sha256(r.stdout.encode()).hexdigest(),'output':parsed if parsed is not None else r.stdout,'stderr':r.stderr})
  print(relative,r.returncode,flush=True)
  if r.returncode:break
 receipt={'schema':'round8-archive-extract-and-recheck-v1','archive_name':a.archive.name,'archive_size_bytes':a.archive.stat().st_size,'archive_sha256':hashlib.sha256(a.archive.read_bytes()).hexdigest(),'zip_crc':'PASS','fresh_extraction':str(a.extract_to),'checks':checks,'status':'PASS' if len(checks)==len(files) and all(x['returncode']==0 for x in checks) else 'FAIL','seconds':time.time()-t,'scope':'Byte identity and the listed exact finite checkers, not formal verification or human independent peer review'}
 a.receipt.write_text(json.dumps(receipt,indent=2)+'\n')
 print(receipt['status'],receipt['archive_sha256'],flush=True)
 if receipt['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()
