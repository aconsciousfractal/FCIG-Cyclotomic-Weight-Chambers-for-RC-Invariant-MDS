#!/usr/bin/env python3
"""Replay all Round 9 checks; writes receipts outside the frozen input/log set."""
from pathlib import Path
import os,sys,json,subprocess,time,argparse,hashlib
ROOT=Path(__file__).resolve().parents[1]

def manifest_check():
    p=ROOT/'MANIFEST_SHA256.txt'
    if not p.exists():return {'status':'NOT_YET_BUILT'}
    count=0
    for line in p.read_text().splitlines():
        digest,name=line.split('  ',1);q=ROOT/name
        if not q.is_file() or hashlib.sha256(q.read_bytes()).hexdigest()!=digest:raise ValueError('MANIFEST:'+name)
        count+=1
    return {'status':'PASS','files':count}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'run_output');ap.add_argument('--skip-walk',action='store_true');args=ap.parse_args();args.out=args.out.resolve();args.out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['R9_CHECK_OUTPUT']=str(args.out)
    commands=[['verify_round9.py','--out',str(args.out/'full_verification.json')],['test_theorems.py'],['stable_oracle.py','--self-test'],['test_mutations.py'],['audit_master.py'],['audit_half_depth.py']]
    if not args.skip_walk:commands.append(['audit_walk.py'])
    receipt={'status':'RUNNING','manifest':manifest_check(),'checks':[],'python':sys.version,'scope':'Round 9 only; no external optimization solver'};start=time.time()
    for cmd in commands:
        t=time.time();name=cmd[0][:-3];command=[sys.executable,'-B',str(ROOT/'code'/cmd[0])]+cmd[1:]
        with open(args.out/(name+'.stdout.log'),'w') as f:
            try:result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=420);status='PASS' if result.returncode==0 else 'FAIL';code=result.returncode
            except subprocess.TimeoutExpired:status='TIMEOUT';code=None
        row={'name':name,'status':status,'returncode':code,'seconds':time.time()-t};receipt['checks'].append(row);print(row,flush=True)
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2))
        if status!='PASS':receipt['status']=status;break
    else:receipt['status']='PASS'
    receipt['seconds']=time.time()-start;receipt['manifest_after']=manifest_check();(args.out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(receipt['status'],flush=True)
    if receipt['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()
