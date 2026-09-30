"""Replay every integrated scientific capsule in an isolated temporary tree."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1]
def need(ok,msg):
    if not ok: raise ValueError(msg)
def clean(value):
    if isinstance(value,dict): return {k:clean(v) for k,v in value.items() if k not in {'seconds','elapsed_seconds'}}
    if isinstance(value,list): return [clean(v) for v in value]
    return value
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--optimized',action='store_true')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args();checks=[]
    with tempfile.TemporaryDirectory(prefix='fcig-cwc-extended-') as td:
        w=Path(td); e=w/'evidence/extended'
        shutil.copytree(ROOT/'evidence',w/'evidence',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if args.optimized else [])
        def run(name,path,*a):
            p=subprocess.run(cmd+[str(path),*map(str,a)],cwd=w,text=True,capture_output=True,
                             env={**os.environ,'PYTHONUTF8':'1'},check=False)
            if p.returncode:
                print(p.stdout);print(p.stderr,file=sys.stderr)
                raise ValueError(f'{name}: exit {p.returncode}')
            checks.append({'check':name,'passed':True})
            print('PASS',name,flush=True)
        for name in ('run_checks','additional_checks','check_additive_arrangement','check_two_family'):
            run('algebra_'+name,e/'algebra/code'/f'{name}.py')
        for name in ('new_results','additional_checks','three_magnitude_arrangement','two_family_and_residual'):
            old=json.loads((ROOT/'evidence/extended/algebra/results'/f'{name}.json').read_text(encoding='utf-8'))
            fresh=json.loads((e/'algebra/results'/f'{name}.json').read_text(encoding='utf-8'))
            need(clean(old)==clean(fresh),'algebra semantic drift '+name)
        checks.append({'check':'algebra_frozen_payloads','passed':True})
        run('atlas_complete',e/'atlas/code/run_all.py',*(['--optimized'] if args.optimized else []))
        run('walls_complete',e/'walls/code/run_all.py','--parent',e/'atlas',*(['--optimized'] if args.optimized else []))
        run('retiming_certificates',e/'retiming/code/verify.py','--root',e/'retiming','--output',w/'retiming.json')
        run('rounding_42',e/'retiming/code/replay_rounding.py')
        run('four_patches',e/'retiming/code/construct_131.py','--out',w/'constructed.bits')
        run('independent_geometry',e/'checks/verify_independent.py',e/'retiming','--output',w/'geometry.json')
        run('independent_cyclotomic',e/'checks/verify_cyclotomic.py',e/'algebra','--output',w/'cyclotomic.json')
        for key,file in [('retiming','retiming.json'),('geometry','geometry.json'),('cyclotomic','cyclotomic.json')]:
            need(json.loads((w/file).read_text())['status']=='PASS',key+' status')
        need((w/'constructed.bits').read_bytes()==(ROOT/'evidence/extended/retiming/certificates/escape_131.bits').read_bytes(),'final bit identity')
        result={'schema':'fcig-cwc-extended-verification-v1','status':'PASS','optimized':args.optimized,
                'checks':checks,'final_bits_sha256':hashlib.sha256((w/'constructed.bits').read_bytes()).hexdigest(),
                'scope':'Exact finite claims; written uniform proofs are separate.'}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,sort_keys=True))
    return 0
if __name__=='__main__':raise SystemExit(main())
