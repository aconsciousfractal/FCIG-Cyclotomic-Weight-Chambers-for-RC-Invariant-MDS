"""Verify the v3 finite claims in a short, isolated directory, without solvers.

Every child receives the requested optimization flag explicitly. Certificate
checks use explicit exceptions; the source adaptations are in PROVENANCE.json.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--optimized',action='store_true');ap.add_argument('--output',type=Path);ap.add_argument('--logs',type=Path);args=ap.parse_args()
    optimized=args.optimized or bool(sys.flags.optimize);rows=[]
    if os.name=='nt':
        for key in ['output','logs']:
            value=getattr(args,key)
            if value is not None:
                absolute=str(value.resolve())
                if not absolute.startswith('\\\\?\\'):absolute='\\\\?\\'+absolute
                setattr(args,key,Path(absolute))
    with tempfile.TemporaryDirectory(prefix='cwc-v3-') as td:
        w=Path(td);d=w/'data';out=w/'out';out.mkdir()
        shutil.copytree(ROOT/'evidence/advances',d,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        env={**os.environ,'PYTHONUTF8':'1','PYTHONDONTWRITEBYTECODE':'1','R9_CHECK_OUTPUT':str(out/'r9')}
        env.pop('PYTHONOPTIMIZE',None)
        prefix=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])
        def run(name,path,*params):
            command=prefix+[str(d/path),*map(str,params)]
            p=subprocess.run(command,cwd=w,env=env,capture_output=True,text=True,encoding='utf-8',timeout=1800)
            if args.logs:
                args.logs.mkdir(parents=True,exist_ok=True);(args.logs/(name+'.txt')).write_text(p.stdout+p.stderr,encoding='utf-8',newline='\n')
            if p.returncode:
                print(p.stdout);print(p.stderr,file=sys.stderr);raise ValueError(name+' failed: '+str(p.returncode))
            rows.append({'check':name,'passed':True});print('PASS '+name,flush=True)
        run('independent_endpoints','verify_endpoints.py','--out',out/'endpoints.json')
        run('completion_and_surplus','verify_completion_excess.py','--out',out/'completion.json')
        for r in [2,4,5]:run('local_'+str(r),'verify_local.py',r,'--out',out/f'local_{r}.json')
        for name in ['verify_round8','verify_model8','test_mutations']:run('three_site_'+name,'r8/code/'+name+'.py','--out',out/(name+'.json'))
        run('stable_fibre','r9/code/verify_round9.py','--out',out/'stable_fibre.json')
        run('stable_small_carriers','r9/code/test_theorems.py')
        run('two_sat_truth_tables','r9/code/stable_oracle.py','--self-test')
        run('stable_mutations','r9/code/test_mutations.py')
        run('height_theorem','r10/code/test_theorem.py','--out',out/'r10')
        for name in ['verify_round11','test_regressions','test_mutations']:run('excess_'+name,'r11/code/'+name+'.py','--out',out/'r11')
        for name in ['verify_domains','test_regressions','test_mutations']:run('regions_'+name,'r12/code/'+name+'.py','--out',out/'r12')
        for name in ['verify_certificates','test_regressions','test_mutations']:run('affine_'+name,'r13/code/'+name+'.py','--out',out/'r13')
        result={'schema':'fcig-cwc-advances-verification-v1','status':'PASS','optimized':optimized,'checks':rows,
                'scope':'Packaged positive witnesses and scoped finite proofs; general theorems require the manuscript proofs.',
                'provenance_sha256':hashlib.sha256((ROOT/'evidence/advances/PROVENANCE.json').read_bytes()).hexdigest(),
                'endpoints':json.loads((out/'endpoints.json').read_text()),'completion':json.loads((out/'completion.json').read_text())}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':result['status'],'checks':len(rows),'optimized':optimized}))
    return 0
if __name__=='__main__':raise SystemExit(main())
