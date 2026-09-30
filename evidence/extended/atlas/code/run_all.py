"""Regenerate in an isolated temporary tree and compare exact product bytes.
Usage: python -B code/run_all.py [--optimized]
The release snapshot is not overwritten. Standard library only.
"""
from pathlib import Path
import sys,subprocess,tempfile,shutil,hashlib,json
ROOT=Path(__file__).resolve().parents[1]

def main():
    optimized='--optimized' in sys.argv[1:]
    if any(x!='--optimized' for x in sys.argv[1:]):raise SystemExit('Only --optimized is accepted')
    cmd=[sys.executable,'-B']+(['-O'] if optimized else [])
    with tempfile.TemporaryDirectory(prefix='fcig_r2_') as td:
        work=Path(td)
        shutil.copytree(ROOT/'code',work/'code',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(ROOT/'sources',work/'sources')
        for script in ('atlas_engine.py','verify_atlas.py','derived_checks.py','mutation_checks.py'):
            p=subprocess.run(cmd+[str(work/'code'/script)],cwd=work,text=True,capture_output=True)
            print(p.stdout,end='',flush=True)
            if p.returncode:
                print(p.stderr,file=sys.stderr);raise SystemExit(p.returncode)
        products=sorted((work/'results').rglob('*'));compared=0;h=hashlib.sha256()
        for fresh in products:
            if not fresh.is_file():continue
            rel=fresh.relative_to(work/'results');old=ROOT/'results'/rel
            if not old.exists() or fresh.read_bytes()!=old.read_bytes():raise SystemExit('BYTE_MISMATCH '+str(rel))
            compared+=1;h.update(rel.as_posix().encode()+b'\0'+fresh.read_bytes())
        print(json.dumps({'status':'PASS','optimized':optimized,'product_files_compared':compared,'exact_product_sha256':h.hexdigest()},sort_keys=True))
if __name__=='__main__':main()
