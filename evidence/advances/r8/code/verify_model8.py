#!/usr/bin/env python3
"""Independent reconstruction of the full model from literal DNA strings.
Imports only archived standard-library checkers, not the numeric compiler.
"""
from pathlib import Path
import sys,json,time,argparse
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'parent_round7/code'))
import verify_round7 as p
import verify_model as prior

def verify(model):
 base=p.read_words(ROOT/'inputs/M107.words.txt');stable=p.pcr_data(base)
 p.require(model['schema']=='complete-bounded-RC-chart-strengthened-v1','MODEL_SCHEMA')
 p.require(model['T']==107 and model['q']==4 and model['k']==8 and model['reference_letter']=='T' and model['base_sha256']==p.digest(base),'MODEL_SCOPE')
 cycles=model['strengthening']['cycles'];vs,rows,fixed=prior.expected(base,stable,107,cycles)
 p.require(model['variables']==vs and model['rows']==json.loads(json.dumps(rows)),'GLOBAL_MODEL_ROWS')
 p.require(model['strengthening']['fixed_unary_edges']==fixed and model['strengthening']['globally_valid'] is True,'GLOBAL_MODEL_UNARY')
 p.require(model['stable_pairs']==[list(x) for x in stable],'GLOBAL_MODEL_STABLE')
 p.require(model['round8_metadata']=={'site_cap':None,'gamma_bits_fixed':0,'local_gcd_rows_copied':False,'global_target_solved':False},'GLOBAL_MODEL_BOUNDARY')
 q=p.load(ROOT/'certificates/triple_refutations.json.gz');all_cycles=q['cycles'].copy();seen={tuple(sorted(Counter(w).items())) for w in all_cycles}
 for path in sorted((ROOT/'discovery/results').glob('adaptive*/result.json')):
  run=p.load(path)
  for cy in run['cycles']:
   key=tuple(sorted(Counter(cy).items()))
   if key not in seen:seen.add(key);all_cycles.append(cy)
 p.require(cycles==all_cycles,'GLOBAL_CYCLE_SOURCE')
 p.require(len(vs)==32896 and len(rows)==34482 and len(cycles)==1330 and len(fixed)==256,'GLOBAL_MODEL_COUNTS')
 control=prior.verify_control(base,stable,p.load(ROOT/'parent_round7/inputs/chart_M107_T108.json.gz'),cycles)
 return {'status':'PASS','variables':len(vs),'rows':len(rows),'raw_walk_cuts':len(cycles),'unary_tightenings':len(fixed),'no_site_cap':True,'all_128_stable_bits_free':True,'target_solved':False,'positive_control':control}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();ts=time.time();r=verify(p.load(ROOT/'models/complete_T107_round8.json.gz'));r['seconds']=time.time()-ts
 if a.out:a.out.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
