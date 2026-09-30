#!/usr/bin/env python3
"""Independent string-based reconstruction of every global model row.
Does not import the compiler, discovery code, NumPy, or a solver.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from collections import Counter
import verify_round7 as v
ROOT=Path(__file__).resolve().parents[1]

def expected(base,stable,T,cycles):
 variables=[{'name':f'z{i:05d}','lower':-2,'upper':2,'kind':'integer'} for i in range(3*len(v.SITES))]
 variables += [{'name':f'b{i:03d}','lower':0,'upper':1,'kind':'binary'} for i in range(128)]
 variables += [{'name':f'h{i:05d}','lower':-T,'upper':2*T-1,'kind':'integer'} for i in range(len(v.SITES))]
 # Edge cochain coefficients, in original units (not multiplied by T).
 ec={w:Counter() for w in v.E}
 for i,u in enumerate(v.SITES):
  for a in range(3):
   for w,k in v.generator(u,v.ALPHABET[a]).items():ec[w][f'z{3*i+a:05d}']+=k
 for i,(a,b) in enumerate(stable):ec[a][f'b{i:03d}']-=1;ec[b][f'b{i:03d}']+=1
 rows=[];tight=[]
 for w in v.E:
  if w>v.rc(w):continue
  co=Counter({n:T*k for n,k in ec[w].items()});constant=0
  for u,sgn in [(w[1:],1),(w[:-1],-1)]:
   if u<v.rc(u):co[f'h{v.SI[u]:05d}']+=sgn
   else:co[f'h{v.SI[v.rc(u)]:05d}']-=sgn;constant+=sgn*(T-1)
  e=v.EI[w] if hasattr(v,'EI') else sum(v.ALPHABET.index(a)*4**(7-i) for i,a in enumerate(w))
  row={'name':f'e{e:05d}','coefficients':sorted((n,k) for n,k in co.items() if k),'rhs':1-T*int(w in base)-constant};rows.append(row)
  if w==v.rc(w):tight.append((e,row))
 for e,r in tight:
  co=dict(r['coefficients']);bits=[(n,k) for n,k in co.items() if n.startswith('b')];v.require(len(bits)==1,'UNARY_UNIQUE_BIT');name,k=bits[0];b=r['rhs'];a0=-(-b//2);a1=-(-(b-k)//2)
  other={n:c for n,c in co.items() if n!=name};v.require(all(c%2==0 for c in other.values()),'UNARY_EVEN_COEFFICIENTS')
  norm={n:c//2 for n,c in other.items()};norm[name]=a0-a1
  rows.append({'name':'tight_'+r['name'],'coefficients':sorted((n,c) for n,c in norm.items() if c),'rhs':a0})
 for i,cy in enumerate(cycles):
  hits=v.walk(cy,base);co=Counter()
  for w in cy:
   for n,k in ec[w].items():co[n]+=k
  rows.append({'name':f'cycle_{i:04d}','coefficients':sorted((n,k) for n,k in co.items() if k),'rhs':(len(cy)+T-1)//T-hits})
 return variables,rows,[e for e,row in tight]

def check_model(base,stable,j):
 v.require(j['schema']=='complete-bounded-RC-chart-strengthened-v1','MODEL_SCHEMA');v.require(j['T']==107 and j['q']==4 and j['k']==8 and j['reference_letter']=='T' and j['base_sha256']==v.digest(base),'MODEL_SCOPE')
 cycles=j['strengthening']['cycles'];vs,rs,fixed=expected(base,stable,107,cycles)
 v.require(j['variables']==vs,'MODEL_VARIABLES');v.require(j['rows']==json.loads(json.dumps(rs)),'MODEL_ROWS');v.require(j['strengthening']['fixed_unary_edges']==fixed,'MODEL_FIXED_CENSUS')
 v.require(j['stable_pairs']==[list(p) for p in stable],'MODEL_STABLE_DICTIONARY')
 v.require(len(vs)==32896 and len(rs)==33345 and len(fixed)==256 and len(cycles)==193,'MODEL_COUNTS')
 return {'status':'PASS','variables':len(vs),'rows':len(rs),'unary_tightenings':len(fixed),'cycle_cuts':len(cycles),'no_site_restriction':True}

def verify_control(base,stable,chart,cycles):
 v.require(chart['T']==108 and chart['base_sha256']==v.digest(base) and chart['target_sha256']==v.digest(base),'CONTROL_BASE_PERIOD')
 T=108;vs,rs,fixed=expected(base,stable,T,cycles)
 vals=chart['z']+chart['gamma']+chart['H_representatives'];v.require(len(vals)==len(vs) and all(v.integer(a) for a in vals),'CONTROL_VALUES')
 sol={r['name']:a for r,a in zip(vs,vals)}
 v.require(all(r['lower']<=sol[r['name']]<=r['upper'] for r in vs),'CONTROL_BOUNDS')
 v.require(all(sum(k*sol[n] for n,k in r['coefficients'])>=r['rhs'] for r in rs),'CONTROL_ROW_VIOLATION')
 return {'period':108,'L_upper':107,'all_strengthened_rows_satisfied':len(rs),'target_106_solution':False}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();t=time.time();base=v.read_words(ROOT/'inputs/M107.words.txt');stable=v.pcr_data(base);j=v.load(ROOT/'models/complete_T107_strengthened.json.gz');r=check_model(base,stable,j)
 r['positive_control']=verify_control(base,stable,v.load(ROOT/'inputs/chart_M107_T108.json.gz'),j['strengthening']['cycles']);r['seconds']=time.time()-t
 if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
