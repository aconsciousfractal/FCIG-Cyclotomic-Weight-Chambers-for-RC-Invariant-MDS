#!/usr/bin/env python3
"""Exact audit of discovery candidates, separate from support refutation proofs.
No NumPy, scipy, numba, solver, or discovery code is imported.
"""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,struct,time
import verify_round8 as v
p=v.p;ROOT=Path(__file__).resolve().parents[1]

def checked_delta(base,stable,sites,z,gamma):
 v.need(len(z)==3*len(sites) and all(v.integer(a) for a in z),'AUDIT_Z_DOMAIN')
 v.need(len(gamma)==128 and all(v.integer(a) and a in (0,1) for a in gamma),'AUDIT_GAMMA_DOMAIN')
 d=Counter()
 for i,u in enumerate(sites):
  v.need(u in p.SI,'AUDIT_SITE_DOMAIN')
  for j,a in enumerate('ACG'):
   if z[3*i+j]:
    for e,k in p.generator(u,a).items():d[e]+=z[3*i+j]*k
 for b,(old,new) in zip(gamma,stable):
  d[old]-=b;d[new]+=b
 d={e:k for e,k in d.items() if k}
 v.need(all(d.get(p.rc(e),0)==k for e,k in d.items()),'AUDIT_SIGNED_RC')
 totals=Counter()
 for e,k in d.items():totals[min(e[t:]+e[:t] for t in range(8))]+=k
 v.need(not any(totals.values()),'AUDIT_SIGNED_PCR')
 return d

def prepare_cycles(cycles,base):
 out=[]
 for cy in cycles:out.append((Counter(cy),p.walk(cy,base),len(cy)))
 return out

def check_rejection(delta,prepared,index,recorded_deficit=None):
 v.need(v.integer(index) and 0<=index<len(prepared),'AUDIT_CYCLE_INDEX')
 for counts,hits,length in prepared[:index]:
  value=hits+sum(k*counts.get(e,0) for e,k in delta.items())
  v.need(107*value>=length,'AUDIT_PREVIOUS_WALK_VIOLATION')
 counts,hits,length=prepared[index];value=hits+sum(k*counts.get(e,0) for e,k in delta.items());gap=length-107*value
 v.need(gap>0,'AUDIT_NONPOSITIVE_OBSTRUCTION')
 if recorded_deficit is not None:v.need(v.integer(recorded_deficit) and recorded_deficit==gap,'AUDIT_DEFICIT_MISMATCH')
 return gap

def fingerprint(delta):return hashlib.sha256(json.dumps(sorted(delta.items()),separators=(',',':')).encode()).hexdigest()

def audit_adaptive(base,stable,q):
 v.need(q['schema']=='adaptive-all-sites-cegar-v1' and q['T']==107 and q['base_sha256']==p.digest(base),'ADAPTIVE_SCOPE')
 v.need(q['all_sites_available']==8192 and q['site_cap'] is None and q['z_bounds']==[-2,2] and q['stable_bits']==128,'ADAPTIVE_DOMAINS')
 prepared=prepare_cycles(q['cycles'],base);initial=q['initial_cycles'];v.need(v.integer(initial) and initial>0,'ADAPTIVE_INITIAL')
 records=[];keys=set();gammas=Counter();sizes=Counter()
 for i,r in enumerate(q['audit']):
  v.need(r['iteration']==i,'ADAPTIVE_SEQUENCE')
  if not r['incumbent']:
   v.need(i==len(q['audit'])-1,'ADAPTIVE_NO_INCUMBENT_NOT_LAST');continue
  z=[0]*24576;gamma=[0]*128;seen=set()
  for j,k in r['operations']:
   v.need(v.integer(j) and 0<=j<24704 and j not in seen and v.integer(k),'ADAPTIVE_OPERATION_DOMAIN');seen.add(j)
   if j<24576:v.need(-2<=k<=2 and k!=0,'ADAPTIVE_Z_RANGE');z[j]=k
   else:v.need(k==1,'ADAPTIVE_BIT_VALUE');gamma[j-24576]=k
  size=sum(any(z[3*j:3*j+3]) for j in range(8192));bits=sum(gamma)
  v.need(size==r['actual_sites'] and bits==r['stable_changes'] and bits>=q['stable_min'],'ADAPTIVE_ACTUAL_ACTIVITY')
  v.need(sum(abs(a) for a in z)==r['l1_z'],'ADAPTIVE_L1')
  d=checked_delta(base,stable,p.SITES,z,gamma)
  vals=[int(e in base)+d.get(e,0) for e in p.E]
  v.need(min(vals)==r['x_min'] and max(vals)==r['x_max'],'ADAPTIVE_COCHAIN_RANGE')
  digest=hashlib.sha256(struct.pack('<65536q',*vals)).hexdigest();v.need(digest==r['cochain_sha256_le_i64'],'ADAPTIVE_NUMERIC_DIGEST')
  v.need(r['clock_status']==0,'ADAPTIVE_OUTCOME');idx=r['cycle_index'];v.need(idx==initial+len(records),'ADAPTIVE_CUT_SEQUENCE')
  gap=check_rejection(d,prepared,idx,r['deficit']);key=fingerprint(d);v.need(key not in keys,'ADAPTIVE_DUPLICATE');keys.add(key)
  records.append({'sha256_delta':key,'actual_sites':size,'stable_changes':bits,'positive_excess':gap,'walk_length':prepared[idx][2]});gammas[bits]+=1;sizes[size]+=1
 v.need(len(prepared)==initial+len(records),'ADAPTIVE_CUT_COUNT')
 v.need(q['status'] in ('WALL_BUDGET','ITERATION_LIMIT','NO_INCUMBENT','SOLVER_UNSAT_UNCERTIFIED'),'ADAPTIVE_TERMINAL_STATUS')
 return {'status':'PASS','search_status':q['status'],'candidates':len(records),'site_histogram':dict(sizes),'gamma_histogram':dict(gammas),'all_sites_available':8192,'stable_min':q['stable_min'],'search_seconds':q['seconds'],'candidate_records':records,'no_106_witness':True}

def audit_triples(base,stable,q):
 v.need(q['schema']=='round8-triple-screen-v1' and q['T']==107 and q['base_sha256']==p.digest(base),'TRIPLE_AUDIT_SCOPE')
 v.need(q['all_gamma_free']==128 and q['site_coefficients']=='arbitrary integers','TRIPLE_AUDIT_DOMAINS')
 prep=prepare_cycles(q['cycles'],base);initial=q['initial_cycles'];seen=set();records=[];gammas=Counter();sizes=Counter()
 for i,r in enumerate(q['candidate_audit']):
  sites=r['sites'];v.need(len(sites)==3 and sites==sorted(set(sites)),'TRIPLE_AUDIT_SITES');d=checked_delta(base,stable,sites,r['z'],r['gamma'])
  actual=sum(any(r['z'][3*j:3*j+3]) for j in range(3));v.need(actual==r['actual_sites'],'TRIPLE_AUDIT_ACTIVITY')
  v.need(r['status']==0 and r['walk']==q['cycles'][initial+i],'TRIPLE_AUDIT_WALK_JOIN')
  gap=check_rejection(d,prep,initial+i,r['deficit']);key=fingerprint(d);seen.add(key)
  vals=[int(e in base)+k for e,k in d.items()]
  # Both unchanged selected and unchanged unselected edges exist: support is tiny.
  v.need(len(d)<8230 and min([0,1]+vals)==r['x_min'] and max([0,1]+vals)==r['x_max'],'TRIPLE_AUDIT_RANGE')
  bits=sum(r['gamma']);gammas[bits]+=1;sizes[actual]+=1
  records.append({'sha256_delta':key,'sites':sites,'actual_sites':actual,'stable_changes':bits,'positive_excess':gap,'walk_length':len(r['walk'])})
 v.need(len(prep)==initial+len(records),'TRIPLE_AUDIT_CUT_COUNT')
 return {'status':'PASS','candidate_records':records,'candidates':len(records),'distinct_cochains':len(seen),'gamma_histogram':dict(gammas),'site_histogram':dict(sizes),'search_seconds':q['seconds'],'discovery_status':q['status'],'no_106_witness':True}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();ts=time.time();base=p.read_words(ROOT/'inputs/M107.words.txt');stable=p.pcr_data(base);out={}
 for path in sorted((ROOT/'discovery/results').glob('adaptive*/result.json')):out[path.parent.name]=audit_adaptive(base,stable,v.load(path))
 out['triple_separation']=audit_triples(base,stable,v.load(ROOT/'discovery/results/triples/result.json.gz'))
 keys=[r['sha256_delta'] for run in out.values() for r in run['candidate_records']]
 out['totals']={'candidate_occurrences':len(keys),'distinct_cochains':len(set(keys)),'all_exactly_rejected':True};out['status']='PASS';out['seconds']=time.time()-ts
 if a.out:a.out.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:({i:j for i,j in x.items() if i!='candidate_records'} if isinstance(x,dict) else x) for k,x in out.items()},indent=2))
if __name__=='__main__':main()
