#!/usr/bin/env python3
"""Exact audit of the two restricted three-site searches. Standard library.
A generated signed cochain is not an MDS. A closed walk is the rejection object.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse,json,time,hashlib
import verify_round7 as v
ROOT=Path(__file__).resolve().parents[1]

def audit_run(base,stable,j):
 v.require(j['schema']=='restricted-three-site-cegar-v1' and j['T']==107,'SEARCH_SCOPE')
 v.require(j['z_bounds']==[-2,2] and j['active_sites_upper']==3 and j['stable_bits_free']==128,'SEARCH_DOMAIN')
 v.require(j['initial_cycles']==193 and j['cycles'][:193]==v.load(ROOT/'inputs/discovery_initial_cycles.json'),'INITIAL_CUT_CORPUS')
 sites=j['sites'];v.require(len(sites)==95 and sites==sorted(set(sites)) and all(v.integer(i) and 0<=i<len(v.V) and v.V[i] in v.SI for i in sites),'SEARCH_SITE_POOL')
 terms=[v.generator(v.V[u],a) for u in sites for a in 'ACG']+[Counter({a:-1,b:1}) for a,b in stable]
 for cy in j['cycles']:v.walk(cy,base)
 accepted=[];gamma_hist=Counter();cut_lengths=[];seen=set()
 for i,rec in enumerate(j['audit']):
  v.require(rec['iteration']==i,'SEARCH_ITERATION')
  if not rec['incumbent']:
   v.require(i==len(j['audit'])-1,'NO_INCUMBENT_NOT_TERMINAL');continue
  zz={}
  for ind,k in rec['operations']:
   v.require(v.integer(ind) and v.integer(k) and 0<=ind<len(terms) and ind not in zz,'OPERATION_INDEX')
   v.require((-2<=k<=2 and k!=0) if ind<285 else k==1,'OPERATION_DOMAIN');zz[ind]=k
  delta=Counter()
  for ind,k in zz.items():
   for w,a in terms[ind].items():delta[w]+=k*a
  delta={w:a for w,a in delta.items() if a}
  v.require(all(delta.get(v.rc(w),0)==a for w,a in delta.items()),'SIGNED_RC')
  per=Counter()
  for w,a in delta.items():per[min(w[h:]+w[:h] for h in range(8))]+=a
  v.require(not any(per.values()),'SIGNED_PCR_BALANCE')
  actual=len({ind//3 for ind in zz if ind<285});bits=sum(k for ind,k in zz.items() if ind>=285)
  v.require(actual<=3 and rec['actual_sites']==actual and rec['stable_changes']==bits,'ACTUAL_SUPPORT_METADATA')
  if j.get('actual_stable_change_required',False):v.require(bits>=1,'ACTUAL_STABLE_CHANGE_REQUIRED')
  v.require(rec['clock_status']==0,'UNEXPECTED_SEARCH_OUTCOME');idx=rec['cycle_index']
  v.require(v.integer(idx) and idx==j['initial_cycles']+i and idx<len(j['cycles']),'SEPARATOR_INDEX')
  def hits(cy):return sum(int(w in base)+delta.get(w,0) for w in cy)
  for cy in j['cycles'][:idx]:v.require(hits(cy)>=(len(cy)+106)//107,'PREVIOUS_CYCLE_CUT')
  cy=j['cycles'][idx];cost=len(cy)-107*hits(cy);v.require(cost>0,'NONPOSITIVE_REJECTION_WALK')
  key=hashlib.sha256(json.dumps(sorted(delta.items()),separators=(',',':')).encode()).hexdigest();v.require(key not in seen,'REPEATED_SIGNED_CANDIDATE');seen.add(key)
  accepted.append({'iteration':i,'sha256_delta':key,'sites':actual,'stable_changes':bits,'witness_length':len(cy),'witness_hits':hits(cy),'positive_excess':cost});gamma_hist[bits]+=1;cut_lengths.append(len(cy))
 v.require(len(j['cycles'])==j['initial_cycles']+len(accepted),'CUT_CANDIDATE_JOIN')
 v.require(j['status'] in ['WALL_BUDGET','ITERATION_LIMIT','NO_INCUMBENT','SOLVER_UNSAT_UNCERTIFIED'],'SEARCH_TERMINAL_STATUS')
 return {'status':'PASS','search_status':j['status'],'candidates':len(accepted),'site_pool':len(sites),'site_upper_bound':3,'coefficient_bounds':[-2,2],'gamma_histogram':dict(gamma_hist),'all_rejected_by_exact_closed_walk':True,'candidate_records':accepted,'search_seconds':j['seconds'],'not_global_unsat':True}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();t=time.time();base=v.read_words(ROOT/'inputs/M107.words.txt');stable=v.pcr_data(base);out={}
 for name in ['three_sites','three_sites_stable']:out[name]=audit_run(base,stable,v.load(ROOT/'discovery/results'/(name+'.json')))
 keys=[r['sha256_delta'] for name in ['three_sites','three_sites_stable'] for r in out[name]['candidate_records']]
 v.require(len(keys)==len(set(keys)),'CROSS_RUN_DUPLICATE');out['joint_distinct_candidates']=len(keys)
 out['seconds']=time.time()-t;out['status']='PASS'
 if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:({a:b for a,b in r.items() if a!='candidate_records'} if isinstance(r,dict) else r) for k,r in out.items()},indent=2))
if __name__=='__main__':main()
