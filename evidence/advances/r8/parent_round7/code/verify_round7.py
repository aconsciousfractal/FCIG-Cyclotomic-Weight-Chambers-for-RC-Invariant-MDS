#!/usr/bin/env python3
"""Independent finite checker: literal DNA strings and Python integer/Fraction.
No solver, numerical array library, chart compiler, or discovery module is imported.
The general telescoping and rounding arguments are given in docs/REPORT.md.
"""
from __future__ import annotations
import argparse,gzip,hashlib,itertools,json,math,time
from collections import Counter,deque
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ALPHABET='ACGT';TRANS=str.maketrans('ACGT','TGCA')
V=[''.join(w) for w in itertools.product(ALPHABET,repeat=7)];VI={v:i for i,v in enumerate(V)}
E=[''.join(w) for w in itertools.product(ALPHABET,repeat=8)]
def rc(w):return w.translate(TRANS)[::-1]
SITES=[v for v in V if v<rc(v)];SI={v:i for i,v in enumerate(SITES)}
class Invalid(ValueError):pass
def require(b,why):
 if not b:raise Invalid(why)
def integer(x):return type(x) is int
def load(path):
 with (gzip.open(path,'rt',encoding='utf8') if str(path).endswith('.gz') else open(path,encoding='utf8')) as f:return json.load(f)
def digest(S):return hashlib.sha256(''.join(w+'\n' for w in sorted(S)).encode('ascii')).hexdigest()
def read_words(path):
 lines=Path(path).read_text(encoding='ascii').splitlines()
 require(len(lines)==len(set(lines)),'DUPLICATE_WORD')
 require(all(type(w) is str and len(w)==8 and not set(w)-set(ALPHABET) for w in lines),'WORD_DOMAIN')
 return set(lines)
def pcr_data(base):
 seen=set();stable=[];count=paired=0
 for w in E:
  if w in seen:continue
  P={w[i:]+w[:i] for i in range(8)};Q={rc(e) for e in P};seen|=P|Q
  require(len(P&base)==len(Q&base)==1,'PCR_CARDINALITY')
  count+=1 if P==Q else 2
  old=next(iter(P&base))
  if P==Q:
   fs=sorted(e for e in P if rc(e)==e);require(len(fs)==2 and old in fs,'STABLE_FIXED_DOMAIN')
   stable.append((old,next(e for e in fs if e!=old)))
  else:paired+=1;require(next(iter(Q&base))==rc(old),'PAIRED_SELECTION')
 require(count==8230 and paired==4051 and len(stable)==128,'PCR_CENSUS')
 return stable

def walk(w,base,closed=True):
 require(type(w) is list and len(w)>0,'EMPTY_WALK')
 require(all(type(e) is str and len(e)==8 and not set(e)-set(ALPHABET) for e in w),'WALK_EDGE_DOMAIN')
 require(all(a[1:]==b[:-1] for a,b in zip(w,w[1:])),'WALK_ADJACENCY')
 if closed:require(w[-1][1:]==w[0][:-1],'WALK_CLOSURE')
 return sum(e in base for e in w)
def verify_base(base,j):
 require(len(base)==8230 and {rc(w) for w in base}==base,'BASE_RC_CARDINALITY')
 stable=pcr_data(base)
 require(j['sha256_words']==digest(base),'BASE_DIGEST')
 rank=j['rank'];require(len(rank)==len(V) and all(integer(x) and 0<=x<=107 for x in rank),'BASE_RANK_DOMAIN')
 residual=0
 for e in E:
  if e not in base:residual+=1;require(rank[VI[e[1:]]]>rank[VI[e[:-1]]],'BASE_RANK_EDGE')
 require(residual==57306 and j['residual_edges']==residual,'RESIDUAL_COUNT')
 require(walk(j['longest_path'],base,False)==0 and len(j['longest_path'])==107,'BASE_LONGEST_PATH')
 require(len({e[:-1] for e in j['longest_path']}|{j['longest_path'][-1][1:]})==108,'BASE_PATH_SIMPLE')
 require(walk(j['critical_cycle'],base)==1 and len(j['critical_cycle'])==108,'BASE_CRITICAL_CYCLE')
 require(len({e[:-1] for e in j['critical_cycle']})==108,'BASE_CYCLE_SIMPLE')
 require(j['L']==107 and j['critical_cycle_hits']==1 and j['critical_cycle_length']==108 and j['component_lower_bound']==107 and j['vertices_processed']==16384 and j['necklaces']==8230,'BASE_METADATA')
 return stable

def generator(u,a):
 d=Counter()
 for e,k in [(u+a,1),(a+u,-1)]:d[e]+=k;d[rc(e)]+=k
 return {e:k for e,k in d.items() if k}
def cycle_data(cycles,base,stable):
 data=[]
 for cy in cycles:
  hits=walk(cy,base);co=Counter(cy);poss={min(v,rc(v)) for e in cy for v in [e[:-1],e[1:]]};sig={}
  for u in poss:
   z=tuple(sum(k*co[e] for e,k in generator(u,a).items()) for a in ALPHABET)
   require(sum(z)==0,'CYCLE_FULL_GRADIENT_SUM')
   if any(z):sig[u]=z
  gamma=[co[b]-co[a] for a,b in stable]
  data.append({'signature':sig,'gamma':gamma,'rhs':(len(cy)+106)//107-hits,'length':len(cy),'hits':hits})
 return data

def pair_system(pair,data):
 u,v=(V[i] for i in pair);A=[];b=[];gcds=[]
 for row in data:
  aa=list(row['signature'].get(u,(0,0,0,0))[:3])+list(row['signature'].get(v,(0,0,0,0))[:3])+row['gamma']
  g=math.gcd(*aa) or 1;A.append([x//g for x in aa]);b.append(-(-row['rhs']//g));gcds.append(g)
 for k in range(128):
  aa=[0]*134;aa[6+k]=1;A.append(aa);b.append(0);gcds.append(1)
  aa=[0]*134;aa[6+k]=-1;A.append(aa);b.append(-1);gcds.append(1)
 return A,b,gcds

def combination(A,b,mults):
 require(type(mults) is list and len(mults)>0,'EMPTY_DUAL')
 lhs=[Fraction(0)]*len(A[0]);rhs=Fraction(0);seen=set()
 for term in mults:
  require(type(term) is list and len(term)==3,'DUAL_FORMAT');i,n,d=term
  require(all(integer(x) for x in term) and 0<=i<len(A) and i not in seen and n>=0 and d>0,'DUAL_COEFFICIENT')
  seen.add(i);q=Fraction(n,d);rhs+=q*b[i]
  for k,a in enumerate(A[i]):
   if a:lhs[k]+=q*a
 return lhs,rhs

def verify_pair_proof(proof,data):
 pair=proof['sites'];require(type(pair) is list and len(pair)==2 and all(integer(i) and 0<=i<len(V) for i in pair) and pair[0]<pair[1],'PAIR_DOMAIN')
 require(all(V[i] in SI for i in pair),'PAIR_SITE_ORIENTATION')
 A,b,div=pair_system(pair,data)
 if proof['kind']=='gcd_farkas':
  lhs,rhs=combination(A,b,proof['multipliers']);require(not any(lhs) and rhs==1,'FARKAS_IDENTITY')
  return {'kind':'gcd_farkas','terms':len(proof['multipliers']),'gcd_steps':sum(div[i]>1 for i,n,d in proof['multipliers'])}
 require(proof['kind']=='integer_strip','UNKNOWN_PAIR_PROOF')
 form=[0]*134;seen=set()
 for k,a in proof['form']:
  require(integer(k) and integer(a) and 0<=k<134 and k not in seen and a!=0,'STRIP_FORM');seen.add(k);form[k]=a
 require(any(form),'STRIP_ZERO_FORM')
 vals=[]
 for key,sgn in [('lower',1),('upper_negative',-1)]:
  p=proof[key];v=p['value'];require(type(v) is list and len(v)==2 and all(integer(x) for x in v) and v[1]>0,'STRIP_VALUE')
  rhs_expected=Fraction(*v);lhs,rhs=combination(A,b,p['multipliers']);require(lhs==[sgn*x for x in form] and rhs==rhs_expected,'STRIP_DUAL_IDENTITY');vals.append(rhs)
 require(math.ceil(vals[0])>math.floor(-vals[1]),'STRIP_NOT_INTEGER_EMPTY')
 return {'kind':'integer_strip','form':proof['form'],'lower':str(vals[0]),'upper':str(-vals[1]),'sites':[V[i] for i in pair]}

def verify_proof(base,stable,j):
 require(j['schema']=='round7-at-most-two-sites-all-stable-v1','PROOF_SCHEMA')
 require(j['q']==4 and j['k']==8 and j['T']==107 and j['alphabet']=='ACGT' and j['complement']=='TGCA' and j['reference_letter']=='T','PROOF_SCOPE')
 require(j['gamma_domain']==[0,1] and j['site_coefficients']=='arbitrary integers' and j['gradient']=='arbitrary integer potential','COEFFICIENT_SCOPE')
 require(j['base_sha256']==digest(base),'PROOF_BASE')
 data=cycle_data(j['cycles'],base,stable);cover=j['cover_walks'];require(len(cover)==len(set(cover)) and cover,'COVER_INDICES');supports=[]
 for i in cover:
  require(integer(i) and 0<=i<len(data),'COVER_INDEX');r=data[i]
  require(r['rhs']>0 and all(a<=0 for a in r['gamma']),'COVER_GAMMA_MONOTONICITY');supports.append(set(r['signature']))
 # Every missing pair fails one verified closed walk for every allowed gamma.
 first=min(supports,key=len);pairs=set()
 for u in first:
  for v in SITES:
   if u!=v and all(u in s or v in s for s in supports):pairs.add(tuple(sorted([VI[u],VI[v]])))
 require(len(pairs)==519==j['expected_candidate_pairs'],'CANDIDATE_PAIR_CENSUS')
 require(len(SITES)==8192==j['expected_sites'] and math.comb(len(SITES),2)==33550336==j['expected_total_pairs'],'TOTAL_PAIR_CENSUS')
 require(len(stable)==j['expected_stable_pairs']==128,'STABLE_CENSUS')
 seen=set();summary=[]
 for p in j['pair_proofs']:
  pair=tuple(p['sites']);require(pair in pairs and pair not in seen,'PAIR_COVERAGE_KEY');seen.add(pair);summary.append(verify_pair_proof(p,data))
 require(seen==pairs,'PAIR_COVERAGE_INCOMPLETE')
 return {'status':'PASS','base_L':107,'target_L':106,'sites':len(SITES),'stable_bits':len(stable),'total_pairs':math.comb(len(SITES),2),'support_survivors':len(pairs),'cover_walks':len(cover),'certificate_walks':len(data),'farkas_cases':sum(x['kind']=='gcd_farkas' for x in summary),'integer_strip_cases':[x for x in summary if x['kind']=='integer_strip'],'unresolved_pairs':0,'coefficient_bounds':'none','gamma_domain':[0,1],'gradient':'arbitrary integer potential','retains_full_128_bit_cube':True}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();t=time.time();base=read_words(ROOT/'inputs/M107.words.txt');stable=verify_base(base,load(ROOT/'inputs/M107.json'));r=verify_proof(base,stable,load(ROOT/'certificates/two_sites_all_stable.json'));r['seconds']=time.time()-t
 if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()
