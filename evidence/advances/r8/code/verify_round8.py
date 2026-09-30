#!/usr/bin/env python3
"""Exact, solver-free Round 8 checker. All arithmetic is Python int/Fraction.
Imports the archived parent verifier, never any discovery/search module.
The parent handles the previously certified two-site implications. This file
reconstructs the new support lifting and the triple proof arithmetic separately.
"""
from __future__ import annotations
import argparse,gzip,hashlib,itertools,json,math,sys,time
from collections import Counter
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'parent_round7/code'))
import verify_round7 as p
class Invalid(ValueError):pass
def need(test,msg):
 if not test:raise Invalid(msg)
def integer(n):return type(n) is int
def load(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return json.load(f)

def rows_from_words(cycles,base,stable):
 """Recompute cycle incidence with literal strings, no numeric chart compiler."""
 out=[]
 for walk in cycles:
  hits=p.walk(walk,base);counts=Counter(walk)
  sites={min(w,p.rc(w)) for e in walk for w in (e[:-1],e[1:])};signature={}
  for u in sites:
   co=tuple(counts[u+a]-counts[a+u]+counts[p.rc(u+a)]-counts[p.rc(a+u)] for a in 'ACGT')
   need(sum(co)==0,'CLOSED_FULL_GRADIENT')
   if any(co):signature[u]=co
  out.append({'signature':signature,'gamma':[counts[b]-counts[a] for a,b in stable],'rhs':-(-len(walk)//107)-hits,'length':len(walk),'hits':hits})
 return out

def raw_row(case_sites,datum):
 return [a for u in case_sites for a in datum['signature'].get(u,(0,0,0,0))[:3]]+datum['gamma']

def forbidden_support(proof,data):
 old_sites=[p.V[i] for i in proof['sites']];groups=[proof['multipliers']] if proof['kind']=='gcd_farkas' else [proof['lower']['multipliers'],proof['upper_negative']['multipliers']]
 bad=set()
 for group in groups:
  divisors={i:(math.gcd(*raw_row(old_sites,data[i])) or 1) for i,n,d in group if i<len(data) and n}
  common=math.lcm(*(den*divisors.get(i,1) for i,n,den in group))
  accumulated={}
  for i,num,den in group:
   if i>=len(data) or not num:continue
   g=divisors[i];scale=common*num//(den*g)
   for u,vec in data[i]['signature'].items():
    if any(c%g for c in vec[:3]):bad.add(u)
    acc=accumulated.setdefault(u,[0,0,0])
    for j in range(3):acc[j]+=scale*vec[j]
  bad.update(u for u,vec in accumulated.items() if u not in old_sites and any(vec))
 need(not set(old_sites)&bad,'PARENT_SUPPORT_FORBIDDEN')
 return sorted(bad)

def minimal_hyperedges(edges):
 kept=[]
 for e in sorted(set(tuple(e) for e in edges),key=lambda t:(len(t),t)):
  s=frozenset(e)
  if not any(t<=s for t in kept):kept.append(s)
 return sorted((tuple(sorted(t)) for t in kept),key=lambda t:(len(t),t))

def transversal_triples(edges):
 """Complete recursion: branch on one unhit edge; depth at most three.
 An already hit all-edge state at depth<3 is separately reported, not dropped.
 """
 E=[frozenset(e) for e in edges];found=set();small=set()
 def rec(chosen,remaining):
  if not remaining:
   (found if len(chosen)==3 else small).add(tuple(sorted(chosen)));return
  if len(chosen)==3:return
  branch=min(remaining,key=lambda e:(len(e),tuple(sorted(e))))
  for u in sorted(branch):rec(chosen|{u},[e for e in remaining if u not in e])
 rec(set(),E)
 return sorted(small),sorted(found)

def verify_kernel(base,stable,j,source):
 need(j['schema']=='round8-lifted-support-kernel-v1' and j['T']==107,'KERNEL_SCHEMA_SCOPE')
 need(j['base_sha256']==p.digest(base),'KERNEL_BASE')
 data=rows_from_words(source['cycles'],base,stable)
 edges=[sorted(data[i]['signature']) for i in source['cover_walks']]
 need(len(j['lifted_supports'])==len(source['pair_proofs']),'LIFT_COVERAGE_COUNT')
 for i,(proof,rec) in enumerate(zip(source['pair_proofs'],j['lifted_supports'])):
  need(rec['pair_proof']==i,'LIFT_INDEX');e=forbidden_support(proof,data);need(e==rec['support'],'LIFTED_SUPPORT_MISMATCH');need(e,'EMPTY_SUPPORT_WOULD_BE_GLOBAL_UNSAT');edges.append(e)
 minimal=minimal_hyperedges(edges);need([list(e) for e in minimal]==j['minimal_edges'],'MINIMAL_HYPEREDGES')
 small,trips=transversal_triples(minimal)
 need([list(t) for t in small]==j['supports_below_three'],'LOWER_SUPPORTS');need(not small,'UNEXPECTED_SMALL_SUPPORT')
 need([list(t) for t in trips]==j['candidate_triples'],'TRIPLE_KERNEL_COMPLETENESS')
 need(j['total_sites']==8192 and j['total_triples']==math.comb(8192,3) and j['candidate_triples_count']==len(trips),'KERNEL_COUNTS')
 return trips,{'parent_lifted_proofs':len(source['pair_proofs']),'all_support_clauses':len(edges),'minimal_support_clauses':len(minimal),'total_three_site_supports':math.comb(8192,3),'candidate_triples':len(trips),'sites_in_triple_kernel':len(set().union(*map(set,trips))),'smaller_transversals':len(small)}

def normalized_row(sites,data,count,i):
 if i<count:
  raw=raw_row(sites,data[i]);div=math.gcd(*raw) or 1
  return [a//div for a in raw],-(-data[i]['rhs']//div),div
 j=i-count;need(0<=j<256,'ROW_INDEX');row=[0]*137;row[9+j//2]=1 if j%2==0 else -1
 return row,0 if j%2==0 else -1,1

def combine(sites,data,count,terms):
 need(type(terms) is list and len(terms)>0,'DUAL_EMPTY')
 seen=set();D=1
 for term in terms:
  need(type(term) is list and len(term)==3 and all(integer(x) for x in term),'DUAL_FORMAT')
  i,n,d=term;need(i>=0 and i<count+256 and i not in seen and n>=0 and d>0,'DUAL_DOMAIN');seen.add(i);D=math.lcm(D,d)
 lhs=[0]*137;rhs=0
 for i,n,d in terms:
  row,b,g=normalized_row(sites,data,count,i);k=n*(D//d);rhs+=k*b
  for j,a in enumerate(row):
   if a:lhs[j]+=k*a
 return lhs,rhs,D

def dominance_system(proof,data):
 count=proof['cycle_count'];need(integer(count) and 0<count<=len(data),'DOMINANCE_PREFIX')
 fixed=proof['fixed_gamma'];free=proof['free_gamma'];seen=set()
 need(type(fixed) is list and type(free) is list,'DOMINANCE_FORMAT')
 for item in fixed:
  need(type(item) is list and len(item)==2,'DOMINANCE_ITEM')
  j,b=item;need(integer(j) and 0<=j<128 and j not in seen and integer(b) and b in (0,1),'DOMINANCE_DOMAIN');seen.add(j)
  need(all(r['gamma'][j]<=0 if b==0 else r['gamma'][j]>=0 for r in data[:count]),'DOMINANCE_SIGN')
 need(all(integer(j) and 0<=j<128 and j not in seen for j in free) and free==sorted(set(free)),'FREE_GAMMA_DOMAIN')
 need(seen|set(free)==set(range(128)),'DOMINANCE_PARTITION')
 n=9+len(free);rows=[]
 for r in data[:count]:
  raw=[a for u in proof['sites'] for a in r['signature'].get(u,(0,0,0,0))[:3]]+[r['gamma'][j] for j in free]
  rhs=r['rhs']-sum(r['gamma'][j]*b for j,b in fixed);g=math.gcd(*raw) or 1
  rows.append(([a//g for a in raw],-(-rhs//g)))
 for j in range(len(free)):
  row=[0]*n;row[9+j]=1;rows.append((row,0));row=row.copy();row[9+j]=-1;rows.append((row,-1))
 return rows,n

def row_combination(rows,n,terms):
 need(type(terms) is list and bool(terms),'EXPLICIT_DUAL_EMPTY');D=1;seen=set()
 for t in terms:
  need(type(t) is list and len(t)==3 and all(integer(x) for x in t),'EXPLICIT_DUAL_FORMAT')
  i,a,b=t;need(0<=i<len(rows) and i not in seen and a>=0 and b>0,'EXPLICIT_DUAL_DOMAIN');seen.add(i);D=math.lcm(D,b)
 lhs=[0]*n;rhs=0
 for i,a,b in terms:
  k=a*(D//b);row,val=rows[i];rhs+=k*val
  lhs=[x+k*y for x,y in zip(lhs,row)]
 return lhs,rhs,D

def verify_dominated(proof,data):
 rows,n=dominance_system(proof,data);kind=proof['kind']
 if kind=='dominated_modular_equalities':
  mod=proof['modulus'];need(integer(mod) and mod>=2,'DOMINATED_MODULUS');lhs=[0]*n;rhs=0;seen=set()
  for eq in proof['equalities']:
   i=eq['row'];k=eq['coefficient'];need(integer(i) and 0<=i<len(rows) and i not in seen and integer(k) and 0<k<mod,'DOMINATED_EQUALITY_DOMAIN');seen.add(i)
   row,b=rows[i];co,val,D=row_combination(rows,n,eq['reverse_multipliers'])
   need(co==[-D*x for x in row] and val==-D*b,'DOMINATED_REVERSE_EQUALITY')
   lhs=[x+k*y for x,y in zip(lhs,row)];rhs+=k*b
  need(seen and all(x%mod==0 for x in lhs) and rhs%mod!=0,'DOMINATED_MODULAR_CONTRADICTION');return kind
 need(kind in ('dominated_gcd_farkas','dominated_cg_farkas'),'DOMINATED_KIND')
 if kind=='dominated_cg_farkas':
  need(type(proof['cuts']) is list,'CG_LIST')
  for cut in proof['cuts']:
   lhs,rhs,D=row_combination(rows,n,cut['multipliers'])
   need(all(x%D==0 for x in lhs),'CG_LEFT_NONINTEGRAL')
   row=[0]*n;seen=set()
   for j,a in cut['coefficients']:
    need(integer(j) and 0<=j<n and j not in seen and integer(a) and a!=0,'CG_COEFFICIENT');seen.add(j);row[j]=a
   need(row==[x//D for x in lhs],'CG_LEFT_IDENTITY')
   need(integer(cut['rhs']) and cut['rhs']==-(-rhs//D),'CG_ROUNDING')
   rows.append((row,cut['rhs']))
 lhs,rhs,D=row_combination(rows,n,proof['multipliers']);need(not any(lhs) and rhs==D,'DOMINATED_FINAL_FARKAS');return kind

def verify_case(proof,data):
 sites=proof['sites'];need(type(sites) is list and len(sites)==3 and sorted(set(sites))==sites and all(u in p.SI for u in sites),'TRIPLE_SCOPE')
 count=proof['cycle_count'];need(integer(count) and 0<count<=len(data),'CYCLE_PREFIX')
 kind=proof['kind']
 if kind.startswith('dominated_'):return verify_dominated(proof,data)
 if kind=='gcd_farkas':
  lhs,rhs,D=combine(sites,data,count,proof['multipliers']);need(not any(lhs) and rhs==D,'TRIPLE_FARKAS');return kind
 if kind=='integer_strip':
  form=[0]*137;seen=set()
  for j,a in proof['form']:
   need(integer(j) and 0<=j<137 and j not in seen and integer(a) and a,'STRIP_FORM');seen.add(j);form[j]=a
  need(any(form),'ZERO_STRIP_FORM');values=[]
  for key,sign in [('lower',1),('upper_negative',-1)]:
   r=proof[key];val=r['value'];need(len(val)==2 and all(integer(x) for x in val) and val[1]>0,'STRIP_VALUE')
   lhs,rhs,D=combine(sites,data,count,r['multipliers']);q=Fraction(*val);need(lhs==[sign*D*x for x in form] and rhs*q.denominator==D*q.numerator,'TRIPLE_STRIP_IDENTITY');values.append(q)
  need(math.ceil(values[0])>math.floor(-values[1]),'STRIP_NOT_EMPTY');return kind
 if kind=='modular_equalities':
  modulus=proof['modulus'];need(integer(modulus) and modulus>=2,'MODULUS')
  lhs=[0]*137;rhs=0;seen=set()
  for eq in proof['equalities']:
   i=eq['row'];a=eq['coefficient'];need(integer(i) and 0<=i<count+256 and i not in seen and integer(a) and 0<a<modulus,'EQUALITY_DOMAIN');seen.add(i)
   row,b,g=normalized_row(sites,data,count,i);co,vv,D=combine(sites,data,count,eq['reverse_multipliers']);need(co==[-D*x for x in row] and vv==-D*b,'REVERSE_EQUALITY_PROOF')
   lhs=[x+a*y for x,y in zip(lhs,row)];rhs+=a*b
  need(seen and all(x%modulus==0 for x in lhs) and rhs%modulus!=0,'MODULAR_CONTRADICTION');return kind
 raise Invalid('UNKNOWN_TRIPLE_PROOF')

def verify_cases(base,stable,trips,path):
 j=load(path);need(j['schema']=='round8-complete-three-site-refutations-v1' and j['T']==107 and j['base_sha256']==p.digest(base),'CASE_BASE_SCOPE')
 need(j['site_coefficients']=='arbitrary integers' and j['stable_bits']==128 and j['gamma_domain']==[0,1] and j['gradient']=='arbitrary integer potential','CASE_UNBOUNDED_SCOPE')
 data=rows_from_words(j['cycles'],base,stable);pool={tuple(t) for t in trips};seen=set();kinds=Counter()
 for proof in j['proofs']:
  t=tuple(proof['sites']);need(t in pool and t not in seen,'CASE_DUPLICATION_SCOPE');seen.add(t);kinds[verify_case(proof,data)]+=1
 missing=sorted(pool-seen)
 need([list(t) for t in missing]==j['unresolved_supports'],'RESIDUAL_SUPPORT_ACCOUNTING')
 need(not missing,'TRIPLE_COVERAGE_INCOMPLETE')
 need(dict(kinds)==j['proof_kinds'],'PROOF_KIND_ACCOUNTING')
 return {'certified_cases':len(seen),'remaining_cases':len(missing),'proof_kinds':dict(kinds),'closed_walks':len(data),'all_triples_excluded':not missing}

def verify_readable(base,stable,all_proofs,doc):
 proof=next(q for q in all_proofs['proofs'] if q['sites']==['AAACAAA','ACACAGA','CTACGCC'])
 data=rows_from_words(all_proofs['cycles'],base,stable);rows,n=dominance_system(proof,data);cut=proof['cuts'][0]
 lhs,rhs,D=row_combination(rows,n,cut['multipliers']);indices=sorted({i for i,a,d in cut['multipliers']}|{i for i,a,d in proof['multipliers'] if i<len(rows)})
 expected={'schema':'round8-readable-refutation-v1','sites':proof['sites'],'base_sha256':all_proofs['base_sha256'],'cycle_count':proof['cycle_count'],'row_indices_are_zero_based':True,'fixed_gamma':proof['fixed_gamma'],'free_gamma':proof['free_gamma'],'first_cut':cut,'pre_rounding_rhs':[rhs,D],'listed_rows':[],'final_multipliers':proof['multipliers']}
 for i in indices:
  row,b=rows[i];expected['listed_rows'].append({'index':i,'coefficients':[[j,a] for j,a in enumerate(row) if a],'rhs':b,'walk':all_proofs['cycles'][i] if i<proof['cycle_count'] else None})
 need(doc==expected,'READABLE_CERTIFICATE_MISMATCH')
 return {'status':'PASS','pre_rounding_rhs':f'{rhs}/{D}','rounded_rhs':cut['rhs'],'cuts_in_case':len(proof['cuts'])}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);ap.add_argument('--kernel-only',action='store_true');a=ap.parse_args();ts=time.time()
 base=p.read_words(ROOT/'inputs/M107.words.txt');stable=p.verify_base(base,load(ROOT/'inputs/M107.json'));source=load(ROOT/'parent_round7/certificates/two_sites_all_stable.json')
 parent=p.verify_proof(base,stable,source);trips,kernel=verify_kernel(base,stable,load(ROOT/'certificates/lifted_support_kernel.json'),source)
 result={'status':'PASS','base_L':107,'target_L':106,'parent_recheck':parent,'kernel':kernel}
 if not a.kernel_only:
  result['triple_refutations']=verify_cases(base,stable,trips,ROOT/'certificates/triple_refutations.json.gz')
  result['readable_case']=verify_readable(base,stable,load(ROOT/'certificates/triple_refutations.json.gz'),load(ROOT/'certificates/readable_refutation.json'))
 result['seconds']=time.time()-ts
 if a.out:a.out.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
