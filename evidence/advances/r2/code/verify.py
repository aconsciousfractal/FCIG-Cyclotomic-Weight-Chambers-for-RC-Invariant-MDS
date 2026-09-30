#!/usr/bin/env python3
"""Exact Round 2 checker. Python standard library only; no solver or search imports.

Graph and retiming checks use the literal-string checker from Round 1.
The frozen-PCR domain reconstruction and local refutation checker are new.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
from math import prod
import argparse, json, time
import graph_checker as g
ROOT=Path(__file__).resolve().parents[1]

def load(p:Path):
    return json.loads(p.read_text(encoding='utf-8'))

def local_exclusion(base:set[str], cert:dict) -> dict:
    g.require(cert['schema']=='rc-mds-frozen-PCR-propagation-v1','LOCAL_SCHEMA')
    g.require(g.digest_set(base)==cert['baseline_sha256'],'LOCAL_BASE_DIGEST')
    B=cert['bound'];g.require(type(B) is int and B==120,'LOCAL_BOUND')
    _,_,hits=g.walk_check(cert['trigger_cycle'],base,cycle=True)
    g.require(len(cert['trigger_cycle'])==122 and hits==1,'LOCAL_TRIGGER')
    U=cert['active_edges'];g.require(isinstance(U,list) and U==sorted(set(U)),'LOCAL_ACTIVE_CANONICAL')
    g.require(all(isinstance(w,str) and len(w)==8 and set(w)<=set(g.ALPHABET) for w in U),'LOCAL_ACTIVE_WORD')
    active=set(U)
    closure={z[i:]+z[:i] for w in cert['trigger_cycle'] for z in (w,g.rc(w)) for i in range(8)}
    g.require(active==closure,'LOCAL_ACTIVE_CLOSURE')
    # Rebuild EVERY RC-compatible selection on the complete active PCR family.
    seen=set();orbits=[];lookup={}
    for w in sorted(active):
        if w in seen:continue
        orb={w[i:]+w[:i] for i in range(8)}
        g.require(orb<=active,'LOCAL_ROTATION_CLOSURE')
        for z in orb:lookup[z]=len(orbits)
        orbits.append(orb);seen.update(orb)
    g.require(len(orbits)==cert['active_necklaces']==99,'LOCAL_NECKLACE_COUNT')
    groups=[];covered=set()
    for j,orb in enumerate(orbits):
        if j in covered:continue
        mate=lookup[g.rc(min(orb))];covered.update((j,mate));opts=[]
        for w in sorted(orb):
            if mate==j and g.rc(w)!=w:continue
            opts.append(sorted({w,g.rc(w)}))
        g.require(len(opts)>0,'LOCAL_EMPTY_INITIAL_DOMAIN')
        groups.append(sorted(opts))
    g.require(groups==cert['groups'],'LOCAL_DOMAIN_COMPLETENESS')
    contributions=[]
    for cut in cert['cuts']:
        kind=cut['kind'];ws=cut['words']
        g.require(kind in ('path','cycle'),'LOCAL_CUT_KIND')
        g.walk_check(ws,base,cycle=kind=='cycle')
        if kind=='path':
            g.require(len(ws)==B+1,'LOCAL_PATH_LENGTH');rhs=1
        else:rhs=(len(ws)+B)//(B+1)
        g.require(type(cut['rhs']) is int and cut['rhs']==rhs,'LOCAL_CUT_RHS')
        counts=Counter(ws)
        fixed=sum(n for w,n in counts.items() if w not in active and w in base)
        terms=[[sum(counts[w] for w in opt) for opt in opts] for opts in groups]
        contributions.append((fixed,terms,rhs))
    domains=[set(range(len(opts))) for opts in groups]
    for step in cert['eliminations']:
        c,h,o=step['cut'],step['group'],step['option']
        g.require(all(type(v) is int for v in (c,h,o)),'LOCAL_STEP_TYPE')
        g.require(0<=c<len(contributions) and 0<=h<len(groups),'LOCAL_STEP_INDEX')
        g.require(o in domains[h],'LOCAL_OPTION_ALREADY_ABSENT')
        fixed,terms,rhs=contributions[c]
        upper=fixed+terms[h][o]+sum(max(vals[x] for x in domains[j]) for j,vals in enumerate(terms) if j!=h)
        g.require(type(step['upper']) is int and upper==step['upper'],'LOCAL_UPPER_RECOMPUTATION')
        g.require(upper<rhs,'LOCAL_UNJUSTIFIED_ELIMINATION')
        domains[h].remove(o)
        g.require(len(domains[h])>0,'LOCAL_EMPTY_DOMAIN_NEEDS_EXPLICIT_CONTRADICTION')
    final=cert['contradiction'];c=final['cut']
    g.require(type(c) is int and 0<=c<len(contributions),'LOCAL_FINAL_INDEX')
    fixed,terms,rhs=contributions[c]
    upper=fixed+sum(max(vals[x] for x in domains[j]) for j,vals in enumerate(terms))
    g.require(type(final['upper']) is int and final['upper']==upper and final['rhs']==rhs,'LOCAL_FINAL_RECOMPUTATION')
    g.require(upper<rhs,'LOCAL_NO_CONTRADICTION')
    return {'status':'EXACT_FROZEN_NEIGHBORHOOD_UNSAT','bound':B,
            'active_PCRs':len(orbits),'active_edges':len(active),'RC_domain_groups':len(groups),
            'initial_options':sum(map(len,groups)),'initial_assignments':prod(map(len,groups)),
            'validated_cuts':len(contributions),'eliminations':len(cert['eliminations']),
            'terminal_upper':upper,'terminal_rhs':rhs,'unrestricted_global_refutation':False}

def parity_obstruction(base:set[str],cert:dict)->dict:
    T,eps=cert['T'],cert['epsilon']
    g.require(type(T) is int and T>=1 and type(eps) is int and eps in (0,1),'PARITY_PARAMETERS')
    P,Q=cert['path_to_mate'],cert['path_back']
    u,v,hp=g.walk_check(P,base);v2,u2,hq=g.walk_check(Q,base)
    g.require(u==cert['vertex'] and v==g.rc(u) and (v2,u2)==(v,u),'PARITY_ENDPOINTS')
    p,q=len(P)-T*hp,len(Q)-T*hq;c=T-1-T*eps
    g.require(p+q==0,'PARITY_TIGHT_RETURN')
    g.require((c-p)%2==1,'PARITY_INTEGER_CONTRADICTION')
    return {'T':T,'epsilon':eps,'P_length':len(P),'P_hits':hp,'Q_length':len(Q),'Q_hits':hq,
            'forced_difference':p,'impossible_2H_value':c-p,'status':'EXACT_SECTOR_OBSTRUCTION'}

def separation(root:Path)->dict:
    baseline=g.read_set(root/'witnesses/M121.words.txt');old=g.read_set(root/'witnesses/M131.words.txt');final=g.read_set(root/'witnesses/M114.words.txt')
    counts={}
    for name,base in [('M121',baseline),('M131',old)]:
        cy=load(root/f'certificates/{name}.json')['critical_cycle']
        _,_,x=g.walk_check(cy,base,cycle=True);_,_,y=g.walk_check(cy,final,cycle=True)
        g.require(x!=y,'NO_COMPONENT_SEPARATION')
        counts[name]={'cycle_length':len(cy),'baseline_hits':x,'M114_hits':y}
    active=set(load(root/'certificates/local_exclusion_99PCR.json')['active_edges'])
    diff=baseline^final
    def necklace(w):return min(w[i:]+w[:i] for i in range(8))
    return {'cycle_separators':counts,'symmetric_difference_words':len(diff),
            'changed_PCRs':len({necklace(w) for w in diff}),
            'changed_PCRs_outside_local_neighborhood':len({necklace(w) for w in diff if w not in active})}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--quick',action='store_true',help='Skip only the 1145-step chain replay.');p.add_argument('--manifest',action='store_true');a=p.parse_args();root=a.root;start=time.time()
    summary={'status':'PASS','sets':{}}
    for name,entry in load(root/'certificates/catalog.json').items():
        summary['sets'][name]=g.set_certificate_check(g.read_set(root/entry['words']),load(root/entry['certificate']))
    ret=load(root/'certificates/retiming_catalog.json')
    for entry in ret:g.retiming_check(g.read_set(root/entry['base']),g.read_set(root/entry['target']),load(root/entry['certificate']))
    summary['retiming_certificates_verified']=len(ret)
    summary['parity_obstructions']=[parity_obstruction(g.read_set(root/e['base']),load(root/e['certificate'])) for e in load(root/'certificates/parity_catalog.json')]
    summary['local_exclusion']=local_exclusion(g.read_set(root/'witnesses/M121.words.txt'),load(root/'certificates/local_exclusion_99PCR.json'))
    summary['separation']=separation(root)
    if not a.quick:
        chain=load(root/'certificates/construction_chain.json');summary['chain']=g.chain_check(root,chain)
        prev=121;counts={'strict_increases':0,'equal_length_steps':0,'strict_decreases':0}
        for s in chain['steps']:
            L=s['L_after'];counts['strict_increases' if L>prev else 'strict_decreases' if L<prev else 'equal_length_steps']+=1;prev=L
        summary['chain']['length_step_counts']=counts
    else:summary['chain']={'status':'NOT_RUN_QUICK_MODE'}
    cy=load(root/'certificates/M114.json')['critical_cycle'];g.require(len(cy)==115,'FINAL_CYCLE_LENGTH')
    g.require((root/'witnesses/M114_critical_cycle_115.txt').read_text().splitlines()==cy,'CYCLE_TEXT_MISMATCH')
    word=(root/'witnesses/M114_longest_avoiding_word.txt').read_text().strip();path=load(root/'certificates/M114.json')['longest_path']
    g.require(word==path[0]+''.join(w[-1] for w in path[1:]) and len(word)==121,'AVOIDING_WORD_TEXT')
    if a.manifest:summary['manifest_files']=g.manifest_check(root)
    summary['elapsed_seconds']=time.time()-start
    print(json.dumps(summary,indent=2,allow_nan=False))
if __name__=='__main__':main()
