"""Literal-DNA verification of sequential domain eliminations.
No numerical solver or discovery code is imported. The supplied witnesses, not a
failed search, justify each removal. Partial pruning is never a refutation.
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import argparse, gzip, hashlib, json
from literal import (ROOT, WORDS, WORD_SET, rc, pcrs, groups, need, read_set,
                     height, rank_path, walk, save, load)

@lru_cache(None)
def universe():
    gs = {g['representative']: g for g in groups()}
    owner = {w: rep for rep,g in gs.items() for w in g['words']}
    opts={}
    selecting={}
    for p in pcrs():
        g=owner[p[0]]
        if gs[g]['stable'] or g in opts: continue
        options=tuple(sorted({tuple(sorted((w,rc(w)))) for w in p}))
        need(len(options)==len(p), 'OPTION_CENSUS')
        need(set().union(*map(set,options))==set(gs[g]['words']), 'OPTION_PARTITION')
        opts[g]=options
        for o in options:
            for w in o:
                need(w not in selecting,'OPTION_UNIQUENESS')
                selecting[w]=o
    return gs,owner,opts,selecting


def verify_trace(c, selected, raw_digest, partial=False):
    gs,owner,opts,selecting=universe()
    need(c.get('schema')=='R12-DOMAIN-TRACE-1','TRACE_SCHEMA')
    need(c.get('target_L')==106 and type(c.get('target_L')) is int,'TRACE_BOUND')
    need(c.get('base_sha256')==raw_digest,'TRACE_BASE_JOIN')
    free=c.get('free_groups')
    need(isinstance(free,list) and free==sorted(set(free)),'FREE_GROUP_ORDER')
    need(all(g in opts for g in free),'FREE_NONSTABLE_GROUPS')
    free=set(free)
    remaining={g:set(opts[g]) for g in free}
    steps=c.get('steps')
    need(isinstance(steps,list) and len(steps)>0,'TRACE_STEPS')
    boundary=set(); types=Counter(); assumptions=set(); occurrence_count=0
    for i,s in enumerate(steps):
        g=s.get('group'); o=s.get('option')
        need(g in free,'STEP_FREE_GROUP')
        need(isinstance(o,list) and len(o)==2 and all(w in WORD_SET for w in o),
             'STEP_OPTION_FORMAT')
        o=tuple(o)
        need(o in remaining[g], 'STEP_OPTION_CURRENT')
        assumptions.add(g)
        kind=s.get('kind'); path=s.get('walk')
        need(kind in ('cycle','path'),'WITNESS_KIND')
        walk(path,closed=(kind=='cycle'),simple=True)
        if kind=='path':need(len(path)==107,'LONG_PATH_BOUND')
        types[kind]+=1;occurrence_count+=len(path)
        for w in path:
            j=owner[w]
            need(not gs[j]['stable'],'STABLE_EDGE_CANNOT_BE_MANDATORY')
            if j==g:
                need(w not in o,'TESTED_OPTION_CUTS_WITNESS')
            elif j in free:
                # Its unique selecting option has already been refuted. No
                # reference to a future removal or to solver propagation.
                need(selecting[w] not in remaining[j], 'UNPROVED_RESIDUAL_EDGE')
            else:
                need(w not in selected,'FROZEN_SELECTED_WITNESS')
                boundary.add(min(w,rc(w)))
        remaining[g].remove(o)
    if partial:
        need(c.get('result')=='PARTIAL_PRUNING','PARTIAL_STATUS')
        need(all(remaining.values()),'PARTIAL_HAS_EMPTY_DOMAIN')
        rebuilt={g:[list(o) for o in sorted(remaining[g])] for g in sorted(free)}
        got=c.get('remaining_options')
        need(isinstance(got,dict) and set(got)==set(rebuilt),'REMAINING_GROUP_JOIN')
        for g in rebuilt:
            need(isinstance(got[g],list) and all(isinstance(o,list) and len(o)==2
                 and all(w in WORD_SET for w in o) for o in got[g]),'REMAINING_FORMAT')
            need(len(got[g])==len(set(map(tuple,got[g]))),'REMAINING_DUPLICATE')
            need(sorted(got[g])==rebuilt[g],'REMAINING_OPTIONS_JOIN')
        return {'status':'PASS','result':'CERTIFIED_PARTIAL_PRUNING',
                'free_nonstable_groups':len(free),'removed_options':len(steps),
                'remaining_options':sum(map(len,remaining.values())),
                'witness_types':dict(types),
                'claim_boundary':'No refutation; no feasible endpoint; maximality of pruning not asserted.'}
    root=c.get('empty_group')
    need(root in free and not remaining[root], 'ROOT_DOMAIN_NOT_EMPTY')
    need(c.get('result')=='EMPTY_DOMAIN','REFUTATION_STATUS')
    need(c.get('boundary_words')==sorted(boundary),'BOUNDARY_EXACT_JOIN')
    boundary_groups=sorted({owner[w] for w in boundary})
    need(c.get('boundary_groups')==boundary_groups,'BOUNDARY_GROUP_JOIN')
    need(c.get('single_assumptions')==sorted(assumptions),'SINGLE_ASSUMPTIONS_JOIN')
    need(root in assumptions,'ROOT_ASSUMPTION')
    # With the boundary words unselected, the same replay uses no information
    # about the outside indicators. Other free groups used as dependencies have
    # a removal in the trace and hence appear in assumptions.
    for s in steps:
        for w in s['walk']:
            j=owner[w]
            if j in free:need(j in assumptions,'UNDECLARED_SINGLE_ASSUMPTION')
            else:need(min(w,rc(w)) in boundary,'UNDECLARED_BOUNDARY_ASSUMPTION')
    need(all(w not in selected and rc(w) not in selected for w in boundary),'BOUNDARY_BASE_ZERO')
    base_excess=sum(len(set(gs[g]['words']) & set(selected))//2-1 for g in assumptions)
    need(base_excess>=1,'BASE_MUST_PAY_SURPLUS')
    return {'status':'PASS','result':'EMPTY_DOMAIN','free_nonstable_groups':len(free),
            'steps':len(steps),'witness_types':dict(types),'checked_word_occurrences':occurrence_count,
            'boundary_RC_literals':len(boundary),'boundary_nonstable_groups':len(boundary_groups),
            'single_assumptions':len(assumptions),'cycles_only':types['path']==0,
            'base_lifted_lhs':base_excess,'min_cardinality_at_L106_in_fixed_region':8232}


def verify_nonstable_path(c, m):
    gs,owner,_,_=universe()
    need(c.get('kind')==1,'NONSTABLE_PATH_KIND')
    p=c.get('walk');walk(p,simple=True)
    need(len(p)==107,'NONSTABLE_PATH_LENGTH')
    need(not(set(p)&m),'NONSTABLE_PATH_SELECTED')
    need(all(not gs[owner[w]]['stable'] for w in p),'NONSTABLE_PATH_STABLE')
    return {'status':'PASS','length':107,'stable_PCR_edges':0,
            'claim':'All stable-PCR changes alone leave this fixed nonstable path.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    basefile=ROOT/'inputs/S8232.words.txt'
    s=read_set(basefile,8232)
    stable_edges=set().union(*(g['words'] for g in groups() if g['stable']))
    need(len(stable_edges)==988 and sum(rc(w)==w for w in stable_edges)==256,'STABLE_CENSUS')
    _,hresult=height(load(ROOT/'inputs/S8232.height.json'),s)
    rresult=rank_path(load(ROOT/'inputs/S8232.rank_path.json'),s)
    digest=hashlib.sha256(basefile.read_bytes()).hexdigest()
    rows={};data={}
    for p in sorted((ROOT/'certificates/regions').glob('*.json')):
        c=load(p);need(c.get('label')==p.stem,'CERTIFICATE_LABEL');data[p.stem]=c
        rows[p.stem]=verify_trace(c,s,digest)
    need(len(rows)==8,'REGION_COUNT')
    for name,n in [('adaptive_1',1),('adaptive_13',13),('adaptive_139',139),
                   ('geometric_64',64),('geometric_128',128),('geometric_256',256),
                   ('geometric_512',512),('geometric_1024',1024)]:
        need(rows[name]['free_nonstable_groups']==n,'REGION_SIZE')
    for old,new in [('adaptive_1','adaptive_13'),('adaptive_13','adaptive_139')]:
        need(set(data[new]['free_groups']) == set(data[old]['free_groups']) |
             set(data[old]['boundary_groups']),'ADAPTIVE_EXPANSION_JOIN')
    c=json.loads(gzip.decompress((ROOT/'certificates/adaptive_686_pruning.json.gz').read_bytes()))
    need(set(c['free_groups']) == set(data['adaptive_139']['free_groups']) |
         set(data['adaptive_139']['boundary_groups']),'ADAPTIVE_FINAL_JOIN')
    pruning=verify_trace(c,s,digest,partial=True)
    need(pruning['free_nonstable_groups']==686 and pruning['removed_options']==2083,'PRUNING_COUNTS')
    m=read_set(ROOT/'inputs/M107.words.txt',8230)
    fixed_path=verify_nonstable_path(load(ROOT/'certificates/nonstable_path.json'),m)
    clauses=[set(c['boundary_words']) for c in data.values()]
    result={'status':'PASS','stable_PCR_arcs_deleted':988,'RC_fixed_words':256,'base':rresult,'height':hresult,'regions':rows,
            'partial':pruning,'frozen_nonstable_obstruction':fixed_path,
            'distinct_boundary_clauses':len({tuple(sorted(b)) for b in clauses}),
            'boundary_antichain_count':sum(not any(t<b for t in clauses) for b in clauses)}
    save(a.out/'domain_verification.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
