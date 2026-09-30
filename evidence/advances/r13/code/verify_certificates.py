"""Exact positive witnesses, complete cost-two augmentation census, and cycles.
Only literal checker primitives are imported. No numerical status is a premise.
"""
from literal import *
from itertools import combinations
import argparse

def endpoint_ref(p):return check_endpoint(load(safe_path(p)))

def check_panel(c):
    need(c.get('schema')=='R13-AUGMENTATION-PANEL-1','PANEL_SCHEMA')
    items=c.get('items');need(type(items) is list and len(items)==14,'PANEL_SIZE')
    need([r.get('source_index') for r in items]==list(range(14)),'PANEL_INDEX')
    hashes=set();out=[]
    for r in items:
        s,d,b=endpoint_ref(r['endpoint']);need(len(s)==8230 and max(d.values())==107,'PANEL_SOURCE_MDS')
        need(digest(s)==r['source'].get('sha256'),'SOURCE_JOIN');need(digest(s) not in hashes,'SOURCE_DUPLICATE');hashes.add(digest(s))
        need(r.get('target_L')==106 and r.get('max_additions')==2,'CENSUS_SCOPE')
        K={w for w in W if w not in s and d[w[:-1]]+1+b[w[1:]]==107}
        need(r.get('critical_words')==sorted(K),'CRITICAL_SUPPORT')
        need({rc(w) for w in K}==K,'CRITICAL_RC')
        fix=sorted(w for w in K if w==rc(w));pairs=sorted(w for w in K if w<rc(w))
        need((len(fix),len(pairs))==(r.get('critical_fixed'),r.get('critical_pairs')),'CRITICAL_COUNTS')
        candidates=[frozenset()]+[frozenset([w]) for w in fix]+[frozenset(p) for p in combinations(fix,2)]+[frozenset([w,rc(w)]) for w in pairs]
        need(len(candidates)==len(set(candidates))==r.get('candidate_count_including_empty'),'CANDIDATE_CENSUS')
        paths=r.get('negative_paths');need(type(paths) is list and paths,'NEGATIVE_PATH_LIST');ps=[]
        for p in paths:
            walk(p);need(len(p)==107 and not(set(p)&s),'NEGATIVE_PATH_VALID');ps.append(frozenset(p))
        need(len(ps)==len(set(ps)),'NEGATIVE_PATH_DUPLICATE')
        uncovered={a for a in candidates if not any(a.isdisjoint(p) for p in ps)}
        good=set()
        for a in r['successes']:
            words=a.get('additions');need(type(words) is list and words==sorted(set(words)),'SUCCESS_FORMAT');A=frozenset(words)
            need(A in candidates and A not in good,'SUCCESS_CANDIDATE')
            n,nd,nb=endpoint_ref(a['endpoint']);need(n==s|A and max(nd.values())<=106,'SUCCESS_ENDPOINT');good.add(A)
        need(uncovered==good,'CENSUS_COMPLETE_COVERAGE')
        out.append({'index':r['source_index'],'critical_edges':len(K),'critical_fixed':len(fix),'critical_pairs':len(pairs),'candidates_including_empty':len(candidates),'negative_paths':len(paths),'successes':len(good)})
    return out

def check_bridge(c):
    need(c.get('schema')=='R13-RETIMING-BRIDGE-1','BRIDGE_SCHEMA')
    m0,_,_=endpoint_ref(c['source0']);m1,_,_=endpoint_ref(c['source1'])
    f=c.get('f');need(type(f) is list and len(f)==len(V) and all(type(v) is int for v in f),'GRADIENT_FORMAT');f=dict(zip(V,f))
    for w in W:need(int(w in m1)-int(w in m0)==f[w[:-1]]-f[w[1:]],'GRADIENT_EQUATION')
    need(all(f[v]+f[rc(v)]==c.get('mate_sum') for v in V),'GRADIENT_RC_SUM')
    p=c.get('cycle');walk(p,closed=True)
    need(len(p)==108 and c.get('target_L')==106 and c.get('expected_cycle_weight')==1,'BRIDGE_TARGET')
    need(not any(w==rc(w) for w in p),'CYCLE_HAS_FIXED_ARC')
    need(len(set(p)&m0)==len(set(p)&m1)==1,'BRIDGE_CYCLE_INTERSECTION')
    A=frozenset(c.get('attaining_additions',[]));need(len(A)==2 and {rc(w) for w in A}==A,'BRIDGE_ADDITION_RC')
    for m in [m0,m1]:need(not(m&A) and max(graph(m|A)[0].values())==105,'BRIDGE_ADDITION_POSITIVE')
    need(c.get('minimum_RC_additions_over_gradient_class')==2,'BRIDGE_MINIMUM')
    return {'gradient_checked_edges':len(W),'gradient_range':[min(f.values()),max(f.values())],'cycle_length':len(p),'cycle_fixed_arcs':0,'augmentation_minimum_at_bound_106':2}

def check_affine(c):
    need(c.get('schema')=='R13-AFFINE-PRUNING-1','AFFINE_SCHEMA')
    items=c.get('items');need(type(items) is list and [z.get('tag') for z in items]==['old','alt2','alt3'],'AFFINE_SCOPE');out=[]
    for r in items:
        s,_,_=endpoint_ref(r['endpoint']);need(len(s)==8232,'AFFINE_BASE_COST')
        over=[g for g in groups() if not g['stable'] and len(g['words']&s)==4]
        need(len(over)==1 and r.get('group')==over[0]['representative'],'AFFINE_OVERGROUP')
        for g in groups():need(len(g['words']&s)==(4 if g==over[0] else 1 if g['stable'] else 2),'AFFINE_ONLY_SURPLUS')
        options={tuple(sorted([w,rc(w)])) for w in over[0]['pcr']}
        cases=r.get('cases');need(type(cases) is list and len(cases)==len(options),'AFFINE_CASE_COUNT')
        need({tuple(z['subtract_pair']) for z in cases}==options,'AFFINE_OPTION_COVERAGE')
        need(r.get('T')==107,'AFFINE_T');weights=[];lengths=[]
        for z in cases:
            removed=z['subtract_pair'];need(tuple(removed) in options,'AFFINE_PAIR');P=z['cycle'];walk(P,closed=True)
            weight=sum(int(w in s)-int(w in removed) for w in P)
            need(type(z.get('cycle_weight')) is int and weight==z['cycle_weight'],'AFFINE_WEIGHT')
            violation=len(P)-107*weight
            need(violation>0 and violation==z.get('violation_at_T107'),'AFFINE_VIOLATION')
            weights.append(weight);lengths.append(len(P))
        if r['tag']=='old':need(r.get('claims_no_acyclic_boolean_endpoint') is True and max(weights)<=0,'AFFINE_OLD_STRONG')
        else:
            need(r.get('claims_no_acyclic_boolean_endpoint') is False,'AFFINE_NEW_SCOPE')
            need(sum(w==1 for w in weights)==1 and all(w<=1 for w in weights),'AFFINE_NEW_WEIGHTS')
            z=next(z for z in cases if z['cycle_weight']==1)
            need(len(z['cycle'])==108 and set(z['subtract_pair'])<=s,'AFFINE_POSITIVE_CASE')
            n=s-set(z['subtract_pair']);need(len(n)==8230 and max(graph(n)[0].values())==107,'AFFINE_MINIMUM_ATTAINED')
            ref='certificates/M107_panel_'+('02' if r['tag']=='alt2' else '03')+'.words.txt'
            need(n==read_set(safe_path(ref)),'AFFINE_SOURCE_CLASS_JOIN')
        out.append({'tag':r['tag'],'cases':len(cases),'cycle_lengths':lengths,'cycle_weights':weights,'all_acyclic_endpoints_excluded':r['tag']=='old','minimum_L_in_affine_pruning_family':None if r['tag']=='old' else 107,'smaller_RC_decycling_endpoints_class':'EMPTY' if r['tag']=='old' else 'GRADIENT_CLASS_OF_PANEL_'+('02' if r['tag']=='alt2' else '03')})
    return out

def check_joins(c):
    need(c.get('schema')=='R13-ESCAPE-JOINS-1','JOINS_SCHEMA')
    paths=sorted(p.relative_to(ROOT).as_posix() for p in (ROOT/'upstream/round12/certificates/regions').glob('*.json'))
    need([r['parent_certificate'] for r in c.get('items',[])]==paths,'JOINS_PARENTS')
    pmap={w:p for p in pcrs() for w in p};old=read_set(ROOT/'inputs/S8232.words.txt');out=[]
    for r in c['items']:
        parent=load(safe_path(r['parent_certificate']));sums=[]
        need([z['endpoint'] for z in r['endpoints']]==['certificates/S8232_reuse_02.rank_path.json','certificates/S8232_reuse_03.rank_path.json'],'JOINS_ENDPOINTS')
        for z in r['endpoints']:
            s,_,_=endpoint_ref(z['endpoint']);hits=[w for w in parent['boundary_words'] if w in s]
            surplus=sum(len(s&set(pmap[g]))-1 for g in parent['single_assumptions'])
            need(hits==z['boundary_hits'] and len(hits)==z['boundary_count'],'JOINS_HITS')
            need(surplus==z['surplus_term'] and len(hits)+surplus>=1,'JOINS_LIFTED')
            need(len(hits)>=1,'JOINS_NAKED_CLAUSE')
            sums.append(len(hits))
        out.append({'cut':Path(r['parent_certificate']).stem,'boundary_counts':sums})
    rebuilt=[]
    for name in ['S8232_reuse_02','S8232_reuse_03']:
        s,_,_=endpoint_ref('certificates/'+name+'.rank_path.json')
        rebuilt.append({'endpoint':'certificates/'+name+'.rank_path.json','changed_words':len(s^old),'changed_groups':sum(bool((s^old)&g['words']) for g in groups())})
    need(c.get('comparisons')==rebuilt,'JOINS_DIFFERENCES')
    return out

def check_rank_data(c,s):
    need(c.get('cardinality')==len(s),'RANK_DATA_CARDINALITY')
    L=c.get('L');r=c.get('rank')
    need(type(L) is int and type(r) is list and len(r)==len(V),'RANK_DATA_FORMAT')
    need(all(type(v) is int and 0<=v<=L for v in r),'RANK_DATA_DOMAIN')
    d=dict(zip(V,r))
    for w in W:
        if w not in s:need(d[w[1:]]>=d[w[:-1]]+1,'RANK_DATA_EDGE')
    p=c.get('path');walk(p);need(len(p)==L and not(set(p)&s),'RANK_DATA_PATH')
    dd,bb,pp=graph(s);need(r==[dd[v] for v in V] and len(pp)==L,'RANK_DATA_REPLAY')
    return L

def check_comparison(c,old):
    need(c.get('schema')=='R13-ORDINARY-AUGMENTATION-1','ORDINARY_SCHEMA')
    expected={(f'certificates/M107_panel_{i:02}.rank_path.json',w) for i in [2,3] for w in ['GAAAAAAA','TTTTTTTC']}
    need({(r['source'],r['added_word']) for r in c['items']}==expected and len(c['items'])==4,'ORDINARY_COMPLETE')
    for r in c['items']:
        m,_,_=endpoint_ref(r['source']);w=r['added_word'];need(w in WS and w not in m,'ORDINARY_NEW_WORD');n=m|{w}
        need(r.get('is_RC') is False and {rc(w) for w in n}!=n,'ORDINARY_NOT_RC')
        need(r.get('selected_sha256')==digest(n),'ORDINARY_HASH')
        need(check_rank_data(r['rank_path'],n)==106,'ORDINARY_LENGTH')
    need(old.get('schema')=='R13-ORIGINAL-COST-FOUR-1','COST4_SCHEMA')
    m,d,b=endpoint_ref(old['source']);need(max(d.values())==107 and len(m)==8230,'COST4_SOURCE')
    need(old.get('target_L')==106 and old.get('max_negative_additions')==3,'COST4_SCOPE')
    K={w for w in W if w not in m and d[w[:-1]]+1+b[w[1:]]==107};need(old.get('critical_words')==sorted(K),'COST4_SUPPORT')
    F=sorted(w for w in K if w==rc(w));P=sorted(w for w in K if w<rc(w))
    candidates=[frozenset(a) for j in range(4) for a in combinations(F,j)]
    for w in P:
        candidates.append(frozenset([w,rc(w)]));candidates.extend(frozenset([w,rc(w),f]) for f in F)
    need(len(candidates)==len(set(candidates))==old.get('candidate_count_including_empty'),'COST4_CENSUS')
    paths=[]
    for p in old['negative_paths']:
        walk(p);need(len(p)==107 and not(set(p)&m),'COST4_PATH');paths.append(frozenset(p))
    need(all(any(a.isdisjoint(p) for p in paths) for a in candidates),'COST4_COVERAGE')
    A=frozenset(old['attaining_additions']);need(len(A)==4 and A<=WS-m and {rc(w) for w in A}==A,'COST4_ADDITIONS')
    need(digest(m|A)==old.get('upper_sha256'),'COST4_UPPER_HASH');need(check_rank_data(old['upper_rank_path'],m|A)==105,'COST4_UPPER_LENGTH')
    return {'ordinary_singleton_witnesses':4,'ordinary_L':106,'ordinary_outputs_are_RC':False,'original_fixed_base_RC_cost':4,'original_negative_candidates':len(candidates),'original_cover_paths':len(paths)}

def check_global_cuts(c):
    need(c.get('schema')=='R13-GLOBAL-CYCLE-CUTS-1' and c.get('target_L')==106 and c.get('all_groups_free') is True,'CUT_SCHEMA_OR_SCOPE')
    cuts=c.get('cuts');need(type(cuts) is list and c.get('count')==len(cuts),'CUT_COUNT')
    seen=set();hist=Counter()
    for r in cuts:
        p=r['cycle'];walk(p,closed=True);q=[rc(w) for w in reversed(p)]
        key=min(tuple(v[i:]+v[:i]) for v in [p,q] for i in range(len(v)))
        need(tuple(p)==key and key not in seen,'CUT_CANONICAL_OR_DUPLICATE');seen.add(key)
        terms=Counter(min(w,rc(w)) for w in p)
        need(r.get('terms')==[[w,n] for w,n in sorted(terms.items())],'CUT_MULTIPLICITIES')
        need(r.get('lower')==(len(p)+106)//107,'CUT_CEILING');hist[r['lower']]+=1
    return {'globally_valid_RC_cycle_cuts':len(cuts),'lower_bound_histogram':dict(hist),'global_solve_claimed':False}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    report={'status':'RUNNING'}
    endpoints=[]
    for p in sorted((ROOT/'certificates').glob('*.rank_path.json')):
        c=load(p);s,d,b=check_endpoint(c);endpoints.append({'name':c['name'],'cardinality':len(s),'L':max(d.values()),'sha256':digest(s)})
    report['endpoints']=endpoints;report['augmentation_panel']=check_panel(load(ROOT/'certificates/augmentation_panel.json.gz'))
    report['bridge']=check_bridge(load(ROOT/'certificates/retiming_bridge.json'));report['affine_pruning']=check_affine(load(ROOT/'certificates/affine_pruning.json'));report['escape_joins']=check_joins(load(ROOT/'certificates/escape_joins.json'))
    report['comparison']=check_comparison(load(ROOT/'certificates/ordinary_augmentation.json'),load(ROOT/'certificates/original_augmentation_cost4.json.gz'))
    report['cycle_cuts']=check_global_cuts(load(ROOT/'certificates/new_cycle_cuts.json.gz'))
    report['status']='PASS';save(a.out/'certificates.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
