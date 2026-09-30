"""Check global cuts from literal walks and, when used, 2-SAT contradictions.
No optimization solver, numerical tolerance, or discovery import.
"""
from verify_round9 import *
from collections import Counter

def path_between(path,s,t):
    if not path:req(s==t,'EMPTY_PATH_ENDPOINTS');return
    walk_check(path);req(path[0][:-1]==s and path[-1][1:]==t,'PATH_ENDPOINTS')

def validate_cut(c,selectable,stable_pairs):
    words={decode(i) for i in c['cut_words']};req(all(rc(w)!=w for w in words),'NONFIXED_CUT')
    kind=c['type']
    if kind in ('residual_cycle','long_fixed_path'):
        w=[decode(i) for i in c['walk']];walk_check(w,kind=='residual_cycle');req(words==set(w),'CUT_WALK_SUPPORT')
        if kind=='long_fixed_path':req(len(w)==107,'FORBIDDEN_PATH_LENGTH')
    elif kind=='stable_2sat_conflict':
        ew=[decode(i) for i in c['fixed_edges']];req(len(ew)==len(set(ew))==256 and all(rc(w)==w for w in ew),'SAT_FIXED_ARCS')
        req({frozenset(ew[i:i+2]) for i in range(0,256,2)}==stable_pairs,'SAT_PAIRS')
        valid=set();underlying=set()
        for row in c['implications']:
            i,j=row['from'],row['to'];r=row['reason'];p=[decode(k) for k in row['path']];underlying.update(p)
            req(all(rc(w)!=w for w in p),'SAT_NONFIXED_WITNESS')
            if r[0]=='length':
                a=r[1];req((i,j)==(a,a^1),'LENGTH_ARROW');req(bool(p),'LENGTH_PREFIX');walk_check(p);req(p[-1][1:]==ew[a][:-1] and 2*len(p)+1>106,'HALF_DEPTH_OBSTRUCTION')
            elif r[0]=='fixed-return':
                a=r[1];req((i,j)==(a,a^1),'RETURN_ARROW');path_between(p,ew[a][1:],ew[a][:-1])
            elif r[0]=='pair':
                a,b=r[1:];req((i,j) in [(a,b^1),(b,a^1)],'PAIR_ARROW');path_between(p,ew[a][1:],ew[b][:-1])
            else:raise ValueError('SAT_REASON')
            valid.add((i,j))
        proof=c['contradiction'];a=2*proof['variable'];f=proof['forward'];b=proof['backward'];req(f[0]==a and f[-1]==(a^1) and b[0]==(a^1) and b[-1]==a,'SAT_CONTRADICTION_ENDPOINTS')
        req(all(pair in valid for p in [f,b] for pair in zip(p,p[1:])),'SAT_CONTRADICTION_ARROWS');req(words==underlying,'SAT_CUT_SUPPORT')
    else:raise ValueError('CUT_TYPE')
    return words&selectable

def main():
    t=time.time();base=read_words(ROOT/'inputs/M107.words.txt');stable=necklaces(base);stablewords=set().union(*stable);selectable=EDGESET-stablewords;pairs={frozenset(w for w in p if rc(w)==w) for p in stable}
    certs=json.loads((ROOT/'certificates/boolean_master_cuts.json').read_text());profile=Counter();iterations=set();rowsizes=[];points={}
    for c in certs:
        support=validate_cut(c,selectable,pairs);req(bool(support),'NONEMPTY_GLOBAL_CUT');profile[c['type']]+=1;rowsizes.append(len(support));iterations.add(c['iteration'])
        outer=(base-stablewords)-{decode(i) for i in c['outer_removed_selected']};outer|={decode(i) for i in c['outer_changed_selected']}
        req(len(outer)==8102 and {rc(w) for w in outer}==outer,'OUTER_CARDINALITY_RC')
        req(not (support&outer),'CUT_REJECTS_CURRENT_OUTER')
        # One selected word per nonstable PCR; stable choices not instantiated here.
        checked=set()
        for w in outer:
            p=frozenset(w[i:]+w[:i] for i in range(8));req(p not in checked,'OUTER_DUPLICATE_PCR');checked.add(p)
        point=outer|{w for w in EDGES if rc(w)==w};digest=hashlib.sha256(bytes(int(w in point) for w in EDGES)).hexdigest()
        if c['iteration'] in points:req(points[c['iteration']]==digest,'CUT_POINT_COHERENCE')
        points[c['iteration']]=digest
    calls=json.loads((ROOT/'logs/boolean_master_calls.json').read_text())
    for row in calls:req(points[row['iteration']+1]==row['point_sha256'],'SOLVER_INCUMBENT_JOIN')
    result={'status':'PASS','globally_valid_cuts':len(certs),'types':dict(profile),'outer_configurations_rejected_including_base':len(points),'distinct_outer_hashes':len(set(points.values())),'solver_incumbents_joined':len(calls),'cut_support_size_range':[min(rowsizes),max(rowsizes)],'seconds':time.time()-t,'boundary':'individual valid cuts and point rejections, not a global UNSAT certificate; 2-SAT cut branch is covered here only when it occurs in saved objects'}
    (OUT/'master_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
