"""Complete single-RC-PCR-group neighborhood, checked without any search code."""
from collections import Counter
import json
import graph_checker as g


def all_domains(base):
    seen=set();orbits=[];lookup={}
    for w in g.EDGES:
        if w in seen:continue
        orbit=sorted({w[i:]+w[:i] for i in range(g.K)})
        for z in orbit:lookup[z]=len(orbits)
        seen.update(orbit);orbits.append(orbit)
    used=set();groups=[]
    for j,orbit in enumerate(orbits):
        if j in used:continue
        mate=lookup[g.rc(orbit[0])];used.update((j,mate));options=[]
        for w in orbit:
            if mate==j and g.rc(w)!=w:continue
            options.append(sorted({w,g.rc(w)}))
        current=sorted(set(orbit+orbits[mate])&base)
        g.require(current in options,'RADIUS_INVALID_BASE_SELECTION')
        groups.append({'necklaces':[j,mate],'options':options,'base_option':options.index(current)})
    return groups


def verify_radius(base,cert,baseline_L):
    g.require(cert['schema']=='rc-mds-radius-one-exclusion-v1','RADIUS_SCHEMA')
    g.require(g.digest_set(base)==cert['baseline_sha256'],'RADIUS_BASE_DIGEST')
    B=cert['bound'];g.require(type(B) is int and B==baseline_L-1,'RADIUS_BOUND')
    groups=all_domains(base)
    g.require(json.dumps(groups,sort_keys=True)==json.dumps(cert['groups'],sort_keys=True),'RADIUS_COMPLETE_DOMAINS')
    counts=[];hits=[];rhss=[]
    for c in cert['cuts']:
        kind=c['kind'];g.require(kind in ('cycle','path'),'RADIUS_CUT_KIND')
        _,_,h=g.walk_check(c['words'],base,cycle=kind=='cycle')
        if kind=='path':g.require(len(c['words'])==B+1,'RADIUS_PATH_LENGTH');rhs=1
        else:rhs=(len(c['words'])+B)//(B+1)
        g.require(type(c['rhs']) is int and c['rhs']==rhs,'RADIUS_CUT_RHS')
        counts.append(Counter(c['words']));hits.append(h);rhss.append(rhs)
    g.require(any(h<r for h,r in zip(hits,rhss)),'RADIUS_BASE_NOT_EXCLUDED')
    wanted={(i,j) for i,gr in enumerate(groups) for j in range(len(gr['options'])) if j!=gr['base_option']}
    seen=set()
    for row in cert['reasons']:
        i,j,c=row['group'],row['option'],row['cut']
        g.require(all(type(z) is int for z in (i,j,c)),'RADIUS_REASON_TYPE')
        g.require((i,j) in wanted and (i,j) not in seen,'RADIUS_REASON_COVERAGE')
        g.require(0<=c<len(counts),'RADIUS_CUT_INDEX')
        gr=groups[i];old=gr['options'][gr['base_option']];new=gr['options'][j]
        updated=hits[c]-sum(counts[c][w] for w in old)+sum(counts[c][w] for w in new)
        g.require(updated<rhss[c],'RADIUS_INVALID_REJECTION')
        seen.add((i,j))
    g.require(seen==wanted,'RADIUS_MISSING_ALTERNATIVES')
    g.require(cert['tested_alternatives']==len(wanted) and cert['successes']==[],'RADIUS_CONCLUSION')
    return {'status':'EXACT_RADIUS_ONE_EXCLUSION','bound':B,'RC_PCR_groups':len(groups),
       'alternative_assignments':len(wanted),'checked_cuts':len(counts),'minimum_in_neighborhood':baseline_L,
       'unrestricted_global_refutation':False}
