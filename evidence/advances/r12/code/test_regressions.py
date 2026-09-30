"""Small exhaustive and seeded finite regressions, not substitutes for proofs."""
from pathlib import Path
from collections import deque,Counter
from itertools import product
import argparse,json,random
from literal import ROOT,need,save


def badwalk(n,edges,B):
    adj=[[] for _ in range(n)];ind=[0]*n
    for e in edges:adj[e[0]].append(e);ind[e[1]]+=1
    q=deque(i for i in range(n) if not ind[i]);order=[];d=[0]*n;par={}
    while q:
        u=q.popleft();order.append(u)
        for e in adj[u]:
            v=e[1]
            if d[v]<d[u]+1:d[v]=d[u]+1;par[v]=e
            ind[v]-=1
            if not ind[v]:q.append(v)
    if len(order)<n:
        color=[0]*n;stack=[];position={}
        def visit(u):
            color[u]=1;position[u]=len(stack)
            for e in adj[u]:
                v=e[1]
                if color[v]==1:return stack[position[v]:]+[e]
                if not color[v]:
                    stack.append(e);ans=visit(v)
                    if ans:return ans
                    stack.pop()
            color[u]=2;position.pop(u);return None
        for u in range(n):
            if not color[u]:
                z=visit(u)
                if z:return ('cycle',z,d)
        raise ValueError('CYCLE_EXTRACTION')
    if max(d,default=0)>B:
        u=max(range(n),key=lambda i:d[i]);p=[]
        while u in par:e=par[u];p.append(e);u=e[0]
        return ('path',list(reversed(p))[:B+1],d)
    return ('ok',[],d)


def tiny_DNA():
    total=acyclic=near=drop_checks=rank_checks=0;records=[]
    for q,k in [(2,2),(2,4),(4,2)]:
        ws=list(product(range(q),repeat=k));vs=list(product(range(q),repeat=k-1));vi={v:i for i,v in enumerate(vs)}
        rc=lambda w:tuple(q-1-a for a in reversed(w))
        orbits=[];seen=set()
        for w in ws:
            if w not in seen:o={w,rc(w)};orbits.append(o);seen|=o
        pcr=[];seen=set()
        for w in ws:
            if w not in seen:p={w[i:]+w[:i] for i in range(k)};pcr.append(p);seen|=p
        pcowner={w:i for i,p in enumerate(pcr) for w in p};local=Counter()
        for bits in product((0,1),repeat=len(orbits)):
            s=set().union(*(o for b,o in zip(bits,orbits) if b));total+=1;local['sets']+=1
            es={(vi[w[:-1]],vi[w[1:]]) for w in ws if w not in s}
            kind,_,d=badwalk(len(vs),es,len(vs))
            if kind!='ok':continue
            acyclic+=1;local['acyclic']+=1
            c=[len(p&s) for p in pcr];need(min(c)>=1,'TINY_PCR_HIT')
            stable=[i for i,p in enumerate(pcr) if {rc(w) for w in p}==p]
            paired=[i for i,p in enumerate(pcr) if i<pcowner[rc(next(iter(p)))]]
            need(len(s)-len(pcr)==2*sum(c[i]-1 for i in paired)+sum(c[i]-1 for i in stable),'TINY_SURPLUS')
            constraint=(all(c[i]==1 for i in range(len(pcr)) if i not in stable)
                        and all(1<=c[i]<=2 for i in stable) and sum(c[i]-1 for i in stable)<=1)
            need((len(s)<=len(pcr)+1)==constraint,'TINY_NEAR_EQUIVALENCE')
            if constraint:near+=1
            B=max(d);r2=[B+d[vi[v]]-d[vi[rc(v)]] for v in vs]
            need(all(0<=r<=2*B for r in r2),'TINY_RANK_BOUNDS')
            for w in ws:
                delta=r2[vi[w[1:]]]-r2[vi[w[:-1]]]
                if w not in s:need(delta>=2,'TINY_RANK_EDGE');rank_checks+=1
                elif c[pcowner[w]]==1:
                    need(delta<=-2*(len(pcr[pcowner[w]])-1),'TINY_SINGLE_DROP');drop_checks+=1
        records.append({'q':q,'k':k,**local})
    return {'RC_sets':total,'acyclic_sets':acyclic,'near_cardinality_sets':near,'single_cut_drop_checks':drop_checks,'rank_edge_checks':rank_checks,'records':records}


def generic_filter():
    rng=random.Random(12012);cases=80;removed=refuted=solutions_checked=lift_checks=0
    for t in range(cases):
        n=4;es=rng.sample([(i,j) for i in range(n) for j in range(n)],8)
        gg=[set(es[2*i:2*i+2]) for i in range(4)];owner={e:i for i,g in enumerate(gg) for e in g}
        allopts=[tuple(sorted(g)) for g in gg];frees=set(rng.sample(range(4),rng.randrange(1,5)))
        frozen={i:rng.choice(allopts[i]) for i in range(4) if i not in frees};B=rng.randrange(1,4)
        remaining={i:set(allopts[i]) for i in frees};basecut=set(frozen.values());trace=[];boundary=set();used=set()
        brute=[]
        for choices in product(*(allopts[i] for i in sorted(frees))):
            sel=basecut|set(choices)
            if badwalk(n,set(es)-sel,B)[0]=='ok':brute.append(dict(zip(sorted(frees),choices)))
            solutions_checked+=1
        changed=True
        while changed and all(remaining.values()):
            changed=False
            mandatory={e for e in es if owner[e] not in frees and e not in basecut}
            mandatory|={e for e in es if owner[e] in frees and e not in remaining[owner[e]]}
            for g in sorted(frees):
                for o in sorted(remaining[g]):
                    trial=mandatory|(gg[g]-{o});kind,path,_=badwalk(n,trial,B)
                    if kind=='ok':continue
                    # Compare the logically justified removal against every genuine completion.
                    need(all(sol[g]!=o for sol in brute),'FILTER_REMOVED_SOLUTION')
                    remaining[g].remove(o);trace.append((g,o,kind,path));used.add(g)
                    boundary|={e for e in path if owner[e] not in frees};removed+=1;changed=True
                    mandatory.add(o)
                    if not remaining[g]:break
                if not remaining[g]:break
        if all(remaining.values()):continue
        refuted+=1;need(not brute,'FALSE_GENERIC_REFUTATION')
        # Check the lifted clause over arbitrary nonempty selections in every
        # group, not merely one-per-group endpoints.
        opts=[[{p[0]},{p[1]},set(p)] for p in allopts]
        for choices in product(*opts):
            sel=set().union(*choices)
            if badwalk(n,set(es)-sel,B)[0]!='ok':continue
            lhs=len(sel&boundary)+sum(len(choices[g])-1 for g in used)
            need(lhs>=1,'LIFTED_CLAUSE_COUNTEREXAMPLE');lift_checks+=1
    return {'seed':12012,'instances':cases,'certified_removals':removed,'empty_domains':refuted,
            'transversal_assignments_tested':solutions_checked,'acyclic_lifted_assignments_tested':lift_checks,
            'not_asserted':'Completeness of domain propagation on unrefuted instances.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    result={'status':'PASS','tiny_DNA':tiny_DNA(),'generic_filter':generic_filter()}
    save(a.out/'regressions.json',result);print(json.dumps(result,indent=2))
