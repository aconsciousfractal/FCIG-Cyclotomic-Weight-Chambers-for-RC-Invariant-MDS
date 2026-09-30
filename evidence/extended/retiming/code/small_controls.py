"""Finite exact controls, never used as the general proof.
Exhausts all subsets in three small de Bruijn carriers and all 4-vertex
loopless RC-anti-invariant digraphs. Original large-search code is not used.
"""
from __future__ import annotations
from itertools import product
from fractions import Fraction
from collections import deque
import retime as rt

def all_cycles(q,k):
    n=q**(k-1);cycles=[]
    for s in range(n):
        def visit(u,vertices,edges):
            for e in range(q*u,q*u+q):
                v=e%n
                if v==s:cycles.append(edges+[e])
                elif v>s and v not in vertices:visit(v,vertices|{v},edges+[e])
        visit(s,{s},[])
    return cycles

def census(q,k):
    n=q**(k-1);E=q*n; pcr=[];seen=set()
    for x in range(E):
        if x in seen:continue
        digits=[];t=x
        for _ in range(k):digits.append(t%q);t//=q
        digits=digits[::-1];orbit=set()
        for j in range(k):
            y=0
            for a in digits[j:]+digits[:j]:y=q*y+a
            orbit.add(y)
        pcr.append(orbit);seen.update(orbit)
    N=len(pcr);sets={}
    for mask in range(1<<E):
        if mask.bit_count()!=N:continue
        if any(sum((mask>>e)&1 for e in C)!=1 for C in pcr):continue
        m=[(mask>>e)&1 for e in range(E)]
        try:z=rt.dag(m,q,k)
        except rt.Failure:continue
        sets[mask]=(m,z['L'])
    # Actual F/RF adjacency: no cut-normal classification is assumed.
    adj={mask:set() for mask in sets}
    for mask,(m,L) in sets.items():
        for u in range(n):
            inc={a*n+u for a in range(q)};out=set(range(q*u,q*u+q))
            if all(m[e] for e in inc):
                new=(mask & ~sum(1<<e for e in inc)) | sum(1<<e for e in out)
                rt.require(new in sets,'SMALL_F_CLOSURE');adj[mask].add(new);adj[new].add(mask)
    groups=[];seen=set()
    for mask in sets:
        if mask in seen:continue
        group=[];queue=[mask];seen.add(mask)
        while queue:
            x=queue.pop();group.append(x)
            for y in adj[x]:
                if y not in seen:seen.add(y);queue.append(y)
        groups.append(group)
    cycles=all_cycles(q,k);reports=[];tests=0
    for group in groups:
        src=group[0];m=sets[src][0]
        ratio=max(Fraction(len(c),sum(m[e] for e in c)) for c in cycles)
        T=(ratio.numerator+ratio.denominator-1)//ratio.denominator
        rt.require(min(sets[x][1] for x in group)==T-1,'SMALL_CYCLE_RATIO_FORMULA')
        rcsets=[x for x in group if all(sets[x][0][e]==sets[x][0][rt.theta(e,k,q)] for e in range(E))]
        rep={'component_size':len(group),'ordinary_min':T-1,'ratio':[ratio.numerator,ratio.denominator],'rc_count':len(rcsets)}
        if rcsets:
            m=sets[rcsets[0]][0];sectors={0:[],1:[]};th0=rt.theta(0,k-1,q)
            for mask in rcsets:
                b=sets[mask][0];f=[0]*n
                for v in range(1,n):f[v]=f[v//q]+m[v]-b[v]
                rt.require(all(b[e]-m[e]==f[e//q]-f[e%n] for e in range(E)),'SMALL_CUT')
                C=f[0]+f[th0];rt.require(all(f[v]+f[rt.theta(v,k-1,q)]==C for v in range(n)),'SMALL_PARITY_CONSTANT')
                sectors[C%2].append(mask)
            rt.require(all(sectors.values()),'TWO_SECTORS_NONEMPTY')
            rep['sector_minima']=[min(sets[x][1] for x in sectors[ep]) for ep in (0,1)]
            for horizon in range(max(1,T-1),T+3):
                ordinary=rt.bellman_ford(m,q,k,horizon)
                for ep in (0,1):
                    expected=any(sets[x][1]<horizon for x in sectors[ep])
                    if ordinary['status']=='FEASIBLE':
                        answer=rt.rounded_rc(m,q,k,horizon,ordinary['potential'],ep)
                        actual=answer['status']=='FEASIBLE_RC'
                        if actual:
                            mask=sum(x<<e for e,x in enumerate(answer['indicator']))
                            rt.require(mask in sectors[ep] and sets[mask][1]<horizon,'SMALL_MODEL_REALIZATION')
                    else:actual=False
                    rt.require(actual==expected,'SMALL_RC_FEASIBILITY');tests+=1
        reports.append(rep)
    return {'q':q,'k':k,'subsets_scanned':1<<E,'ordinary_mds':len(sets),'components':len(groups),'directed_simple_cycles':len(cycles),'rc_feasibility_comparisons':tests,'rows':reports}

def symmetric_rank_controls():
    n=4;th=[3,2,1,0];all_edges=[(u,v) for u in range(n) for v in range(n) if u!=v];orbits=[];seen=set()
    for e in all_edges:
        if e in seen:continue
        pair={e,(th[e[1]],th[e[0]])};orbits.append(pair);seen|=pair
    dags=0;checks=0
    for bits in range(1<<len(orbits)):
        es=set().union(*(orbits[i] for i in range(len(orbits)) if bits>>i&1)) if bits else set()
        indeg=[0]*n;d=[0]*n;adj=[[] for _ in range(n)]
        for u,v in es:indeg[v]+=1;adj[u].append(v)
        queue=deque(u for u in range(n) if indeg[u]==0);num=0
        while queue:
            u=queue.popleft();num+=1
            for v in adj[u]:
                d[v]=max(d[v],d[u]+1);indeg[v]-=1
                if indeg[v]==0:queue.append(v)
        if num<n:continue
        dags+=1
        for h0 in range(max(d),max(d)+4):
            r2=[d[v]+h0-d[th[v]] for v in range(n)];h=[]
            for v in range(n):
                if r2[v]%2==0:h.append(r2[v]//2)
                elif r2[v]<h0:h.append(r2[v]//2)
                elif r2[v]>h0:h.append((r2[v]+1)//2)
                else:h.append((r2[v]+(1 if v>th[v] else -1))//2)
            rt.require(all(0<=h[v]<=h0 and h[v]+h[th[v]]==h0 for v in range(n)),'SELF_DUAL_RANK_SUM')
            rt.require(all(h[v]-h[u]>=1 for u,v in es),'SELF_DUAL_RANK_EDGES');checks+=1
    return {'anti_invariant_graphs':1<<len(orbits),'dags':dags,'integer_rank_checks':checks}

def run():
    return {'status':'PASS','scope':'finite exact controls only','censuses':[census(2,2),census(2,4),census(4,2)],'self_dual_ranks':symmetric_rank_controls()}
if __name__=='__main__':
    import json
    print(json.dumps(run(),sort_keys=True,indent=2))
