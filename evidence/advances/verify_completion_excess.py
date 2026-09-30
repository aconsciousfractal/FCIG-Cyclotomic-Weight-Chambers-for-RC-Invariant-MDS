"""Separate finite check of the R9/R11 fixed-exterior, excess-at-most-one corollary.
No submitted checker is imported; exact graph, exhaustive truth, 2-SAT via SCC.
"""
from collections import deque, defaultdict
from itertools import product
from pathlib import Path
import argparse,json

def need(x,m):
    if not x: raise ValueError(m)

def carrier(q,k):
    V=q**(k-1);E=q**k
    def rc(e,length=k):
        r=0
        for _ in range(length):r=q*r+q-1-e%q;e//=q
        return r
    orbits=[];owner={}
    for e in range(E):
        if e in owner:continue
        p=[];f=e
        while f not in p:p.append(f);f=(f%V)*q+f//V
        idx=len(orbits)
        for f in p:owner[f]=idx
        orbits.append(set(p))
    stable=[];groups=[]
    for i,p in enumerate(orbits):
        j=owner[rc(next(iter(p)))]
        if i==j:
            fixed=sorted(e for e in p if rc(e)==e);need(len(fixed)==2,'TWO_FIXED')
            stable.append((p,fixed))
        elif i<j:groups.append([frozenset([e,rc(e)]) for e in sorted(p)])
    rco=[frozenset([e,rc(e)]) for e in range(E) if e<=rc(e)]
    return V,E,rc,orbits,stable,groups,rco

def dag(V,edges):
    adj=[[] for _ in range(V)];deg=[0]*V
    for e,q in edges:adj[e].append(q);deg[q]+=1
    Q=deque(v for v in range(V) if not deg[v]);dist=[0]*V;seen=0
    while Q:
        u=Q.popleft();seen+=1
        for v in adj[u]:
            dist[v]=max(dist[v],dist[u]+1);deg[v]-=1
            if not deg[v]:Q.append(v)
    return (max(dist) if seen==V else None),dist,adj

def sat(n,clauses):
    # literal 2*i means restored; its xor-one negation means not restored.
    adj=[[] for _ in range(2*n)];rev=[[] for _ in range(2*n)]
    for a,b in clauses:
        adj[a^1].append(b);rev[b].append(a^1);adj[b^1].append(a);rev[a].append(b^1)
    seen=set();order=[]
    def dfs(v):
        seen.add(v)
        for w in adj[v]:
            if w not in seen:dfs(w)
        order.append(v)
    for v in range(2*n):
        if v not in seen:dfs(v)
    comp={}
    for v in reversed(order):
        if v in comp:continue
        todo=[v];comp[v]=v
        while todo:
            u=todo.pop()
            for w in rev[u]:
                if w not in comp:comp[w]=v;todo.append(w)
    return all(comp[2*i]!=comp[2*i+1] for i in range(n))

def run(q,k):
    V,E,rc,pcrs,stable,groups,rco=carrier(q,k);N=len(pcrs)
    fixed=sorted(e for p,fs in stable for e in fs);fi={e:i for i,e in enumerate(fixed)}
    nonstable=set(range(E))-set().union(*(p for p,fs in stable))
    truth=defaultdict(list);masks=0;all_dags=[]
    for bits in product((0,1),repeat=len(rco)):
        selected=set().union(*(p for p,b in zip(rco,bits) if b));masks+=1
        L,_,_=dag(V,[(e//q,e%V) for e in range(E) if e not in selected])
        if L is not None:all_dags.append((selected,L))
        if not N<=len(selected)<=N+1 or not all(p&selected for p in pcrs):continue
        key=frozenset(selected&nonstable)
        if L is not None:truth[key].append(L)
    branches=[(None,set())]
    for i,(p,fs) in enumerate(stable):
        branches.append((i,set(fs)))
        branches.extend((i,{e,rc(e)}) for e in sorted(p) if e<rc(e))
    stats=defaultdict(int);comparisons=0
    for options in product(*groups):
        external=set().union(*options);key=frozenset(external)
        for ell in range(V):
            success=False
            for exception,selection in branches:
                D=set(range(E))-external-selection-set(fixed)
                L,dist,adj=dag(V,[(e//q,e%V) for e in D])
                if L is None:stats['cyclic_D']+=1;continue
                if L>ell:stats['long_D']+=1;continue
                clauses=[]
                for i,(p,fs) in enumerate(stable):
                    a,b=[2*fi[e] for e in fs]
                    if i==exception:
                        for e in fs:
                            literal=2*fi[e]+int(e in selection);clauses.append((literal,literal))
                    else:clauses.extend([(a,b),(a^1,b^1)])
                reach=[]
                for u in range(V):
                    r={u};todo=[u]
                    while todo:
                        for v in adj[todo.pop()]:
                            if v not in r:r.add(v);todo.append(v)
                    reach.append(r)
                for i,e in enumerate(fixed):
                    u=e//q;v=e%V
                    if u in reach[v] or 2*dist[u]+1>ell:clauses.append((2*i+1,2*i+1))
                    for j,f in enumerate(fixed):
                        if j>i and f//q in reach[v]:clauses.append((2*i+1,2*j+1))
                feasible=sat(len(fixed),clauses);stats['sat' if feasible else 'unsat']+=1
                success|=feasible
            expected=any(L<=ell for L in truth[key]);need(success==expected,f'COMPLETION:{q}:{k}:{sorted(key)}:{ell}')
            comparisons+=1
    need(stats['sat']>0 and stats['unsat']>0,'BOTH_SAT_BRANCHES_EXERCISED')
    liftchecks=0;unlifted_failures=0;example=None
    for h in range(1,min(3,(V-1)//2)+1):
        ending=defaultdict(list)
        def paths(u,visited,edges):
            if len(edges)==h:ending[u].append(set(edges));return
            for a in range(q):
                e=q*u+a;v=e%V
                if e==rc(e) or v in visited:continue
                paths(v,visited|{v},edges+[e])
        for u in range(V):paths(u,{u},[])
        for p,fs in stable:
            for left,right in product(ending[fs[0]//q],ending[fs[1]//q]):
                U=left|right
                for selected,L in all_dags:
                    if L>2*h:continue
                    lhs=len(selected&U)+len(selected&p)-1
                    need(lhs>=1,'LIFTED_TWO_ARM');liftchecks+=1
                    if not selected&U:
                        unlifted_failures+=1
                        if example is None:example=dict(selected=sorted(selected),stable_pcr=sorted(p),left=sorted(left),right=sorted(right),target_L=2*h,actual_L=L)
    return dict(q=q,k=k,rc_masks=masks,branches_per_exterior=len(branches),comparisons=comparisons,branch_outcomes=dict(stats),lifted_two_arm_checks=liftchecks,unlifted_failures=unlifted_failures,first_unlifted_counterexample=example)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    rows=[run(2,2),run(2,4),run(4,2)]
    data={'status':'PASS','scope':'Exhaustive small carriers, fixed nonstable choices, excess at most one; no DNA k8 global search','rows':rows}
    a.out.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n');print(json.dumps(data))
if __name__=='__main__':main()
