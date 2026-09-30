"""Exhaustive small-graph checks of critical supports and scoped cycle arguments.
These tests do not replace the all-size proofs in the report.
"""
from literal import need,save
from itertools import product,combinations
from collections import deque
from pathlib import Path
import argparse,random

def small():
    A='01';V=[''.join(v) for v in product(A,repeat=3)];E=[v+a for v in V for a in A]
    def rc(w):return ''.join('1' if x=='0' else '0' for x in w[::-1])
    def dag(s):
        adj={u:[] for u in V};indeg=dict.fromkeys(V,0)
        for w in E:
            if w not in s:adj[w[:-1]].append(w);indeg[w[1:]]+=1
        q=deque(v for v in V if indeg[v]==0);order=[];d=dict.fromkeys(V,0)
        while q:
            u=q.popleft();order.append(u)
            for w in adj[u]:
                v=w[1:];d[v]=max(d[v],d[u]+1);indeg[v]-=1
                if indeg[v]==0:q.append(v)
        if len(order)!=len(V):return None
        b=dict.fromkeys(V,0)
        for u in reversed(order):b[u]=max([0]+[b[w[1:]]+1 for w in adj[u]])
        return max(d.values()),d,b,adj
    orbits=[frozenset([w,rc(w)]) for w in E if w<=rc(w)]
    # Enumerate simple directed cycles, identifying rotations but not RC.
    cycles=set()
    def extend(start,u,verts,edges):
        for w in E:
            if w[:-1]!=u:continue
            v=w[1:]
            if v==start:
                p=edges+[w];cycles.add(min(tuple(p[i:]+p[:i]) for i in range(len(p))))
            elif v not in verts:extend(start,v,verts|{v},edges+[w])
    for v in V:extend(v,v,{v},[])
    rc_graphs=1<<len(orbits);dags=0;augtests=0;cycletests=0;supporttests=0
    for mask in range(rc_graphs):
        s=frozenset().union(*(o for j,o in enumerate(orbits) if (mask>>j)&1));data=dag(s)
        if data is None:continue
        dags+=1;L,d,b,adj=data
        paths=[]
        def visit(u,p):
            if len(p)==L:paths.append(frozenset(p));return
            for w in adj[u]:visit(w[1:],p+[w])
        for v in V:visit(v,[])
        K={w for w in E if w not in s and d[w[:-1]]+1+b[w[1:]]==L}
        need(K==frozenset().union(*paths),'SMALL_CRITICAL_SUPPORT');supporttests+=1
        need({rc(w) for w in K}==K,'SMALL_CRITICAL_RC')
        F=[w for w in E if w not in s and w==rc(w)];P=[w for w in E if w not in s and w<rc(w)]
        candidates=[frozenset()]+[frozenset([w]) for w in F]+[frozenset(p) for p in combinations(F,2)]+[frozenset([w,rc(w)]) for w in P]
        for a in candidates:
            final=dag(s|a)[0];reduced=dag(s|(a&K))[0]
            need((final<L)==all(bool(a&p) for p in paths),'SMALL_HIT_EQUIVALENCE')
            need((final<L)==(reduced<L),'SMALL_OUTSIDE_IRRELEVANCE');augtests+=1
        for Q in cycles:
            r=len(set(Q)&s);need(r>=1 and len(Q)<=(L+1)*r,'SMALL_CYCLE_BOUND');cycletests+=1
    # Canonical cochain normalization on a directed spanning tree.
    q=deque(['000']);seen={'000'};tree=[]
    while q:
        u=q.popleft()
        for a in A:
            w=u+a;v=w[1:]
            if v not in seen:seen.add(v);q.append(v);tree.append(w)
    need(seen==set(V),'SMALL_TREE_REACHES')
    def norm(x):
        f={'000':0}
        for w in tree:f[w[1:]]=f[w[:-1]]-x[w]
        return tuple(x[w]-f[w[:-1]]+f[w[1:]] for w in E)
    rng=random.Random(131313);gradienttests=0
    for _ in range(256):
        x={w:rng.randrange(-3,4) for w in E};f={v:rng.randrange(-3,4) for v in V}
        y={w:x[w]+f[w[:-1]]-f[w[1:]] for w in E};need(norm(x)==norm(y),'SMALL_GRADIENT_NORMALIZATION')
        for Q in cycles:need(sum(x[w] for w in Q)==sum(y[w] for w in Q),'SMALL_CYCLE_TELESCOPING')
        z=x.copy();z['0000']+=1;need(norm(x)!=norm(z),'SMALL_NONGRADIENT_DIFFERENCE');gradienttests+=1
    return {'status':'PASS','binary_order':4,'RC_masks':rc_graphs,'DAG_masks':dags,'critical_support_checks':supporttests,'augmentation_checks':augtests,'simple_cycles':len(cycles),'cycle_inequality_checks':cycletests,'gradient_samples':gradienttests,'scope':'Finite regressions; not premises replacing general proofs.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);r=small();save(a.out/'regressions.json',r);print(r)
if __name__=='__main__':main()
