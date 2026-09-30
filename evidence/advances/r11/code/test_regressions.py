#!/usr/bin/env python3
"""Finite regressions, not substitutes for the proofs in the report."""
from exact import ROOT, Path, need, save
from itertools import product
from collections import deque
import argparse


def depths(n, edges):
    a=[[] for _ in range(n)]; deg=[0]*n; d=[0]*n
    for u,v in edges:a[u].append(v);deg[v]+=1
    q=deque(i for i in range(n) if deg[i]==0); order=[]
    while q:
        u=q.popleft();order.append(u)
        for v in a[u]:
            d[v]=max(d[v],d[u]+1);deg[v]-=1
            if not deg[v]:q.append(v)
    if len(order)!=n:return None
    b=[0]*n
    for u in reversed(order):b[u]=max([0]+[b[v]+1 for v in a[u]])
    return d,b


def main(out):
    graph_count=acyclic_count=rank_cases=boundary_cases=rounding_cases=0
    for n in (3,4):
        theta=lambda v:n-1-v
        seen=set();orbits=[]
        for u in range(n):
            for v in range(n):
                if (u,v) not in seen:
                    o={(u,v),(theta(v),theta(u))};seen.update(o);orbits.append(o)
        for mask in range(1<<len(orbits)):
            graph_count+=1
            selected=[orbits[j] for j in range(len(orbits)) if mask>>j&1]
            edges=set().union(*selected) if selected else set()
            data=depths(n,edges)
            if data is None:continue
            acyclic_count+=1;d,b=data
            for B in range(5):
                if max(d)>B:continue
                rank_cases+=1;r2=[B+d[v]-d[theta(v)] for v in range(n)]
                need(all(0<=x<=2*B for x in r2),'SMALL_RANK_DOMAIN')
                need(all(r2[v]+r2[theta(v)]==2*B for v in range(n)),'SMALL_RANK_RC')
                need(all(r2[v]-r2[u]>=2 for u,v in edges),'SMALL_RANK_EDGE')
                for sm in range(1<<len(selected)):
                    small=set().union(*(selected[j] for j in range(len(selected)) if sm>>j&1))
                    d0,b0=depths(n,small)
                    need(all(2*d0[v]<=r2[v]<=2*(B-b0[v]) for v in range(n)),'SMALL_BOUNDARY')
                    boundary_cases+=1
                if B%2==0:
                    H=B//2;p=[d[v]-d[theta(v)] for v in range(n)]
                    h=[((x+1)//2 if x>0 else -((-x+1)//2) if x<0 else 0) for x in p]
                    need(all(-H<=x<=H for x in h),'SMALL_INTEGER_HEIGHT_DOMAIN')
                    need(all(h[v]==-h[theta(v)] for v in range(n)),'SMALL_INTEGER_HEIGHT_RC')
                    need(all(h[v]>=h[u]+1 for u,v in edges),'SMALL_INTEGER_HEIGHT_EDGE')
                    rounding_cases+=1
    truth=0
    for hu,hv,gu,gv in product(range(-2,3),repeat=4):
        f00=int(hv<=hu);f10=int(hv<=gu);f01=int(gv<=hu);f11=int(gv<=gu)
        for a,b in product((0,1),repeat=2):
            direct=int((gv if b else hv)<=(gu if a else hu))
            expr=f00+(f10-f00)*a+(f01-f00)*b+(f11-f01-f10+f00)*a*b
            need(expr==direct,'FUSION_TRUTH_TABLE');truth+=1
    surplus_cases=0
    # Two nonstable pairs and three stable PCRs, every count between one and three.
    for a,b,c,d,e in product(range(1,4),repeat=5):
        eps=2*(a-1)+2*(b-1)+(c-1)+(d-1)+(e-1)
        expected=(a==b==1 and (c-1)+(d-1)+(e-1)<=1)
        need((eps<=1)==expected,'SMALL_SURPLUS_FORM')
        if eps==1:need(sorted((c,d,e))==[1,1,2],'SMALL_ODD_SURPLUS_LOCATION')
        surplus_cases+=1
    result={'status':'PASS','RC_digraphs_including_loops':graph_count,'acyclic_RC_digraphs':acyclic_count,
            'symmetric_real_rank_cases':rank_cases,'fixed_subgraph_bound_cases':boundary_cases,
            'even_integer_roundings':rounding_cases,'binary_fusion_truth_cases':truth,
            'small_surplus_cases':surplus_cases,
            'scope':'All involution-closed digraphs on n=3,4 with theta(v)=n-1-v; B=0..4 when feasible. Fusion labels -2..2. Surplus test is a finite abstract model.'}
    save(out/'regressions.json',result);print(result)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);main(a.out)
