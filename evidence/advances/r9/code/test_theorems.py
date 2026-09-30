"""Exact finite regressions, not the proof of the general statements."""
from itertools import product
from collections import Counter
from pathlib import Path
import json,time,os
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('R9_CHECK_OUTPUT',ROOT/'run_output'));OUT.mkdir(parents=True,exist_ok=True)

def dag(n,edges):
    adj=[[] for _ in range(n)];ind=[0]*n;a=[0]*n
    for u,v in edges:adj[u].append(v);ind[v]+=1
    q=[i for i in range(n) if not ind[i]]
    for u in q:
        for v in adj[u]:
            a[v]=max(a[v],a[u]+1);ind[v]-=1
            if not ind[v]:q.append(v)
    if len(q)!=n:return None
    b=[0]*n;reach=[1<<i for i in range(n)]
    for u in reversed(q):
        for v in adj[u]:b[u]=max(b[u],1+b[v]);reach[u]|=reach[v]
    return max(a),a,b,reach,q,adj

def fixed_arc_regression(n):
    mate=lambda u:u^1;allpairs=[];seen=set()
    for u in range(n):
        for v in range(n):
            if u==v or v==mate(u) or (u,v) in seen:continue
            pair={(u,v),(mate(v),mate(u))};seen|=pair;allpairs.append(pair)
    data={'n':n,'fixed_D_candidates':2**len(allpairs),'fixed_D_acyclic':0,'completions':0,'acyclic_completions':0}
    special=[(u,mate(u)) for u in range(n)]
    for mask in range(1<<len(allpairs)):
        edges=set().union(*(allpairs[i] for i in range(len(allpairs)) if mask>>i&1))
        dd=dag(n,edges)
        if dd is None:continue
        data['fixed_D_acyclic']+=1;ld,a,b,reach,*_=dd
        for rm in range(1<<n):
            R=[i for i in range(n) if rm>>i&1];blocked=any(reach[mate(i)]>>i&1 for i in R)
            blocked |= any(reach[mate(i)]>>j&1 for i in R for j in R if i!=j)
            full=dag(n,edges|{special[i] for i in R});data['completions']+=1
            if (full is None)!=blocked:raise ValueError('TWO_SAT_EQUIVALENCE')
            if full is not None:
                predicted=max([ld]+[a[i]+1+b[mate(i)] for i in R]);
                if full[0]!=predicted:raise ValueError('EXACT_STABLE_LENGTH')
                data['acyclic_completions']+=1
    return data

def feasible_difference(n,rows):
    d=[0]*n
    for k in range(n):
        changed=False
        for u,v,c in rows:
            if d[v]<d[u]+c:d[v]=d[u]+c;changed=True
        if not changed:return True
    return False

def projection_regression():
    n=5;chain=[(i,i+1) for i in range(4)];F=[(4,0),(1,3),(3,1),(0,4)];B=sorted({v for e in F for v in e});idx={v:i for i,v in enumerate(B)};count=0
    for dm in range(16):
        D={e for i,e in enumerate(chain) if dm>>i&1};ld,a,b,reach,q,adj=dag(n,D);fixed=[]
        for s in B:
            ds={s:0}
            for u in q:
                if u not in ds:continue
                if u!=s and u in idx:fixed.append((idx[s],idx[u],ds[u]));continue
                for v in adj[u]:ds[v]=max(ds.get(v,-1),ds[u]+1)
        for rm in range(16):
            R={e for i,e in enumerate(F) if rm>>i&1};full=dag(n,D|R)
            for L in range(7):
                origin=len(B);rows=list(fixed)+[(idx[u],idx[v],1) for u,v in R]
                for u in B:rows.extend([(origin,idx[u],a[u]),(idx[u],origin,b[u]-L)])
                proj=ld<=L and feasible_difference(len(B)+1,rows)
                truth=full is not None and full[0]<=L
                if proj!=truth:raise ValueError('BOUNDARY_EQUIVALENCE')
                count+=1
    # No boundary: the fixed longest path is the whole condition.
    if dag(3,{(0,1),(1,2)})[0]!=2:raise ValueError('EMPTY_BOUNDARY')
    # Missing RC hypothesis: two fixed arcs can coexist in this nonsymmetric DAG.
    if dag(4,{(0,1),(1,2),(2,3)}) is None:raise ValueError('MISSING_RC_CONTROL')
    # Missing prefix bounds gives a false acceptance, while full graph has length 3.
    if not feasible_difference(3,[(0,1,1),(2,0,0),(0,2,-2),(2,1,0),(1,2,-2)]):raise ValueError('PREFIX_CONTROL')
    if dag(4,{(0,1),(1,2),(2,3)})[0]!=3:raise ValueError('PREFIX_FULL')
    return {'comparisons':count,'omitted_boundary_vertex':2,'empty_boundary_checked':True,'missing_RC_counterexample_checked':True,'missing_prefix_counterexample_checked':True}

if __name__=='__main__':
    t=time.time();out={'status':'PASS','scope':'complete labelled fixed-D subfamilies on 4 and 6 paired vertices; all subsets of their fixed arcs','fixed_arc_tests':[fixed_arc_regression(4),fixed_arc_regression(6)],'projection':projection_regression(),'seconds':time.time()-t}
    (OUT/'theorem_regressions.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
