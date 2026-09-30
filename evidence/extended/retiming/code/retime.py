"""Exact unit-retiming and RC rounding. Standard library; no numerical solver.

Conventions: q-ary words in lexicographic base q, edges length k, lower
vertices length k-1. Complement is digit a -> q-1-a. The feasibility
algorithm uses integer Bellman--Ford. Large discovery runs used the separate
C++ queue implementation; neither implementation is trusted by the checker.
"""
from __future__ import annotations
from collections import deque
from pathlib import Path
import argparse, json

class Failure(ValueError):
    pass

def require(ok: bool, message: str) -> None:
    if not ok:
        raise Failure(message)

def theta(x: int, length: int, q: int) -> int:
    y = 0
    for _ in range(length):
        y = q*y + q-1-x%q
        x //= q
    return y

def read_bits(path: Path, size: int) -> list[int]:
    b = path.read_bytes()
    require(len(b) == (size+7)//8, 'BIT_LENGTH')
    require(size%8 == 0 or b[-1] >> (size%8) == 0, 'NONZERO_PADDING')
    return [(b[e//8] >> (e%8)) & 1 for e in range(size)]

def write_bits(path: Path, m: list[int]) -> None:
    b = bytearray((len(m)+7)//8)
    for e,x in enumerate(m):
        require(x in (0,1), 'NONBOOLEAN')
        b[e//8] |= x << (e%8)
    path.write_bytes(bytes(b))

def pcr_check(m: list[int], q: int, k: int) -> int:
    n = q**(k-1); seen = bytearray(q*n); count = 0
    require(len(m) == q*n and all(x in (0,1) for x in m), 'INPUT_INDICATOR')
    for e in range(q*n):
        if seen[e]:
            continue
        x=e; total=0
        while not seen[x]:
            seen[x]=1; total+=m[x]; x=(x%n)*q+x//n
        require(x==e and total==1, 'PCR_NOT_EXACTLY_ONE')
        count+=1
    return count

def dag(m: list[int], q: int, k: int) -> dict:
    n=q**(k-1); incoming=[0]*n; level=[0]*n; parent=[-1]*n
    for e,x in enumerate(m):
        if not x: incoming[e%n]+=1
    queue=deque(v for v in range(n) if incoming[v]==0); scanned=0
    while queue:
        u=queue.popleft(); scanned+=1
        for e in range(q*u,q*u+q):
            if m[e]: continue
            v=e%n
            if level[v]<level[u]+1:
                level[v]=level[u]+1; parent[v]=e
            incoming[v]-=1
            if incoming[v]==0: queue.append(v)
    require(scanned==n, 'RESIDUAL_CYCLE')
    end=max(range(n),key=level.__getitem__); path=[];v=end
    while parent[v]>=0:
        e=parent[v];path.append(e);v=e//q
    path.reverse()
    return {'L':level[end], 'levels':level, 'path_edges':path}

def bellman_ford(m: list[int], q: int, k: int, T: int) -> dict:
    """O(|V||E|) exact feasibility, with an independently checkable output."""
    require(isinstance(T,int) and T>=1, 'HORIZON')
    n=q**(k-1); d=[0]*n; parent=[-1]*n; changed=-1
    for _ in range(n):
        changed=-1
        for e,x in enumerate(m):
            u=e//q; v=e%n; value=d[u]+1-T*x
            if d[v]<value:
                d[v]=value; parent[v]=e; changed=v
        if changed<0:
            return {'status':'FEASIBLE', 'potential':d}
    v=changed
    for _ in range(n):
        require(parent[v]>=0, 'PARENT_MISSING'); v=parent[v]//q
    root=v; backwards=[]
    while True:
        e=parent[v]; require(e>=0,'CYCLE_PARENT'); backwards.append(e);v=e//q
        if v==root: break
        require(len(backwards)<=n, 'CYCLE_EXTRACTION')
    cycle=list(reversed(backwards))
    require(sum(1-T*m[e] for e in cycle)>0, 'CYCLE_NOT_POSITIVE')
    return {'status':'CYCLE_OBSTRUCTION','cycle_edges':cycle}

def rounded_rc(m: list[int], q: int, k: int, T: int,
               d: list[int], epsilon: int=0) -> dict:
    """Round a feasible ordinary potential in either RC parity sector.
    The tight-SCC obstruction is the UTVPI mechanism of Lahiri--Musuvathi.
    The reduction to MDS and the offset T-1-T*epsilon are proved in the note.
    """
    n=q**(k-1); require(k>=2 and k%2==0 and q>=2 and q%2==0,'RC_DOMAIN')
    require(epsilon in (0,1) and len(d)==n and all(type(x) is int for x in d),'RC_INPUT')
    require(all(m[e]==m[theta(e,k,q)] for e in range(q*n)), 'INPUT_NOT_RC')
    require(all(d[e%n]-d[e//q]>=1-T*m[e] for e in range(q*n)), 'INFEASIBLE_INPUT_POTENTIAL')
    th=[theta(v,k-1,q) for v in range(n)]; c=T-1-T*epsilon
    D=[d[v]-d[th[v]]+c for v in range(n)]
    adj=[[] for _ in range(n)]; rev=[[] for _ in range(n)]
    for e,x in enumerate(m):
        u=e//q;v=e%n
        require(D[v]-D[u]>=2*(1-T*x),'HALF_FEASIBILITY')
        if D[u]%2 and D[v]%2 and D[v]-D[u]==2*(1-T*x):
            adj[u].append(e);rev[v].append(e)
    seen=bytearray(n);order=[]
    for root in range(n):
        if not D[root]%2 or seen[root]: continue
        seen[root]=1;stack=[(root,0)]
        while stack:
            u,j=stack[-1]
            if j<len(adj[u]):
                e=adj[u][j];stack[-1]=(u,j+1);v=e%n
                if not seen[v]:seen[v]=1;stack.append((v,0))
            else:
                order.append(u);stack.pop()
    comp=[-1]*n;next_id=0
    for root in reversed(order):
        if comp[root]>=0:continue
        comp[root]=next_id;stack=[root]
        while stack:
            v=stack.pop()
            for e in rev[v]:
                u=e//q
                if comp[u]<0:comp[u]=next_id;stack.append(u)
        next_id+=1
    bad=next((v for v in range(n) if D[v]%2 and comp[v]==comp[th[v]]),None)
    if bad is not None:
        def path(s,t):
            previous={s:-1};queue=deque([s])
            while queue and t not in previous:
                u=queue.popleft()
                for e in adj[u]:
                    v=e%n
                    if v not in previous:previous[v]=e;queue.append(v)
            require(t in previous,'PARITY_PATH_MISSING');out=[];v=t
            while v!=s:e=previous[v];out.append(e);v=e//q
            return out[::-1]
        return {'status':'PARITY_OBSTRUCTION','T':T,'epsilon':epsilon,'vertex':bad,
                'path_to_mate':path(bad,th[bad]),'path_back':path(th[bad],bad)}
    H=[D[v]//2 if not D[v]%2 else
       (D[v]+(1 if comp[v]>comp[th[v]] else -1))//2 for v in range(n)]
    require(all(H[v]+H[th[v]]==c for v in range(n)), 'RC_INTEGER_SUM')
    require(all(H[e%n]-H[e//q]>=1-T*m[e] for e in range(q*n)), 'RC_INTEGER_DIFF')
    f=[-(x//T) for x in H]; target=[m[e]+f[e//q]-f[e%n] for e in range(q*n)]
    require(all(x in (0,1) for x in target),'RETIMING_NOT_BOOLEAN')
    require(all(f[v]+f[th[v]]==epsilon for v in range(n)),'RETIMING_SECTOR')
    require(all(target[e]==target[theta(e,k,q)] for e in range(q*n)), 'RETIMING_NOT_RC')
    rank=[x%T for x in H]
    require(all(target[e] or rank[e%n]>rank[e//q] for e in range(q*n)),'RETIMING_RANK')
    return {'status':'FEASIBLE_RC','T':T,'epsilon':epsilon,'H':H,
            'cut_potential':f,'indicator':target,'bounded_rank':rank}

def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bits',type=Path);p.add_argument('--q',type=int,required=True)
    p.add_argument('--k',type=int,required=True);p.add_argument('--T',type=int,required=True)
    p.add_argument('--epsilon',type=int,choices=(0,1),default=0)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    require(a.q>=2 and a.k>=2,'DOMAIN');m=read_bits(a.bits,a.q**a.k)
    pcr_check(m,a.q,a.k);dag(m,a.q,a.k)
    r=bellman_ford(m,a.q,a.k,a.T)
    if r['status']=='FEASIBLE':r=rounded_rc(m,a.q,a.k,a.T,r['potential'],a.epsilon)
    a.out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(r['status'])
if __name__=='__main__':
    main()
