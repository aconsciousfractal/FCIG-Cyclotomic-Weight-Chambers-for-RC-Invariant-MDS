"""Exact arithmetic and lower de Bruijn data. Standard library only.
No floating point; quadratic numbers (a+b*sqrt(2))/d, d>0.
The fixed trig tables are elementary coordinates, not empirical data.
"""
from __future__ import annotations
from pathlib import Path
from itertools import product
from collections import deque
import json, math, hashlib, zipfile, tempfile
class CheckError(RuntimeError): pass
def check(ok,code):
    if not ok: raise CheckError(code)
def dump(path:Path,obj):path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
def sha(b):return hashlib.sha256(b).hexdigest()
def sign(a:int,b:int)->int:
    if a==0:return (b>0)-(b<0)
    if b==0:return (a>0)-(a<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    v=a*a-2*b*b
    return ((v>0)-(v<0))*((a>0)-(a<0))
def value(form,r):
    a,b,c,d=form;x,y,z=r
    check(z>0,'DENOMINATOR_NOT_POSITIVE')
    return a*x+2*c*y+b*z,a*y+c*x+d*z
def ev(form,r):return sign(*value(form,r))
def der(form):return sign(form[0],form[2])
def rc(u,n,q=4):
    r=0
    for _ in range(n):r=q*r+q-1-u%q;u//=q
    return r
def word(u,n,q=4):
    w=[0]*n
    for i in range(n-1,-1,-1):w[i]=u%q;u//=q
    return w
def text(u,n,q=4):return ''.join(('ACGT' if q==4 else '01')[a] for a in word(u,n,q))
def outedges(u,n,q=4):return range(q*u,q*u+q)
def inedges(u,n,q=4):return [a*n+u for a in range(q)]
def unpack(bits,n):
    check(len(bits)==(n+7)//8,'BITS_SIZE')
    return bytearray((bits[i//8]>>(i%8))&1 for i in range(n))
def pack(m):
    b=bytearray((len(m)+7)//8)
    for i,x in enumerate(m):
        check(x in (0,1),'PACK_NOT_BOOLEAN')
        b[i//8]|=x<<(i%8)
    return bytes(b)
V=((-1,0),(0,-1),(0,1),(1,0))
TRIG={
4: ([(0,0),(-2,0),(0,0)],[(1,0),(0,0),(-1,0)]),
6: ([(1,0),(-1,0),(-2,0),(-1,0),(1,0)],[(1,0),(1,0),(0,0),(-1,0),(-1,0)]),
8: ([(0,1),(0,0),(0,-1),(-2,0),(0,-1),(0,0),(0,1)],[(1,0),(0,1),(1,0),(0,0),(-1,0),(0,-1),(-1,0)])}
def form(w,coeff):
    a=b=c=d=0
    for ch,(x,y) in zip(w,coeff):
        h,l=V[ch];a+=h*x;b+=l*x;c+=h*y;d+=l*y
    return a,b,c,d
def forms(k):
    check(k in TRIG,'UNSUPPORTED_ORDER')
    C,S=TRIG[k];ww=list(product(range(4),repeat=k-1))
    return [form(w,S) for w in ww],[form(w,[(x-2,y) for x,y in C]) for w in ww],[form(w,C) for w in ww]
def state(k,r,ff):
    lf,rf,_=ff;n=4**(k-1)
    l=[ev(f,r) for f in lf];s=[ev(f,r) for f in rf]
    a=[int(l[u]>0 and (s[u]<0 or s[u]==0 and u<rc(u,k-1))) for u in range(n)]
    return l,s,a
def singular_cycles(k,Z):
    n=4**(k-1);ad={};inc={}
    for u in sorted(Z):
        es=[e for e in outedges(u,n) if e%n in Z]
        check(len(es)<=1,'SINGULAR_BRANCH')
        ies=[e for e in inedges(u,n) if e//4 in Z]
        check(len(ies)==len(es),'SINGULAR_UNBALANCED')
        if es:ad[u]=es[0];inc[u]=ies[0]
    seen=set();cycles=[]
    for u in sorted(ad):
        if u in seen:continue
        c=[];x=u
        while x not in seen:
            seen.add(x);e=ad[x];c.append(e);x=e%n
        check(x==u,'SINGULAR_TAIL')
        check(len(c)<=k and k%len(c)==0,'SINGULAR_PERIOD')
        check(all(c[(i+1)%len(c)]==(e%n)*4+e//n for i,e in enumerate(c)),'SINGULAR_NOT_PCR')
        cycles.append(c)
    return cycles

def dag(m,n,q=4,allowed=None):
    if allowed is None:allowed=[True]*n
    ind=[0]*n
    for e,x in enumerate(m):
        if not x and allowed[e//q] and allowed[e%n]:ind[e%n]+=1
    todo=deque(u for u in range(n) if allowed[u] and ind[u]==0);order=[];dist=[0]*n;pred=[-1]*n
    while todo:
        u=todo.popleft();order.append(u)
        for e in outedges(u,n,q):
            v=e%n
            if m[e] or not allowed[v]:continue
            if dist[v]<dist[u]+1:dist[v]=dist[u]+1;pred[v]=u
            ind[v]-=1
            if ind[v]==0:todo.append(v)
    check(len(order)==sum(allowed),'DAG_CYCLE')
    if not order:return [],0,[],dist
    t=max(order,key=lambda u:dist[u]);path=[t]
    while pred[path[-1]]>=0:path.append(pred[path[-1]])
    path.reverse();return order,dist[t],path,dist

def parent_data(parent:Path,k):
    a=json.loads((parent/f'results/atlas_k{k}.json').read_text());n=4**k
    ts=[unpack((parent/f"results/certificates/{r['id']}.bits").read_bytes(),n) for r in a['rows']]
    return a,ts

def locate_parent(path:Path|None=None):
    parent=path if path is not None else Path(__file__).resolve().parents[2]/'atlas'
    check((parent/'results/atlas_k8.json').is_file(),'MISSING_ATLAS')
    return parent,None
