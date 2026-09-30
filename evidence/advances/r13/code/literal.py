"""Independent literal-word finite checker primitives. Standard library only.
No imports from discovery, NumPy, SciPy, Numba, or any optimization solver.
The lower graph retains loops and all isolated vertices after deletion.
"""
from pathlib import Path
from itertools import product
from collections import Counter, deque
from functools import lru_cache
import gzip, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
ABC='ACGT';TR=str.maketrans(ABC,'TGCA')
V=tuple(''.join(t) for t in product(ABC,repeat=7));W=tuple(v+a for v in V for a in ABC)
WS=frozenset(W)
def rc(w):return w.translate(TR)[::-1]
def need(c,m):
    if not c:raise ValueError(m)
def load(p):
    p=Path(p);b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def save(p,obj):Path(p).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def digest(s):return hashlib.sha256(('\n'.join(sorted(s))+'\n').encode('ascii')).hexdigest()
def safe_path(p):
    q=(ROOT/p).resolve();need(q.is_relative_to(ROOT),'PATH_ESCAPE');return q

def read_set(p):
    data=Path(p).read_text(encoding='ascii').splitlines()
    need(data and all(w in WS for w in data),'SET_WORDS')
    need(len(data)==len(set(data)),'SET_DUPLICATES')
    s=frozenset(data);need({rc(w) for w in s}==s,'SET_RC');return s

@lru_cache(None)
def pcrs():
    seen=set();out=[]
    for w in W:
        if w in seen:continue
        p=tuple(sorted({w[i:]+w[:i] for i in range(8)}));out.append(p);seen.update(p)
    need(len(out)==8230 and seen==WS,'PCR_CENSUS');return tuple(out)
@lru_cache(None)
def groups():
    out=[];seen=set()
    for p in pcrs():
        if p[0] in seen:continue
        q=tuple(sorted(rc(w) for w in p));full=frozenset(p)|frozenset(q);seen.update(full)
        out.append({'representative':min(full),'words':full,'pcr':p,'stable':p==q})
    need(len(out)==4179 and sum(g['stable'] for g in out)==128,'GROUP_CENSUS');return out

def walk(p,closed=False,simple=True):
    need(type(p) is list and p and all(w in WS for w in p),'WALK_WORDS')
    need(all(x[1:]==y[:-1] for x,y in zip(p,p[1:])),'WALK_OVERLAP')
    if closed:need(p[-1][1:]==p[0][:-1],'CYCLE_CLOSE')
    vertices=[w[:-1] for w in p]+([] if closed else [p[-1][1:]])
    if simple:need(len(vertices)==len(set(vertices)),'WALK_SIMPLE')

def graph(s):
    adj={v:[] for v in V};indeg=dict.fromkeys(V,0)
    for w in W:
        if w not in s:adj[w[:-1]].append(w);indeg[w[1:]]+=1
    q=deque(v for v in V if indeg[v]==0);order=[];depth=dict.fromkeys(V,0);parent={}
    while q:
        u=q.popleft();order.append(u)
        for w in adj[u]:
            v=w[1:]
            if depth[v]<depth[u]+1:depth[v]=depth[u]+1;parent[v]=w
            indeg[v]-=1
            if indeg[v]==0:q.append(v)
    need(len(order)==len(V),'RESIDUAL_DIRECTED_CYCLE')
    suffix=dict.fromkeys(V,0)
    for u in reversed(order):suffix[u]=max([0]+[1+suffix[w[1:]] for w in adj[u]])
    end=max(V,key=depth.get);path=[]
    while end in parent:
        w=parent[end];path.append(w);end=w[:-1]
    return depth,suffix,path[::-1]

def check_endpoint(c):
    need(c.get('schema')=='R13-ENDPOINT-1','ENDPOINT_SCHEMA')
    s=read_set(safe_path(c['selected_file']))
    need(c.get('selected_sha256')==digest(s),'ENDPOINT_HASH')
    need(type(c.get('cardinality')) is int and c['cardinality']==len(s),'ENDPOINT_CARDINALITY')
    L=c.get('L');need(type(L) is int and L>=0,'L_TYPE')
    rank=c.get('rank');need(type(rank) is list and len(rank)==len(V),'RANK_SIZE')
    need(all(type(x) is int and 0<=x<=L for x in rank),'RANK_DOMAIN')
    d=dict(zip(V,rank))
    for w in W:
        if w not in s:need(d[w[1:]]>=d[w[:-1]]+1,'RANK_EDGE')
    path=c.get('path');walk(path)
    need(len(path)==L and not(set(path)&s),'PATH_LENGTH_OR_SELECTED')
    text=path[0][:-1]+''.join(w[-1] for w in path)
    need(c.get('avoiding_word')==text and len(text)==L+7,'AVOIDING_WORD')
    need(all(text[i:i+8] not in s for i in range(L)),'AVOIDING_WINDOWS')
    hist=Counter(str(len(set(p)&s)) for p in pcrs())
    need('0' not in hist and dict(hist)==c.get('PCR_histogram'),'PCR_HISTOGRAM')
    rd,suffix,rpath=graph(s)
    need(rank==[rd[v] for v in V] and len(rpath)==L,'RANK_RECOMPUTE')
    return s,rd,suffix
