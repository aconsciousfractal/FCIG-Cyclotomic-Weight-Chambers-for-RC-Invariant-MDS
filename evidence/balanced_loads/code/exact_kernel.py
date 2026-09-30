"""Certificate kernel: standard library; no discovery or inherited graph import.
Base-q indexing independently reconstructs the graph, PCRs and row coefficients.
"""
from pathlib import Path
from collections import Counter, deque
from fractions import Fraction
from math import lcm
import gzip, hashlib, json, struct
ROOT=Path(__file__).resolve().parents[1]
def need(ok, why):
    if not ok: raise ValueError(why)
def load(path):
    path=Path(path)
    if path.suffix=='.gz':
        with gzip.open(path,'rt',encoding='utf8') as f:return json.load(f)
    return json.loads(path.read_text(encoding='utf8'))
def save(path,obj):
    path=Path(path)
    data=(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()
    if path.suffix=='.gz':
        with path.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as z:z.write(data)
    else:path.write_bytes(data)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class Net:
    def __init__(self,q=4,k=8):
        need(q in (2,4) and k>=2 and k%2==0,'carrier hypothesis')
        self.q,self.k=q,k;self.n=q**(k-1);self.m=q**k
        self.alphabet='ACGT' if q==4 else '01'
        self.rc=[self.rev(i,k) for i in range(self.m)]
        self.theta=[self.rev(i,k-1) for i in range(self.n)]
        seen=set();self.orbits=[];self.pof=[-1]*self.m
        for e in range(self.m):
            if e in seen:continue
            p=set();u=e
            for _ in range(k):p.add(u);u=q*(u%self.n)+u//self.n
            p=sorted(p);idx=len(self.orbits);self.orbits.append(p)
            for u in p:self.pof[u]=idx
            seen.update(p)
        self.stable_p=[p for p in self.orbits if self.rc[p[0]] in p]
        self.S=set(e for p in self.stable_p for e in p)
        self.E0=[e for e in range(self.m) if e not in self.S]
        self.stable_of={e:i for i,p in enumerate(self.stable_p) for e in p}
        self.groups=[p for p in self.orbits if p[0] not in self.S and p[0]<min(self.rc[e] for e in p)]
        self.options=[e for p in self.groups for e in p]
        self.opt={e:i for i,u in enumerate(self.options) for e in {u,self.rc[u]}}
        self.group_of={e:i for i,p in enumerate(self.groups) for u in p for e in {u,self.rc[u]}}
        self.greps=sorted(e for e in self.S if e<=self.rc[e]);self.gi={e:i for i,u in enumerate(self.greps) for e in {u,self.rc[u]}}
        self.tree=self.make_tree()
    def rev(self,x,k):
        r=0
        for _ in range(k):r=r*self.q+self.q-1-x%self.q;x//=self.q
        return r
    def word(self,x,k=None):
        k=self.k if k is None else k
        return ''.join(self.alphabet[(x//self.q**i)%self.q] for i in range(k-1,-1,-1))
    def encode(self,w):
        x=0
        for a in w:x=self.q*x+self.alphabet.index(a)
        return x
    def make_tree(self):
        parent=[None]*self.n;order=[]
        for root in range(self.n):
            if parent[root] is not None:continue
            parent[root]=(-1,-1,0);todo=deque([root])
            while todo:
                u=todo.popleft();order.append(u)
                edges=[(self.q*u+a,1) for a in range(self.q)]+[(a*self.n+u,-1) for a in range(self.q)]
                for e,sign in edges:
                    v=e%self.n if sign==1 else e//self.q
                    if e not in self.S and parent[v] is None:parent[v]=(u,e,sign);todo.append(v)
        return order,parent
    def gauge(self,x):
        order,parent=self.tree;pot=[0]*self.n
        for v in order:
            u,e,sgn=parent[v]
            if u>=0:pot[v]=pot[u]-sgn*x[e]
        return [x[e]+pot[e%self.n]-pot[e//self.q] for e in self.E0]
    def class_key(self,x):
        y=self.gauge(x)
        return hashlib.sha256(struct.pack('<'+'q'*len(y),*y)).hexdigest()
    def indicator(self,selection):
        x=[0]*self.m
        for e in selection:need(type(e) is int and 0<=e<self.m,'arc ID');x[e]=1
        return x
    def wordsfile(self,path):
        ws=Path(path).read_text(encoding='utf8').split();need(len(ws)==len(set(ws)),'duplicate word')
        need(all(len(w)==self.k for w in ws),'word length')
        return {self.encode(w) for w in ws}
    def transfer(self,u,mask):
        need(type(u) is int and 0<=u<self.n and 0<mask<2**self.q,'transfer parameters')
        d=Counter()
        for a in range(self.q):
            if mask>>a&1:
                incoming=a*self.n+u;outgoing=self.q*u+a
                need(incoming not in self.S and outgoing not in self.S,'transfer stable support')
                for e,c in ((incoming,-1),(outgoing,1),(self.rc[incoming],-1),(self.rc[outgoing],1)):d[e]+=c
        return {e:c for e,c in d.items() if c}
    def moved(self,selection,move):
        x=self.indicator(selection)
        for e,c in self.transfer(*move).items():x[e]+=c
        return x
    def cochain(self,x):
        need(len(x)==self.m and all(type(a) is int for a in x),'integer cochain')
        need(all(x[e]==x[self.rc[e]] for e in self.E0),'RC cochain')
        need(all(sum(x[e] for e in p)==1 for p in self.orbits if p[0] not in self.S),'PCR cochain sums')
    def walk(self,edges,closed=False):
        need(bool(edges) and all(type(e) is int and 0<=e<self.m for e in edges),'walk IDs')
        need(all(a%self.n==b//self.q for a,b in zip(edges,edges[1:])),'walk adjacency')
        if closed:need(edges[-1]%self.n==edges[0]//self.q,'walk closure')
    def core(self,sel,rank=None,path=None,L=None):
        need(isinstance(sel,list) and len(sel)==len(set(sel)),'unique selection')
        selected=set(sel);need(not(selected&self.S),'nonstable selection')
        x=self.indicator(sel);self.cochain(x)
        need(selected=={self.rc[e] for e in selected},'RC selection')
        inc=[0]*self.n;out=[[] for _ in range(self.n)]
        for e in self.E0:
            if e not in selected:out[e//self.q].append(e%self.n);inc[e%self.n]+=1
        queue=deque(v for v in range(self.n) if inc[v]==0);dist=[0]*self.n;seen=0
        while queue:
            u=queue.popleft();seen+=1
            for v in out[u]:
                dist[v]=max(dist[v],dist[u]+1);inc[v]-=1
                if inc[v]==0:queue.append(v)
        need(seen==self.n,'cyclic core')
        length=max(dist)
        if rank is not None:
            need(len(rank)==self.n and all(type(t) is int and 0<=t<=length for t in rank),'rank bounds')
            need(all(rank[v]>=rank[u]+1 for u in range(self.n) for v in out[u]),'rank inequality')
            need(rank==dist,'rank is exact longest incoming rank')
        if path is not None:
            self.walk(path);need(len(path)==length and all(e not in selected and e not in self.S for e in path),'attaining path')
        if L is not None:need(length==L,'core length')
        return dist,length
    def clock(self,x,H,T,C):
        self.cochain(x);need(type(T) is int and T>=1,'positive T')
        need(type(C) is int and C%T==T-1,'RC clock congruence')
        need(len(H)==self.n and all(type(v) is int for v in H),'integral clock')
        need(all(H[v]+H[self.theta[v]]==C for v in range(self.n)),'RC clock sum')
        sel=[]
        for e in self.E0:
            a=H[e//self.q];b=H[e%self.n]
            need(b-a+T*x[e]>=1,'clock inequality')
            t=x[e]+b//T-a//T;need(t in (0,1),'clock Boolean reconstruction')
            if t:sel.append(e)
        return sel

def cycle_key(w):return tuple(sorted(Counter(w).items()))
def pool(n,root=ROOT):
    raw=load(root/'inputs/return_sources.json.gz')['families'];lib=load(root/'inputs/returns.json.gz')['families']
    out=[];seen=set()
    def add(w):
        n.walk(w,True);need(any(e in n.S for e in w),'mixed walk')
        key=cycle_key(w)
        if key not in seen:seen.add(key);out.append(w)
    def perm(e,p):
        z=0
        for i in range(n.k-1,-1,-1):z=n.q*z+p[(e//n.q**i)%n.q]
        return z
    for f in lib:
        for row in raw[f['source_index']]['returns']:
            w=[perm(e,f['alphabet_permutation']) for e in [row['stable_edge']]+row['path']]
            add(w);add([n.rc[e] for e in reversed(w)])
    for row in load(root/'inputs/blocks.json.gz')['joined']['walks']:add(row['walk'])
    for w in list(out):add([n.rc[e] for e in reversed(w)])
    return out

def compile_rows(n,walks,T=107):
    rows=[];seen={};mapping=[]
    def add(ext,sta,rhs,sense,kind):
        ext=sorted([i,c] for i,c in ext.items() if c);sta=sorted([i,c] for i,c in sta.items() if c)
        key=json.dumps([ext,sta,rhs,sense],separators=(',',':'))
        if key not in seen:seen[key]=len(rows);rows.append(dict(external_terms=ext,stable_terms=sta,rhs=rhs,sense=sense,kind=kind))
        return seen[key]
    for w in walks:
        n.walk(w,True);need(any(e in n.S for e in w),'mixed primitive')
        mapping.append(add(Counter(n.opt[e] for e in w if e not in n.S),Counter(n.gi[e] for e in w if e in n.S),(len(w)+T-1)//T,'>=','rounded_walk'))
    for p in n.stable_p:add(Counter(),Counter(n.gi[e] for e in p),1,'>=','PCR_lower')
    add(Counter(),Counter(n.gi[e] for e in n.S),len(n.stable_p)+1,'<=','stable_budget')
    return rows,mapping

def sparse_ratios(items,N):
    need(len({a[0] for a in items})==len(items),'duplicate rational coordinate')
    den=1
    for i,num,d in items:
        need(type(i) is int and 0<=i<N and type(num) is int and type(d) is int and d>0,'rational encoding')
        den=lcm(den,d)
    out=[0]*N
    for i,num,d in items:out[i]=num*(den//d)
    return out,den

def lp_certificate(n,rows,x,cert):
    nr=len(rows)-1;G=len(n.greps)
    gi,gd=sparse_ratios(cert['gamma'],G);la,ld=sparse_ratios(cert['lambda'],nr)
    need(all(v>=0 for v in la),'dual nonnegativity')
    base=[r['rhs']-sum(c*x[n.options[i]] for i,c in r['external_terms']) for r in rows[:-1]]
    need(all(sum(c*gi[i] for i,c in row['stable_terms'])>=gd*b for row,b in zip(rows[:-1],base)),'LP primal infeasible')
    comb=[0]*G
    for lam,row in zip(la,rows[:-1]):
        if lam:
            for i,c in row['stable_terms']:comb[i]+=lam*c
    cost=[0]*G
    for i,c in rows[-1]['stable_terms']:cost[i]=c
    need(comb==[ld*c for c in cost],'LP dual coefficient identity')
    pv=Fraction(sum(c*v for c,v in zip(cost,gi)),gd);dv=Fraction(sum(v*b for v,b in zip(la,base)),ld)
    need(pv==dv==Fraction(*cert['optimum']),'LP exact optimum')
    return pv

def lift_check(n,rows,selection,gamma):
    need(len(gamma)==len(n.greps),'lift dimension')
    xs=n.indicator(selection);gg,den=sparse_ratios([[i,*v] for i,v in enumerate(gamma)],len(gamma))
    for row in rows:
        lhs=den*sum(c*xs[n.options[i]] for i,c in row['external_terms'])+sum(c*gg[i] for i,c in row['stable_terms'])
        need(lhs>=den*row['rhs'] if row['sense']=='>=' else lhs<=den*row['rhs'],'lift row')
    return Fraction(sum(gg[n.gi[e]] for e in n.S),den)

def count_tight(n,selected,r,L):
    """Count all simple E0 cycles of L+1 arcs with exactly one selected arc."""
    selected=set(selected);need(max(r)==L,'tight count rank range')
    heads={}
    for e in selected:
        if r[e//n.q]==L and r[e%n.n]==0:heads.setdefault(e%n.n,[]).append(e//n.q)
    order=sorted(range(n.n),key=lambda v:r[v]);answer=0
    for s,targets in heads.items():
        a=[0]*n.n;a[s]=1
        for u in order:
            if not a[u]:continue
            for b in range(n.q):
                e=n.q*u+b;v=e%n.n
                if e not in n.S and e not in selected and r[v]==r[u]+1:a[v]+=a[u]
        answer+=sum(a[v] for v in targets)
    return answer
