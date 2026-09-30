"""Exact literal-word routines. No discovery or solver imports."""
from pathlib import Path
from itertools import product,combinations
from collections import deque,Counter
import json,hashlib,gzip
ROOT=Path(__file__).resolve().parents[1]
ABC='ACGT';TR=str.maketrans('ACGT','TGCA')
V=[''.join(w) for w in product(ABC,repeat=7)]
WORDS=[u+a for u in V for a in ABC];WSET=set(WORDS);VI={u:i for i,u in enumerate(V)}
def need(p,s):
 if not p:raise ValueError(s)
def rc(w):return w.translate(TR)[::-1]
def decode(i):
 need(type(i)is int and 0<=i<65536,'INTEGER_EDGE');return WORDS[i]
def load(p):
 if str(p).endswith('.gz'):
  with gzip.open(p,'rt') as f:return json.load(f)
 return json.loads(Path(p).read_text())
def read_set(p,card=None):
 a=Path(p).read_text().splitlines();need(len(set(a))==len(a),'DUPLICATE_WORD');need(all(w in WSET for w in a),'DNA_WORDS')
 if card is not None:need(len(a)==card,'CARDINALITY')
 m=set(a);need({rc(w) for w in m}==m,'RC_SET');return m
def pcrs():
 seen=set();out=[]
 for w in WORDS:
  if w in seen:continue
  p=sorted({w[i:]+w[:i] for i in range(8)});seen.update(p);out.append(p)
 need(len(out)==8230,'PCR_CENSUS');return out
def graph(kept):
 a={u:[] for u in V};deg={u:0 for u in V};d={u:0 for u in V};par={}
 for w in sorted(kept):a[w[:-1]].append(w);deg[w[1:]]+=1
 q=deque(u for u in V if deg[u]==0);order=[]
 while q:
  u=q.popleft();order.append(u)
  for w in a[u]:
   v=w[1:]
   if d[v]<d[u]+1:d[v]=d[u]+1;par[v]=w
   deg[v]-=1
   if deg[v]==0:q.append(v)
 need(len(order)==len(V),'DIRECTED_CYCLE')
 b={u:0 for u in V}
 for u in reversed(order):b[u]=max([0]+[b[w[1:]]+1 for w in a[u]])
 t=max(V,key=d.get);path=[]
 while t in par:w=par[t];path.append(w);t=w[:-1]
 return a,order,d,b,path[::-1]
def walk(p,closed=False):
 need(len(p)>0 and all(w in WSET for w in p),'WALK_LABELS');need(all(x[1:]==y[:-1] for x,y in zip(p,p[1:])),'WALK_CONTINUITY')
 if closed:need(p[-1][1:]==p[0][:-1],'WALK_CLOSURE')
def height(c):
 h=c['heights'];H=c['half_height'];need(type(H)is int and H>=0 and c['target_L']==2*H,'EVEN_BOUND');need(len(h)==len(V),'HEIGHT_SIZE')
 need(all(type(t)is int and -H<=t<=H for t in h),'HEIGHT_INTEGER_BOUND');r=dict(zip(V,h));need(all(r[rc(v)]==-r[v] for v in V),'HEIGHT_ANTISYMMETRY')
 m={w for w in WORDS if r[w[1:]]<=r[w[:-1]]};need(len(m)==c['cardinality'],'HEIGHT_CARDINALITY');need({rc(w) for w in m}==m,'HEIGHT_RC')
 _,_,d,b,p=graph(WSET-m);need(len(p)<=2*H,'HEIGHT_LENGTH');counts=Counter(sum(w in m for w in q) for q in pcrs());need(0 not in counts,'PCR_HIT')
 need(len(m)-8230==sum((j-1)*n for j,n in counts.items()),'SURPLUS_IDENTITY')
 return {'cardinality':len(m),'L':len(p),'surplus':len(m)-8230,'PCR_cut_histogram':dict(sorted(counts.items()))},m,p
