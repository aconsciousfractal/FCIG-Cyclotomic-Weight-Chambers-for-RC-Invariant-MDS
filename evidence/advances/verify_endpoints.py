"""Independent base-four verification of the packaged endpoints.

Usage: python -B verify_endpoints.py --out RESULT.json
No producer/checker imports, solver, floating point or Python assertions.
Only the explicitly listed endpoints are tested as positive witnesses.
"""
from pathlib import Path
from collections import Counter, deque
import argparse, hashlib, json, os, re, zipfile

def require(ok, why):
    if not ok: raise ValueError(why)

def encode(word):
    value=0
    for a in word:
        require(a in 'ACGT','ALPHABET')
        value=4*value+'ACGT'.index(a)
    return value

def reverse_complement(e,k=8):
    r=0
    for _ in range(k):r=4*r+3-e%4;e//=4
    return r

def graph_data(q=4,k=8):
    n=q**(k-1);size=q*n
    pcr=[-1]*size;orbits=[]
    for e in range(size):
        if pcr[e]>=0:continue
        orbit=[];w=e
        while w not in orbit:orbit.append(w);w=(w%n)*q+w//n
        for w in orbit:pcr[w]=len(orbits)
        orbits.append(orbit)
    return n,pcr,orbits

def evaluate(selected,n,pcr,orbits):
    present=bytearray(4*n)
    for e in selected:present[e]=1
    deg=[0]*n;dist=[0]*n
    for e in range(4*n):
        if not present[e]:deg[e%n]+=1
    queue=deque(v for v in range(n) if not deg[v]);count=0
    while queue:
        u=queue.popleft();count+=1
        for a in range(4):
            e=4*u+a
            if present[e]:continue
            v=e%n;dist[v]=max(dist[v],dist[u]+1);deg[v]-=1
            if not deg[v]:queue.append(v)
    counts=Counter(pcr[e] for e in selected)
    return {'cardinality':len(selected),'rc':all(reverse_complement(e) in selected for e in selected),
            'pcr_hit':len(counts),'pcr_histogram':dict(sorted(Counter(counts.values()).items())),
            'acyclic':count==n,'L':max(dist) if count==n else None}

def walk(words,selected,cycle=False):
    require(bool(words),'EMPTY_WALK');es=[encode(w) for w in words]
    require(all(len(w)==8 for w in words),'WORD_LENGTH')
    require(all(es[i]%16384==es[i+1]//4 for i in range(len(es)-1)),'ADJACENCY')
    if cycle:
        require(es[-1]%16384==es[0]//4,'CYCLE_CLOSURE')
        require(len({e//4 for e in es})==len(es),'SIMPLE_CYCLE')
    return sum(e in selected for e in es)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    root=Path(__file__).resolve().parent;n,pcr,orbits=graph_data();require(len(orbits)==8230,'PCR_COUNT')
    targets=[(1,'witnesses/M121.words.txt',8230,121),(2,'witnesses/M114.words.txt',8230,114),(3,'witnesses/M111.words.txt',8230,111),(4,'witnesses/M109.words.txt',8230,109),(5,'witnesses/M107.words.txt',8230,107),(10,'certificates/DS8234_L105.words.txt',8234,105),(11,'certificates/S8232.words.txt',8232,105),(13,'certificates/S8232_reuse_02.words.txt',8232,105),(13,'certificates/S8232_reuse_03.words.txt',8232,105)]
    def words(r,name):
        ls=(root/f'r{r}'/name).read_text().splitlines();require(len(ls)==len(set(ls)) and all(len(w)==8 for w in ls),'WORDS');return {encode(w) for w in ls}
    checked=[]
    for r,name,size,L in targets:
        selected=words(r,name);data=evaluate(selected,n,pcr,orbits)
        require(data['rc'] and data['acyclic'] and data['cardinality']==size and data['pcr_hit']==8230 and data['L']==L,'ENDPOINT')
        if r<=5:
            c=json.loads((root/f'r{r}'/'certificates'/Path(name).name.replace('.words.txt','.json')).read_text());cy=c['critical_cycle']
            require(walk(cy,selected,True)==1 and len(cy)==L+1,'COMPONENT_OPTIMUM');data['component_optimum']=L
        checked.append({'path':f'r{r}/'+name,**data})
    mstar=words(1,'witnesses/Mstar125.words.txt');m124=words(1,'witnesses/M124_sector1.words.txt')
    for selected,L in [(mstar,125),(m124,124)]:
        d=evaluate(selected,n,pcr,orbits);require(d['rc'] and d['acyclic'] and d['cardinality']==8230 and d['L']==L,'PARITY_ENDPOINT')
    c=json.loads((root/'r1/certificates/parity_obstruction.json').read_text());p=c['path_to_mate'];q=c['path_back']
    require(len(p)==71 and len(q)==54 and walk(p,mstar)==0 and walk(q,mstar)==1 and walk(p+q,mstar,True)==1,'PARITY_PATHS')
    require(p[0][:-1]=='AAAAAAA' and p[-1][1:]=='TTTTTTT','PARITY_ENDPOINTS')
    f=json.loads((root/'r1/certificates/retiming_Mstar125_to_M124.json').read_text())['cut_potential']
    require(len(f)==n and all(type(x) is int for x in f),'POTENTIAL')
    require(all(int(e in m124)-int(e in mstar)==f[e//4]-f[e%n] for e in range(4*n)),'GRADIENT')
    sums={f[v]+f[reverse_complement(v,7)] for v in range(n)};require(len(sums)==1 and next(iter(sums))%2==1,'SECTOR')
    result={'status':'PASS','targets':checked,'parity':{'ordinary_optimum':124,'sector0_optimum':125,'sector1_optimum':124,'contradiction':'2H(AAAAAAA)=53'}}
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result))
if __name__=='__main__':main()
