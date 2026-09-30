#!/usr/bin/env python3
"""Independent exact checker. Python standard library only.
Uses literal DNA words, not the numerical discovery graph routines.
The certificate is data. No solver status is accepted as a lower-bound proof.
"""
from __future__ import annotations
from pathlib import Path
from collections import deque,Counter
from itertools import product,combinations
import json,hashlib,time,argparse,os
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('R9_CHECK_OUTPUT',ROOT/'run_output'));OUT.mkdir(parents=True,exist_ok=True)
ABC='ACGT';TRAN=str.maketrans('ACGT','TGCA')
V=[''.join(x) for x in product(ABC,repeat=7)]
EDGES=[v+a for v in V for a in ABC]

def req(ok,msg):
    if not ok:raise ValueError(msg)

def rc(w):return w.translate(TRAN)[::-1]
def decode(i):
    req(type(i)is int and 0<=i<65536,'EDGE_INTEGER')
    s=''
    for _ in range(8):s=ABC[i%4]+s;i//=4
    return s

def read_words(p):
    raw=p.read_text().splitlines();req(len(raw)==len(set(raw))==8230,'CARDINALITY')
    req(all(len(w)==8 and set(w)<=set(ABC) for w in raw),'WORD_ALPHABET')
    m=set(raw);req({rc(w) for w in m}==m,'RC_MEMBERSHIP');return m

def necklaces(m):
    seen=set();stable=[];count=0
    for w in EDGES:
        if w in seen:continue
        p={w[i:]+w[:i] for i in range(8)};seen|=p;count+=1
        req(len(p&m)==1,'ONE_PER_PCR')
        if {rc(x) for x in p}==p:stable.append(p)
    req(count==8230 and len(stable)==128,'PCR_COUNTS')
    req(all(sum(rc(w)==w for w in p)==2 for p in stable),'TWO_FIXED_OPTIONS')
    return stable

def graph_data(kept):
    adj={v:[] for v in V};ind={v:0 for v in V};d={v:0 for v in V};par={}
    for w in sorted(kept):adj[w[:-1]].append(w);ind[w[1:]]+=1
    queue=deque(v for v in V if not ind[v]);order=[]
    while queue:
        u=queue.popleft();order.append(u)
        for w in adj[u]:
            v=w[1:]
            if d[v]<d[u]+1:d[v]=d[u]+1;par[v]=w
            ind[v]-=1
            if not ind[v]:queue.append(v)
    req(len(order)==len(V),'RESIDUAL_DIRECTED_CYCLE')
    suffix={v:0 for v in V}
    for u in reversed(order):suffix[u]=max([0]+[1+suffix[w[1:]] for w in adj[u]])
    end=max(V,key=d.__getitem__);path=[];u=end
    while u in par:w=par[u];path.append(w);u=w[:-1]
    path.reverse();return adj,order,d,suffix,path

def walk_check(path,closed=False):
    req(bool(path) and all(w in EDGESET for w in path),'WALK_LABELS')
    req(all(a[1:]==b[:-1] for a,b in zip(path,path[1:])),'WALK_CONTINUITY')
    if closed:req(path[-1][1:]==path[0][:-1],'WALK_CLOSURE')

EDGESET=set(EDGES)
def verify_endpoint(p):
    m=read_words(p);necklaces(m);adj,order,d,b,path=graph_data(EDGESET-m)
    req(all(d[w[1:]]>=d[w[:-1]]+1 for w in EDGESET-m),'ENDPOINT_RANK')
    walk_check(path)
    return {'file':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'L':len(path),'cardinality':len(m)},m

def verify_cube(c,base,stable,direct=True):
    req(c['base_sha256']==hashlib.sha256((ROOT/'inputs/M107.words.txt').read_bytes()).hexdigest(),'CUBE_BASE_HASH')
    ew=[decode(i) for i in c['edges']];req(len(ew)==len(set(ew))==256,'CUBE_FIXED_EDGES');req(ew==c['edge_words'],'FIXED_WORD_LABELS')
    req(set(ew)=={w for w in EDGES if rc(w)==w},'ALL_FIXED_EDGES')
    PCRmap={w:frozenset(p) for p in stable for w in p};seen=set()
    for i in range(0,256,2):
        p=PCRmap[ew[i]];req(p==PCRmap[ew[i+1]] and p not in seen,'CUBE_VARIABLES');seen.add(p)
    D=EDGESET-base-set(ew);req({rc(w) for w in D}==D,'FIXED_D_RC')
    adj,order,a,b,path=graph_data(D);req(len(path)==c['fixed_L'],'FIXED_L')
    req(c.get('fixed_longest_path')==path,'FIXED_PATH_CERTIFICATE')
    tailindices={v:[] for v in V}
    for i,w in enumerate(ew):tailindices[w[:-1]].append(i)
    reach={v:sum(1<<i for i in tailindices[v]) for v in V}
    for v in reversed(order):
        for w in adj[v]:reach[v]|=reach[w[1:]]
    costs=[a[w[:-1]]+1+b[w[1:]] for w in ew]
    req(costs==c['costs'] and all(cost%2 for cost in costs),'FIXED_ARC_COSTS')
    unary=[];pairs=[]
    for i,w in enumerate(ew):
        if reach[w[1:]]>>i&1:unary.append(i)
        for j in range(i):
            r=bool(reach[w[1:]]>>j&1);s=bool(reach[ew[j][1:]]>>i&1)
            req(r==s,'REFLECTED_REACHABILITY')
            if r:pairs.append([j,i])
    req(unary==c['acyclic_unary_forbidden'],'UNARY_CLAUSES')
    req(pairs==c['acyclic_pairs_forbidden'],'PAIR_CLAUSES')
    # In this concrete cube all forced values already follow from unary cycles.
    forced={i^1 for i in unary};free=sorted(set(range(128))-{i//2 for i in forced})
    req(len(forced)==121 and free==c['free_variables'],'CUBE_FREE_DIMENSION')
    req(sorted(forced)==c['acyclic_forced_literals'],'FORCED_LITERALS')
    census=[];hist=Counter();sets=[];fullrechecks=0;component107=0
    cy=json.loads((ROOT/'inputs/M107.json').read_text())['critical_cycle'];walk_check(cy,True)
    req(len(cy)==108 and len({w[:-1] for w in cy})==108,'BASE_CRITICAL_CYCLE')
    for bits in product([0,1],repeat=len(free)):
        selected=forced|{2*i+bit for i,bit in zip(free,bits)}
        if any(i in selected for i in unary) or any(i in selected and j in selected for i,j in pairs):continue
        m=(base-set(ew))|{w for i,w in enumerate(ew) if i not in selected}
        req(len(m)==8230 and {rc(w) for w in m}==m,'CUBE_ENDPOINT')
        L=max([len(path)]+[costs[i] for i in selected]);req(L>=107,'CUBE_LOWER')
        # A full DAG traversal checks the formula at every one of the 128 survivors.
        if direct:
            _,_,dd,_,pp=graph_data(EDGESET-m);req(len(pp)==L,'CUBE_DIRECT_LONGEST_PATH');fullrechecks+=1
        idx=len(census); claimed=c['complete_census'][idx]; cyc=claimed.get('component_cycle',[]);walk_check(cyc,True)
        req(len({w[:-1] for w in cyc})==len(cyc),'COMPONENT_CYCLE_SIMPLE')
        hit=sum(w in m for w in cyc);req(hit==claimed.get('component_hits')==1 and len(cyc)==L+1,'COMPONENT_CYCLE_LOWER')
        census.append({'bits':list(bits),'L':L,'component_cycle':cyc,'component_hits':hit});hist[L]+=1;sets.append((bits,m))
        if L==107:
            req(sum(w in m for w in cy)==1,'CUBE_COMPONENT_107');component107+=1
    req(census==c['complete_census'],'COMPLETE_CENSUS');req({str(k):v for k,v in sorted(hist.items())}==c['length_profile'],'LENGTH_PROFILE')
    req(c['optimum']==107 and len(census)==128,'CUBE_OPTIMUM')
    req(c['lower_obstruction']=={'reason':'fixed_path','length':107},'CUBE_LOWER_OBJECT')
    chosen=set(c['restored_literals']);req(len(chosen)==128 and all((i in chosen)!=(i^1 in chosen) for i in range(0,256,2)),'OPTIMUM_CHOICES')
    witness=(base-set(ew))|{w for i,w in enumerate(ew) if i not in chosen};req(witness==read_words(ROOT/'certificates/stable_cube_optimum.words.txt'),'OPTIMUM_WITNESS_JOIN')
    # All vectors differ on <=7 bits, and agree on 121, hence are not complements.
    req(len({tuple(bits) for bits,m in sets})==128,'COMPONENT_LABELS_DISTINCT')
    return {'variables':128,'formal_assignments':str(2**128),'unary_forbidden':len(unary),'pair_forbidden':len(pairs),'forced_variables':121,'free_variables':free,'admissible_MDS':128,'length_profile':dict(sorted(hist.items())),'component_minimum_profile':dict(sorted(hist.items())),'component_minimum_107_witnesses':component107,'distinct_ordinary_components':128,'direct_DAG_rechecks':fullrechecks,'fixed_residual_L':len(path)}

def verify_boundary(c,base):
    opts=c['options'];req(len(opts)==len({tuple(o) for o in opts}),'PROJECTION_OPTIONS_UNIQUE')
    F={w for o in opts for w in o};req(all(w in EDGESET for w in F),'PROJECTION_WORDS')
    req({rc(w) for w in F}==F,'PROJECTION_RC')
    D=EDGESET-base-F;adj,order,a,b,path=graph_data(D)
    B=sorted({w[:-1] for w in F}|{w[1:] for w in F});req(B==c['boundary'],'BOUNDARY_SET')
    req(len(path)==c['fixed_L'],'BOUNDARY_FIXED_L')
    Bset=set(B);rows=[]
    for s in B:
        dis={s:0}
        for u in order:
            if u not in dis:continue
            if u!=s and u in Bset:
                if (s,u)<=(rc(u),rc(s)):rows.append([s,u,dis[u]])
                continue
            for w in adj[u]:dis[w[1:]]=max(dis.get(w[1:],-1),dis[u]+1)
    req(sorted(rows)==c['compressed_paths'],'BOUNDARY_PROJECTION_COMPLETE')
    req([[v,a[v],b[v]] for v in B]==c['boundary_bounds'],'BOUNDARY_BOUNDS')
    # Positive control with the exact base at L=107 and doubled half-integral ranks.
    _,_,d,_,_=graph_data(EDGESET-base);h2={v:d[v]+107-d[rc(v)] for v in V}
    req(all(2*a[v]<=h2[v]<=2*(107-b[v]) for v in B),'BOUNDARY_POSITIVE_BOUNDS')
    req(all(h2[v]-h2[u]>=2*ell for u,v,ell in rows),'BOUNDARY_POSITIVE_PATHS')
    req(all(h2[w[1:]]-h2[w[:-1]]>=2-216*(w in base) for w in F),'BOUNDARY_POSITIVE_ARCS')
    return {'boundary_vertices':len(B),'symmetric_rank_variables':len(B)//2,'options':len(opts),'free_arcs':len(F),'compressed_path_rows':len(rows),'positive_control_L':107}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args();start=time.time()
    base=read_words(ROOT/'inputs/M107.words.txt');stable=necklaces(base)
    result={'status':'PASS','schema':'round9-independent-verification-v1','endpoints':[]}
    for p in sorted((ROOT/'inputs').glob('*.words.txt'))+sorted((ROOT/'certificates').glob('*.words.txt')):
        row,_=verify_endpoint(p);result['endpoints'].append(row)
    c=json.loads((ROOT/'certificates/stable_cube.json').read_text());result['stable_cube']=verify_cube(c,base,stable)
    result['boundary_projections']={p.name:verify_boundary(json.loads(p.read_text()),base) for p in sorted((ROOT/'certificates').glob('boundary_*.json'))}
    oc=json.loads((ROOT/'certificates/oracle_107.json').read_text());om=read_words(ROOT/'certificates/oracle_107.words.txt');rank=oc['rank'];req(len(rank)==16384 and all(type(x)is int and 0<=x<=oc['L'] for x in rank),'ORACLE_RANK_DOMAIN');rd=dict(zip(V,rank))
    req(all(rd[w[1:]]>rd[w[:-1]] for w in EDGESET-om),'ORACLE_RANK_ARCS');walk_check(oc['longest_path']);req(len(oc['longest_path'])==oc['L']==107 and not(set(oc['longest_path'])&om),'ORACLE_LONGEST_PATH');req(set(oc['restored_words'])=={w for w in EDGES if rc(w)==w and w not in om},'ORACLE_RESTORED_WORDS')
    result['standalone_oracle_positive_control']={'target':oc['target'],'L':oc['L'],'status':'PASS'}
    result['seconds']=time.time()-start
    if args.out:args.out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
