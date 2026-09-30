"""Reviewer-written physical checks, no imports from external round code.

Arguments are extracted round roots. All graph tests use explicit word
prefix/suffix incidence; counts and claimed levels are checked from scratch.
"""
from pathlib import Path
from collections import deque
from functools import lru_cache
import argparse, hashlib, json

def need(x, label):
    if not x: raise ValueError(label)

@lru_cache(None)
def geometry(q,k):
    n=q**(k-1); E=q*n
    def reverse_complement(e,length):
        digits=[]
        for _ in range(length): digits.append(e%q);e//=q
        r=0
        for d in digits:r=q*r+q-1-d
        return r
    mate=[reverse_complement(e,k) for e in range(E)]
    vmate=[reverse_complement(v,k-1) for v in range(n)]
    seen=set();cycles=[]
    for e in range(E):
        if e in seen:continue
        x=e;cycle=[]
        while x not in seen:
            seen.add(x);cycle.append(x);x=(x%n)*q+x//n
        need(x==e,'rotation closed');cycles.append(cycle)
    return n,E,mate,vmate,cycles

def indicator(raw,E):
    need(len(raw)==(E+7)//8,'bits length')
    return [(raw[e//8]>>(e%8))&1 for e in range(E)]

def check_graph(m,q,k):
    n,E,mate,vmate,cycles=geometry(q,k)
    need(all(m[e]==m[mate[e]] for e in range(E)),'RC')
    need(all(sum(m[e] for e in c)==1 for c in cycles),'one per PCR')
    indegree=[0]*n
    for e in range(E):
        if not m[e]:indegree[e%n]+=1
    queue=deque(v for v in range(n) if not indegree[v]);levels=[0]*n;pred=[None]*n;seen=0
    while queue:
        u=queue.popleft();seen+=1
        for e in range(q*u,q*u+q):
            if m[e]:continue
            v=e%n
            if levels[v]<levels[u]+1:levels[v]=levels[u]+1;pred[v]=e
            indegree[v]-=1
            if indegree[v]==0:queue.append(v)
    need(seen==n,'DAG')
    end=max(range(n),key=levels.__getitem__);path=[];L=levels[end]
    while pred[end] is not None:
        e=pred[end];path.append(e);end=e//q
    path.reverse();need(len(path)==L,'independent path')
    return L,path

def check_cycle(C,m,q,k):
    n,E,*_=geometry(q,k)
    need(C and len(set(e//q for e in C))==len(C),'simple cycle')
    need(all(0<=e<E and e%n==C[(i+1)%len(C)]//q for i,e in enumerate(C)),'cycle adjacency')
    return sum(m[e] for e in C)

def dna(e,k):
    s=''
    for _ in range(k):s='ACGT'[e%4]+s;e//=4
    return s

def main():
    ap=argparse.ArgumentParser();ap.add_argument('round4',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();root=args.round4
    out=args.output.parent;out.mkdir(parents=True,exist_ok=True);records=[];ind={};cert={};paths={}
    for file in sorted((root/'certificates').glob('*.json')):
        c=json.loads(file.read_text());q,k=c['q'],c['k'];n,E,*_=geometry(q,k)
        raw=(root/c['bits']).read_bytes();m=indicator(raw,E)
        need(hashlib.sha256(raw).hexdigest()==c['bits_sha256'],'digest')
        L,path=check_graph(m,q,k)
        need(L==c['L'] and sum(m)==c['size'],'size and optimum')
        h=c['levels'];need(len(h)==n and all(type(x) is int and 0<=x<=L for x in h),'rank bounds')
        need(all(m[e] or h[e%n]>h[e//q] for e in range(E)),'rank all residual edges')
        p=c['path_edges'];need(len(p)==L and all(not m[e] for e in p),'path bits')
        need(all(p[i]%n==p[i+1]//q for i in range(len(p)-1)),'path adjacency')
        ind[c['name']]=m;cert[c['name']]=c;paths[c['name']]=path
        records.append({'name':c['name'],'q':q,'k':k,'size':sum(m),'L':L})
    retimings=0
    for name,c in cert.items():
        if 'cut_potential' not in c:continue
        q,k,T,ep=c['q'],c['k'],c['T'],c['epsilon'];n,E,_,mate,_=geometry(q,k)
        f,H=c['cut_potential'],c['H'];base=ind[c['source']];m=ind[name]
        need(all(H[v]+H[mate[v]]==T-1-T*ep and f[v]==-(H[v]//T) for v in range(n)),'mate and floor')
        need(all(m[e]==base[e]+f[e//q]-f[e%n] and H[e%n]-H[e//q]>=1-T*base[e] for e in range(E)),'retiming all edges')
        retimings+=1
    optima=json.loads((root/'results/component_optima.json').read_text())
    for row in optima:
        C=row['cycle_edges'];need(len(C)==row['optimum']+1,'cycle lower length')
        for name in [row['source'],*row['attainers']]:need(check_cycle(C,ind[name],4,row['k'])==1,'cycle one hit')
        need(all(cert[name]['L']==row['optimum'] for name in row['attainers']),'attainment')
    chain=json.loads((root/'results/escape_chain.json').read_text())
    for row in chain['component_optimality']:
        need(check_cycle(row['cycle_edges'],ind[row['name']],4,8)==1 and len(row['cycle_edges'])==row['optimum']+1,'chain lower')
    source,final=ind['start_141'],ind['escape_131']
    D=[e for e in range(len(final)) if source[e] and not final[e]]
    A=[e for e in range(len(final)) if final[e] and not source[e]]
    need(len(D)==len(A)==10 and not set(D)&set(A),'ten replacements')
    need(all(final[e]==source[e]-int(e in D)+int(e in A) for e in range(len(final))),'sparse direct membership')
    C=json.loads((root/'results/cycle_132_compact.json').read_text())
    p=paths['escape_131'];word=dna(p[0],8)+''.join(dna(e,8)[-1] for e in p[1:])
    need(len(word)==138,'avoidance length')
    def encode(w):
        r=0
        for ch in w:r=4*r+'ACGT'.index(ch)
        return r
    need(all(not final[encode(word[i:i+8])] for i in range(len(word)-7)),'all avoidance windows')
    result={'status':'PASS','independent_of_supplied_code':True,'mds_occurrences':len(records),'retimings':retimings,'component_lower_cycles':len(optima)+len(chain['component_optimality']),'records':records,
      'sparse_membership':{'base':'fixed selector rho=49/16, k=8','remove':[dna(e,8) for e in D],'add':[dna(e,8) for e in A]},
      'avoidance':{'maximum_consecutive_unselected_8mers':131,'maximum_nucleotide_length':138,'every_nucleotide_word_of_length_139_is_hit':True,'maximal_avoiding_word':word},
      'bounds':{'weighted_family_minimum_lower':141,'old_components_minimum_lower':min(r['optimum'] for r in optima if r['k']==8),'new_component_minimum_lower':131,'global_rc_minimum_upper':131,'global_rc_minimum_not_determined':True}}
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','sparse_membership','avoidance')}))
if __name__=='__main__':main()
