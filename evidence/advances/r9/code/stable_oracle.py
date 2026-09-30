#!/usr/bin/env python3
"""Exact stable-PCR completion with frozen nonstable selectors.
Python standard library only. Input: a Boolean RC PCR-transversal of 8230 words.
Example: python code/stable_oracle.py inputs/M107.words.txt --target 107 --out /tmp/oracle
A negative output concerns this frozen nonstable assignment, not the global problem.
"""
from verify_round9 import *

def solve_2sat(n,clauses):
    """Literal i^1 is the negation of literal i; clauses are disjunctions.
    Returns one true literal for each variable, or None. Iterative Kosaraju.
    """
    adj=[set() for _ in range(2*n)];rev=[set() for _ in range(2*n)]
    for a,b in clauses:
        req(0<=a<2*n and 0<=b<2*n,'SAT_LITERAL_RANGE')
        for u,v in [(a^1,b),(b^1,a)]:adj[u].add(v);rev[v].add(u)
    seen=set();order=[]
    for root in range(2*n):
        if root in seen:continue
        seen.add(root);stack=[(root,iter(sorted(adj[root])))]
        while stack:
            u,it=stack[-1]
            try:v=next(it)
            except StopIteration:order.append(u);stack.pop();continue
            if v not in seen:seen.add(v);stack.append((v,iter(sorted(adj[v]))))
    comp=[-1]*(2*n);k=0
    for root in reversed(order):
        if comp[root]>=0:continue
        comp[root]=k;stack=[root]
        while stack:
            u=stack.pop()
            for v in rev[u]:
                if comp[v]<0:comp[v]=k;stack.append(v)
        k+=1
    if any(comp[i]==comp[i^1] for i in range(0,2*n,2)):return None
    chosen={i if comp[i]>comp[i^1] else i^1 for i in range(0,2*n,2)}
    req(all(a in chosen or b in chosen for a,b in clauses),'SAT_ASSIGNMENT_CHECK')
    return chosen

def complete(base,L):
    stable=necklaces(base);pairs=[sorted(w for w in p if rc(w)==w) for p in stable];fixed=sum(pairs,[]);D=EDGESET-base-set(fixed)
    adj,order,a,b,path=graph_data(D)
    if len(path)>L:return None,{'status':'FIXED_PATH_OBSTRUCTION','target':L,'walk':path[:L+1]}
    mask={v:0 for v in V}
    for i,w in enumerate(fixed):mask[w[:-1]]|=1<<i
    for u in reversed(order):
        for w in adj[u]:mask[u]|=mask[w[1:]]
    clauses=[]
    for i,w in enumerate(fixed):
        if a[w[:-1]]+1+b[w[1:]]>L or mask[w[1:]]>>i&1:clauses.append((i^1,i^1))
        for j in range(i):
            if mask[w[1:]]>>j&1:clauses.append((i^1,j^1))
    chosen=solve_2sat(128,clauses)
    if chosen is None:return None,{'status':'TWO_SAT_UNSAT','target':L,'literal_words':fixed,'clauses':clauses,'boundary':'exact SCC decision; use the included clauses to replay, not as a global nonexistence claim'}
    out=(base-set(fixed))|{w for i,w in enumerate(fixed) if i not in chosen};necklaces(out);_,_,d,_,p=graph_data(EDGESET-out);req(len(p)<=L,'ORACLE_ENDPOINT_LENGTH')
    return out,{'status':'SAT','target':L,'L':len(p),'restored_words':[fixed[i] for i in sorted(chosen)],'rank':[d[v] for v in V],'longest_path':p}

def regression():
    from itertools import combinations,product
    count=0
    for n,limit in [(2,8),(3,4)]:
        possibilities=[(i,i) for i in range(2*n)]+[(i,j) for i in range(2*n) for j in range(i+1,2*n) if i//2!=j//2]
        for k in range(limit+1):
            for clauses in combinations(possibilities,k):
                truth=any(all(a in ss or b in ss for a,b in clauses) for bits in product([0,1],repeat=n) for ss in [{2*i+bits[i] for i in range(n)}])
                result=solve_2sat(n,clauses);req((result is not None)==truth,'SCC_REGRESSION');count+=1
    return {'formulas_checked':count,'scope':'all 2-variable non-tautological 2-CNF formulas; all 3-variable formulas with at most 4 distinct non-tautological clauses'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path,nargs='?');ap.add_argument('--target',type=int,default=106);ap.add_argument('--out',type=Path);ap.add_argument('--self-test',action='store_true');args=ap.parse_args()
    if args.self_test:
        out={'status':'PASS',**regression()};(OUT/'oracle_regression.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));return
    req(args.input is not None and args.out is not None and args.target>=0,'ORACLE_ARGUMENTS');args.out.mkdir(parents=True,exist_ok=True)
    base=read_words(args.input);m,record=complete(base,args.target);(args.out/'result.json').write_text(json.dumps(record,indent=2))
    if m is not None:(args.out/'endpoint.words.txt').write_text(''.join(w+'\n' for w in sorted(m)))
    print(record['status'])
if __name__=='__main__':main()
