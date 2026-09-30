"""Cut potentials, short directed-cycle separators and RC-paired F geodesics.
A positive certificate is an explicit sequence of legal F_u / RF_RC(u)
pairs. A negative certificate is a directed cycle with unequal hit counts.
"""
from core import *
from heapq import heapify,heappop,heappush
from itertools import combinations

def tree_potential(m,n,q=4):
    p=[0]*n
    # Root is 0^(k-1); parent of v>0 is v//q, edge-code v.
    for v in range(1,n):p[v]=p[v//q]+m[v]
    return p

def cut_normal(m,n,q=4):
    p=tree_potential(m,n,q)
    return tuple(m[e]+p[e//q]-p[e%n] for e in range(len(m))),p

def potential_difference(m1,m2,n,q=4):
    p1=tree_potential(m1,n,q);p2=tree_potential(m2,n,q)
    f=[p1[u]-p2[u] for u in range(n)]
    check(all(m2[e]-m1[e]==f[e//q]-f[e%n] for e in range(q*n)),'NOT_A_CUT')
    return f

def pair_geodesic(m1,m2,k,q=4):
    n=q**(k-1);theta=[rc(u,k-1,q) for u in range(n)];f=potential_difference(m1,m2,n,q)
    C=f[0]+f[theta[0]]
    check(all(f[u]+f[theta[u]]==C for u in range(n)),'RC_CUT_CONSTANT')
    check(C%2==0,'ODD_RC_PARITY')
    g=[x-C//2 for x in f];distance=sum(abs(x) for x in g)//2
    original_g=g[:];m=bytearray(m1);moves=[]
    while max(g)>0:
        h=max(g);A={u for u in range(n) if g[u]==h};ind={u:0 for u in A}
        for v in A:ind[v]=sum(not m[e] and e//q in A for e in inedges(v,n,q))
        ready=[u for u in A if ind[u]==0];heapify(ready)
        while ready:
            u=heappop(ready)
            if u not in A:continue
            v=theta[u]
            iu=inedges(u,n,q);ou=list(outedges(u,n,q));iv=inedges(v,n,q);ov=list(outedges(v,n,q))
            check(not (set(iu+ou)&set(iv+ov)),'PAIRED_STARS_OVERLAP')
            check(all(m[e]==1 for e in iu),'F_NOT_LEGAL')
            check(all(m[e]==1 for e in ov),'RF_NOT_LEGAL')
            # Save source-level outgoing residual neighbors before the flip.
            affected=[e%n for e in ou if not m[e] and e%n in A and e%n!=u]
            delta={}
            for e in ou:delta[e]=delta.get(e,0)+1
            for e in iu:delta[e]=delta.get(e,0)-1
            for e in ov:delta[e]=delta.get(e,0)-1
            for e in iv:delta[e]=delta.get(e,0)+1
            for e,d in delta.items():
                check(m[e]+d in (0,1),'PAIR_NONBOOLEAN');m[e]+=d
            A.remove(u);g[u]-=1;g[v]+=1;moves.append(u)
            for z in affected:
                if z in A:
                    ind[z]-=1
                    if ind[z]==0:heappush(ready,z)
        check(not A,'MAX_LEVEL_SOURCE_BLOCKED')
    check(all(x==0 for x in g),'PAIRED_POTENTIAL_NOT_EXHAUSTED')
    check(m==m2,'GEODESIC_ENDPOINT');check(len(moves)==distance,'GEODESIC_NOT_SHARP')
    return {'cut_constant':C,'max_abs_antisymmetric_potential':max(map(abs,original_g),default=0),
            'minimum_paired_moves':distance,'moves':moves}

def replay_pair(m1,m2,k,record,q=4):
    n=q**(k-1);f=potential_difference(m1,m2,n,q);C=f[0]+f[rc(0,k-1,q)]
    check(C==record['cut_constant'] and C%2==0,'PAIR_CERT_PARITY')
    g=[x-C//2 for x in f];lower=sum(abs(x) for x in g)//2
    check(len(record['moves'])==record['minimum_paired_moves']==lower,'PAIR_CERT_DISTANCE')
    m=bytearray(m1)
    for u in record['moves']:
        check(0<=u<n,'PAIR_VERTEX_RANGE');v=rc(u,k-1,q)
        iu=set(inedges(u,n,q));ou=set(outedges(u,n,q));iv=set(inedges(v,n,q));ov=set(outedges(v,n,q))
        check(not ((iu|ou)&(iv|ov)),'PAIR_CERT_OVERLAP')
        check(all(m[e]==1 for e in iu) and all(m[e]==1 for e in ov),'PAIR_CERT_LEGALITY')
        # Set-operation implementation, distinct from additive producer.
        for e in (iu-ou)|(ov-iv):m[e]=0
        for e in (ou-iu)|(iv-ov):
            check(m[e]==0,'PAIR_CERT_OVERSELECTED');m[e]=1
    check(m==m2,'PAIR_CERT_ENDPOINT')

def outward(v,q=4):
    es=[]
    while v:es.append(v);v//=q
    return es[::-1]
def homeward(v,n,q=4):
    es=[]
    while v:
        e=q*v;es.append(e);v=e%n
    return es

def simple_nonzero_cycle(es,d,n,q=4):
    check(bool(es),'EMPTY_CLOSED_WALK')
    start=es[0]//q;vertices=[start];stack=[];pos={start:0}
    for e in es:
        check(e//q==vertices[-1],'WALK_DISCONTINUITY')
        v=e%n
        if v in pos:
            i=pos[v];cy=stack[i:]+[e]
            if sum(d[x] for x in cy):return cy
            for w in vertices[i+1:]:del pos[w]
            vertices=vertices[:i+1];stack=stack[:i]
        else:
            stack.append(e);pos[v]=len(vertices);vertices.append(v)
    raise CheckError('NO_NONZERO_SIMPLE_CYCLE')

def separator(m1,m2,k,q=4):
    n=q**(k-1);d=[b-a for a,b in zip(m1,m2)];p=tree_potential(d,n,q)
    for u in range(n):
        es=outward(u,q)+homeward(u,n,q)
        if sum(d[e] for e in es):break
    else:
        e=next((e for e in range(q*n) if d[e]+p[e//q]-p[e%n]),None)
        check(e is not None,'SAME_COMPONENT_NO_SEPARATOR')
        es=outward(e//q,q)+[e]+homeward(e%n,n,q)
    check(len(es)<=2*k-1,'SEPARATOR_WALK_TOO_LONG')
    cy=simple_nonzero_cycle(es,d,n,q)
    result={'edge_codes':cy,'words':[text(e,k,q) for e in cy],
            'length':len(cy),'hits_source':sum(m1[e] for e in cy),'hits_target':sum(m2[e] for e in cy)}
    verify_separator(m1,m2,k,result,q);return result

def verify_separator(m1,m2,k,rec,q=4):
    n=q**(k-1);es=rec['edge_codes']
    check(0<len(es)==rec['length']<=2*k-1,'SEPARATOR_LENGTH')
    check(all(0<=e<q*n for e in es),'SEPARATOR_EDGE_RANGE')
    check(len({e//q for e in es})==len(es),'SEPARATOR_NOT_SIMPLE')
    check(all(es[i]%n==es[(i+1)%len(es)]//q for i in range(len(es))),'SEPARATOR_NOT_CLOSED')
    h1=sum(m1[e] for e in es);h2=sum(m2[e] for e in es)
    check(h1==rec['hits_source'] and h2==rec['hits_target'] and h1!=h2,'SEPARATOR_HITS')

def run(parent,out):
    all_summary=[]
    for k in (4,6,8):
        a,ts=parent_data(parent,k);n=4**(k-1);norms={};classes=[];index=[]
        for r,t in zip(a['rows'],ts):
            norm,_=cut_normal(t,n)
            if norm not in norms:norms[norm]=len(classes);classes.append([])
            z=norms[norm];classes[z].append(len(index));index.append(z)
        positives=[]
        for ids in classes:
            source=ids[0]
            for target in ids:
                path=pair_geodesic(ts[source],ts[target],k)
                replay_pair(ts[source],ts[target],k,path)
                positives.append({'source':a['rows'][source]['id'],'target':a['rows'][target]['id'],**path})
        negatives=[]
        for i,j in combinations(range(len(classes)),2):
            u,v=classes[i][0],classes[j][0]
            negatives.append({'source_class':i,'target_class':j,
                              'source':a['rows'][u]['id'],'target':a['rows'][v]['id'],
                              **separator(ts[u],ts[v],k)})
        reps=[a['rows'][ids[0]] for ids in classes]
        res={'k':k,'component_count':len(classes),'row_count':len(index),
             'classes':[[a['rows'][j]['id'] for j in ids] for ids in classes],
             'representative_ratios':[r['ratio'] for r in reps],
             'row_component':{r['id']:c for r,c in zip(a['rows'],index)},
             'positive_geodesics':positives,'negative_cycles':negatives}
        dump(out/f'components_k{k}.json',res)
        summary={'k':k,'components':len(classes),'strata':len(index),'distinct_literal_terminals':a['distinct_terminals'],
                  'positive_certificates':len(positives),'paired_moves_replayed':sum(r['minimum_paired_moves'] for r in positives),
                  'negative_cycles':len(negatives),'maximum_separator_length':max((r['length'] for r in negatives),default=0),
                  'max_geodesic_length':max(r['minimum_paired_moves'] for r in positives),
                  'class_stratum_sizes':[len(ids) for ids in classes]}
        all_summary.append(summary);print('COMPONENTS_OK',summary,flush=True)
    dump(out/'component_summary.json',{'status':'PASS','rows':all_summary})
    return all_summary

def binary_example(out):
    q=2;k=4;n=q**(k-1);unused=set(range(q*n));cy=[]
    while unused:
        e=min(unused);c=[];v=e
        while v not in c:c.append(v);v=(v%n)*q+v//n
        unused.difference_update(c);cy.append(c)
    results=[];scanned=0
    for choices in product(*cy):
        scanned+=1;m=bytearray(int(e in choices) for e in range(q*n))
        if any(m[e]!=m[rc(e,k,q)] for e in range(q*n)):continue
        try:dag(m,n,q)
        except CheckError:continue
        results.append(m)
    groups={}
    for i,m in enumerate(results):groups.setdefault(cut_normal(m,n,q)[0],[]).append(i)
    odd=[]
    for ids in groups.values():
        for i,j in combinations(ids,2):
            f=potential_difference(results[i],results[j],n,q);C=f[0]+f[rc(0,k-1,q)]
            check(C%2==1,'BINARY_EXPECTED_ODD_PARITY')
            odd.append({'source':i,'target':j,'potential':f,'constant':C})
    res={'q':q,'k':k,'pcr_transversals_scanned':scanned,'rc_mds_count':len(results),
         'F_components':[v for v in groups.values()],
         'selected_edge_codes':[[e for e,x in enumerate(m) if x] for m in results],
         'selected_words':[[text(e,k,q) for e,x in enumerate(m) if x] for m in results],
         'odd_parity_pairs':odd,
         'negative_cycles':[{'source':i[0],'target':j[0],**separator(results[i[0]],results[j[0]],k,q)} for i,j in combinations(groups.values(),2)]}
    check(scanned==128 and len(results)==6 and len(groups)==3 and len(odd)==3,'BINARY_ENUMERATION')
    dump(out/'binary_parity_obstruction.json',res);return res
