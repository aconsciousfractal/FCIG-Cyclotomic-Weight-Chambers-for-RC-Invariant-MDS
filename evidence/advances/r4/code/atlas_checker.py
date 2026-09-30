"""Independent exact checker of the complete tight simple-cycle census.
Python standard library; no NetworkX, Numba, search code or cycle generator.
All cycles are counted up to cyclic starting position (not modulo RC).
"""
from collections import Counter, deque
import graph_checker as g


def strong_components(vertices, edges):
    adj={v:[] for v in vertices};rev={v:[] for v in vertices}
    for e in edges:adj[e[:-1]].append(e[1:]);rev[e[1:]].append(e[:-1])
    seen=set();order=[]
    for root in vertices:
        if root in seen:continue
        seen.add(root);stack=[(root,0)]
        while stack:
            v,j=stack[-1]
            if j==len(adj[v]):order.append(v);stack.pop();continue
            u=adj[v][j];stack[-1]=(v,j+1)
            if u not in seen:seen.add(u);stack.append((u,0))
    used=set();comps=[]
    for root in reversed(order):
        if root in used:continue
        used.add(root);c=[];stack=[root]
        while stack:
            v=stack.pop();c.append(v)
            for u in rev[v]:
                if u not in used:used.add(u);stack.append(u)
        if len(c)>1 or root in adj[root]:comps.append(sorted(c))
    return sorted(comps)


def component_counts(vertices, edges, anchor, T, target=None):
    adj={v:[] for v in vertices};rev={v:[] for v in vertices}
    for e in edges:adj[e[:-1]].append(e);rev[e[1:]].append(e)
    levels={anchor:0};q=deque([anchor])
    while q:
        u=q.popleft()
        for e in adj[u]:
            v=e[1:];ell=(levels[u]+1)%T
            if v not in levels:levels[v]=ell;q.append(v)
            else:g.require(levels[v]==ell,'ATLAS_CIRCULAR_LEVEL')
    g.require(set(levels)==set(vertices),'ATLAS_UNREACHED_VERTEX')
    g.require([v for v in vertices if levels[v]==0]==[anchor],'ATLAS_ANCHOR_NOT_UNIQUE')
    g.require(set(levels.values())==set(range(T)),'ATLAS_LEVEL_COVERAGE')
    order=sorted(vertices,key=lambda v:(levels[v],v))
    F=dict.fromkeys(vertices,0);F[anchor]=1
    for u in order:
        for e in adj[u]:
            if e[1:]!=anchor:F[e[1:]]+=F[u]
    total=sum(F[e[:-1]] for e in rev[anchor]);g.require(total>0,'ATLAS_NO_CYCLE')
    R=dict.fromkeys(vertices,0);R[anchor]=1
    for u in reversed(order):
        if u!=anchor:R[u]=sum(R[e[1:]] for e in adj[u])
    g.require(sum(R[e[1:]] for e in adj[anchor])==total,'ATLAS_TWO_COUNTS_DISAGREE')
    weights={e:F[e[:-1]]*R[e[1:]] for e in edges}
    g.require(all(type(n) is int and n>0 for n in weights.values()),'ATLAS_UNUSED_EDGE')
    g.require(sum(weights.values())==T*total,'ATLAS_INCIDENCE_MASS')
    result={'total':total,'weights':weights}
    if target is not None:
        polys={v:Counter() for v in vertices};polys[anchor][0]=1;Z=Counter()
        for u in order:
            for e in adj[u]:
                dest=Z if e[1:]==anchor else polys[e[1:]]
                for h,n in polys[u].items():dest[h+int(e in target)]+=n
        g.require(sum(Z.values())==total,'ATLAS_POLYNOMIAL_MASS')
        result['histogram']={str(h):n for h,n in sorted(Z.items())}
    return result


def verify_atlas(base, cert, atlas, target=None, expected_histogram=None):
    g.require(atlas['schema']=='rc-mds-tight-cycle-atlas-v1','ATLAS_SCHEMA')
    g.require(atlas['baseline_sha256']==g.digest_set(base),'ATLAS_BASE_DIGEST')
    g.require(type(atlas['L']) is int and type(atlas['tight_edge_count']) is int,'ATLAS_HEADER_TYPES')
    T=cert['L']+1;g.require(atlas['L']==T-1,'ATLAS_BASE_LENGTH')
    rank=dict(zip(g.VERTICES,cert['rank']))
    # The set certificate must already have been checked by the caller.
    edges=[e for e in g.EDGES if rank[e[1:]]-rank[e[:-1]]==1-T*int(e in base)]
    g.require(atlas['tight_edge_count']==len(edges),'ATLAS_TIGHT_EDGE_COUNT')
    comps=strong_components(g.VERTICES,edges)
    g.require(len(atlas['components'])==len(comps),'ATLAS_COMPONENT_COMPLETENESS')
    summaries=[];allhist=Counter()
    for vertices,record in zip(comps,atlas['components']):
        g.require(type(record['T']) is int and type(record['simple_cycle_count']) is int
          and type(record['baseline_weighted_sum']) is int
          and all(type(n) is int for n in record['edge_cycle_counts'].values()),'ATLAS_COUNT_TYPES')
        V=set(vertices);ce=[e for e in edges if e[:-1] in V and e[1:] in V]
        g.require(record['vertices']==vertices,'ATLAS_VERTEX_SET')
        g.require(record['edges']==ce,'ATLAS_EDGE_SET')
        g.require(record['T']==T,'ATLAS_PERIOD')
        levels={}
        for v in vertices:levels.setdefault(rank[v],[]).append(v)
        singles={str(r):vs[0] for r,vs in sorted(levels.items()) if len(vs)==1}
        g.require(record['singleton_levels']==singles and len(singles)>0,'ATLAS_SINGLETON_LEVELS')
        anchor=record['anchor'];g.require(anchor in singles.values(),'ATLAS_BAD_ANCHOR')
        counts=component_counts(vertices,ce,anchor,T,target)
        g.require(record['simple_cycle_count']==counts['total'],'ATLAS_CYCLE_COUNT')
        g.require(record['edge_cycle_counts']==counts['weights'],'ATLAS_EDGE_MULTIPLICITIES')
        selected=[e for e in ce if e in base]
        g.require(record['selected_edges']==selected,'ATLAS_SELECTED_EDGES')
        mass=sum(counts['weights'][e] for e in selected)
        g.require(record['baseline_weighted_sum']==mass==counts['total'],'ATLAS_BASE_WEIGHTED_SUM')
        cut=record['all_cycle_cut'];B=cut['bound']
        g.require(type(B) is int and B==T-2,'ATLAS_CUT_BOUND')
        g.require(type(cut['rhs']) is int,'ATLAS_RHS_TYPE')
        g.require(cut['rhs']==((T+B)//(B+1))*counts['total'],'ATLAS_CUT_RHS')
        summaries.append({'vertices':len(vertices),'edges':len(ce),'simple_cycles':counts['total'],
          'common_vertices':len(singles),'common_edges':sum(n==counts['total'] for n in counts['weights'].values()),
          'selected_edges':len(selected),'cycle_length':T})
        if target is not None:allhist.update(counts['histogram'])
    if target is not None and expected_histogram is not None:
        g.require(all(isinstance(k,str) and type(n) is int for k,n in expected_histogram.items()),'ATLAS_HISTOGRAM_TYPES')
        g.require(dict(allhist)==expected_histogram,'ATLAS_CROSS_HISTOGRAM')
    return {'status':'PASS','tight_edges':len(edges),'cyclic_components':summaries,
      'cross_histogram':dict(allhist) if target is not None else None}
