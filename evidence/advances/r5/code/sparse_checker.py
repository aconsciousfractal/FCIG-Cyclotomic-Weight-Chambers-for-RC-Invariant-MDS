#!/usr/bin/env python3
"""Independent exact verification of sparse-neighborhood cut/branch refutations.
Standard library only. Literal words and complete RC-PCR domains; no solver,
NumPy, Numba, or discovery modules. Local scope is enforced by the radius.
"""
from collections import Counter, defaultdict
from pathlib import Path
import json, sys, time
import graph_checker as g
from radius_checker import all_domains

def verify_sparse(base, cert):
    g.require(cert['schema']=='rc-mds-sparse-neighborhood-refutation-v1','SPARSE_SCHEMA')
    g.require(cert['status']=='EXACT_EXCLUSION','SPARSE_UNFINISHED')
    g.require(cert['baseline_sha256']==g.digest_set(base),'SPARSE_BASE_DIGEST')
    B,R=cert['bound'],cert['radius']
    g.require(type(B) is int and B>=0 and type(R) is int and R>=0,'SPARSE_PARAMETERS')
    groups=all_domains(base)
    alternatives=[]; inverse=defaultdict(list)
    for i,gr in enumerate(groups):
        old=set(gr['options'][gr['base_option']])
        for j,op in enumerate(gr['options']):
            if j==gr['base_option']: continue
            new=set(op); delta={w:int(w in new)-int(w in old) for w in old|new}
            a=len(alternatives);alternatives.append((i,j))
            for w,d in delta.items():
                if d: inverse[w].append((a,d))
    changes=[];slacks=[];positives=[]; kinds=Counter()
    for cut in cert['cuts']:
        kind=cut['kind'];g.require(kind in ('cycle','path'),'SPARSE_CUT_KIND')
        words=cut['words'];_,_,hit=g.walk_check(words,base,cycle=kind=='cycle')
        if kind=='path':
            g.require(len(words)==B+1,'SPARSE_PATH_LENGTH');rhs=1
        else:rhs=(len(words)+B)//(B+1)
        g.require(type(cut['rhs']) is int and cut['rhs']==rhs,'SPARSE_CUT_RHS')
        coeff=Counter()
        for w,c in Counter(words).items():
            for a,d in inverse[w]:coeff[a]+=c*d
        changes.append(coeff);slacks.append(hit-rhs)
        positives.append(sorted(a for a,d in coeff.items() if d>0));kinds[kind]+=1
    nodes=cert['nodes'];root=cert['root']
    g.require(type(root) is int and 0<=root<len(nodes),'SPARSE_ROOT')
    expected={root:()};queue=[root];types=Counter();depths=Counter();branch_edges=0
    for ni in queue:
        node=nodes[ni];g.require(isinstance(node,dict),'SPARSE_MISSING_NODE')
        chosen=expected[ni]
        g.require(node['chosen']==list(chosen),'SPARSE_NODE_ASSIGNMENT')
        g.require(all(type(a) is int and 0<=a<len(alternatives) for a in chosen),'SPARSE_OPTION_RANGE')
        used={alternatives[a][0] for a in chosen}
        g.require(len(used)==len(chosen) and len(chosen)<=R,'SPARSE_GROUP_OR_RADIUS')
        c=node['cut'];g.require(type(c) is int and 0<=c<len(changes),'SPARSE_CUT_INDEX')
        row=changes[c];slack=slacks[c]+sum(row[a] for a in chosen)
        g.require(slack<0,'SPARSE_CUT_NOT_VIOLATED')
        available=[a for a in positives[c] if alternatives[a][0] not in used]
        kind=node['type'];types[kind]+=1;depths[len(chosen)]+=1
        if kind=='capacity_leaf':
            maxima={}
            for a in available:
                gr=alternatives[a][0];maxima[gr]=max(maxima.get(gr,0),row[a])
            remaining=R-len(chosen)
            cap=sum(sorted(maxima.values(),reverse=True)[:remaining])
            g.require(all(type(node[x]) is int for x in ('slack','remaining','capacity')),'SPARSE_LEAF_TYPES')
            g.require(node['slack']==slack and node['remaining']==remaining and node['capacity']==cap,'SPARSE_LEAF_ACCOUNTING')
            g.require(slack+cap<0,'SPARSE_LEAF_INEQUALITY')
        elif kind=='branch':
            children=node['children']
            g.require(len(chosen)<R,'SPARSE_BRANCH_AT_LIMIT')
            g.require([p[0] for p in children]==available,'SPARSE_BRANCH_COVERAGE')
            for a,child in children:
                g.require(type(a) is int and type(child) is int and 0<=child<len(nodes),'SPARSE_CHILD_RANGE')
                nxt=tuple(sorted(chosen+(a,)))
                if child in expected:g.require(expected[child]==nxt,'SPARSE_BAD_MERGE')
                else:expected[child]=nxt;queue.append(child)
                branch_edges+=1
        else:raise g.VerificationError('SPARSE_NODE_KIND')
    g.require(len(expected)==len(nodes),'SPARSE_UNREACHABLE_RECORDS')
    sizes=[len(gr['options'])-1 for gr in groups];counts=[1]+[0]*R
    for x in sizes:
        for j in range(R,0,-1): counts[j]+=counts[j-1]*x
    return {'status':'PASS','bound':B,'radius':R,'RC_PCR_groups':len(groups),
      'alternatives':len(alternatives),'neighborhood_counts_by_exact_distance':counts,
      'assignments_covered':sum(counts),'nodes':len(nodes),'branch_edges':branch_edges,
      'node_types':dict(types),'nodes_by_depth':dict(sorted(depths.items())),
      'cuts':len(changes),'cut_types':dict(kinds),'global_nonexistence':False}

if __name__=='__main__':
    t=time.time();base=g.read_set(Path(sys.argv[1]));cert=json.loads(Path(sys.argv[2]).read_text())
    result=verify_sparse(base,cert);result['elapsed_seconds']=time.time()-t
    print(json.dumps(result,indent=2))
