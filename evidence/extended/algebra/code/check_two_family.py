"""Check the new two-family selector and three-block path formula on all 30 rows.

Reuses the original producer's exact sign arithmetic and deterministic zero
PCR tie rule. It does NOT call its seed reconstruction, terminal construction,
legal firing routine, or residual path routine. Hence this tests the new
algebraic reduction; the original separate validator supplies arithmetic
cross-assurance. The lower graph, not its line graph, is used for edge counts.
"""
from __future__ import annotations
from collections import Counter,deque
from pathlib import Path
import importlib.util,json,time

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('source_producer',ROOT.parents[1]/'finite_assurance/scripts/produce.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

def dag_longest(adj,allowed=None):
    n=len(adj); vertices=list(range(n)) if allowed is None else list(allowed); S=set(vertices)
    indeg={u:0 for u in vertices}
    for u in vertices:
        for v in adj[u]:
            if v in S: indeg[v]+=1
    Q=deque(u for u in vertices if indeg[u]==0); dist={u:0 for u in vertices};order=[]
    while Q:
        u=Q.popleft();order.append(u)
        for v in adj[u]:
            if v not in S: continue
            dist[v]=max(dist[v],dist[u]+1);indeg[v]-=1
            if not indeg[v]: Q.append(v)
    if len(order)!=len(vertices): raise AssertionError('cycle in reconstructed graph')
    return max(dist.values(),default=0),dist,order

def main():
    start=time.time();original=json.loads((ROOT/'results/original_producer.json').read_text())['scientific_payload']
    byid={row['row_id']:row for row in original['rows']}
    skeleton=json.loads((ROOT.parents[1]/'finite_assurance/inputs/critical_ratio_skeleton.json').read_text())['scientific_payload']
    rows=[];cache={};checks=0;decisions=0
    for spec in p.row_specs(skeleton):
        k=spec['k'];n=4**(k-1);reference=byid[spec['row_id']];ratio=p.ExactRatio(spec['weight'])
        if k not in cache:
            words=list(p.words(k-1))
            lforms=[p.expression(w,p.family_coefficients(k,'LOWER_LAMBDA')) for w in words]
            rforms=[p.expression(w,p.family_coefficients(k,'LOWER_MIDPOINT')) for w in words]
            reps={};maps={};struct=0
            for code,w in enumerate(words):
                if lforms[code]==(0,0,0,0): struct+=1;continue
                c=tuple(p.COMPLEMENT[x] for x in w);r=tuple(reversed(w));cr=p.rc(w)
                orbit={w:(1,1),c:(-1,-1),r:(-1,1),cr:(1,-1)}
                if len(orbit)!=4: raise AssertionError('nonfree active orbit')
                rep=min(orbit); reps[p.word_code(rep)]=None
                # The four involutions have their own inverses, hence the same signs.
                maps[code]=(p.word_code(rep),*orbit[rep])
            if len(reps)*4!=n-struct: raise AssertionError('orbit count')
            cache[k]=(words,lforms,rforms,reps,maps,struct)
        words,lforms,rforms,reps,maps,struct=cache[k]
        sl=[0]*n; sr=[0]*n
        evaluated={i:(ratio.sign_expression(lforms[i]),ratio.sign_expression(rforms[i])) for i in reps}
        for code,(rep,ls,rs) in maps.items(): sl[code]=ls*evaluated[rep][0];sr[code]=rs*evaluated[rep][1]
        # Check the restoration independently against all directly evaluated forms.
        if any(sl[i]!=ratio.sign_expression(lforms[i]) for i in range(n)): raise AssertionError('lambda sign restoration')
        if any(sr[i]!=ratio.sign_expression(rforms[i]) for i in maps): raise AssertionError('midpoint sign restoration')
        a=[int(sl[i]>0 and (sr[i]<0 or (sr[i]==0 and words[i]<p.rc(words[i])))) for i in range(n)]
        colour=[0 if sl[i]<=0 else (1 if a[i] else -1) for i in range(n)]
        terminal=[];zero_count=0;adj=[[] for _ in range(n)];seed=[]
        for code in range(4*n):
            u,v=code//4,code%n
            if sl[u]==sl[v]==0:
                zero_count+=1;m0=p.zero_seed_membership(p.decode_word(code,k))
            else: m0=int(sl[u]<=0<sl[v])
            m=m0+a[u]-a[v]
            if m not in (0,1): raise AssertionError('nonboolean two-family selector')
            if colour[u]==-1 and colour[v]==1: raise AssertionError('forbidden P0 to P1 edge')
            if m==0:
                if colour[u]>colour[v]: raise AssertionError('descending colour on residual edge')
                adj[u].append(v)
            terminal.append(m);seed.append(m0)
        selected=[i for i,x in enumerate(terminal) if x]
        if p.digest(selected)!=reference['terminal']['selected_word_sha256']: raise AssertionError('terminal mismatch')
        if p.digest([i for i,x in enumerate(a) if x])!=reference['terminal']['a_star_sha256']: raise AssertionError('ideal mismatch')
        if zero_count!=reference['zero_and_descent']['zero_embedding_word_count']: raise AssertionError('zero mismatch')
        lengths={};groups={c:{i for i,x in enumerate(colour) if x==c} for c in (-1,0,1)}
        distances={};orders={}
        for c in (-1,0,1): lengths[c],distances[c],orders[c]=dag_longest(adj,groups[c])
        L,_,_=dag_longest(adj)
        if L!=reference['residual_dag']['longest_path_edges']: raise AssertionError('path mismatch')
        if lengths[-1]!=lengths[1]: raise AssertionError('RC block mismatch')
        alpha={v:0 for v in groups[0]}
        for u in groups[-1]:
            for v in adj[u]:
                if v in alpha: alpha[v]=max(alpha[v],distances[-1][u]+1)
        beta=dict(alpha)
        for u in orders[0]:
            for v in adj[u]:
                if v in beta: beta[v]=max(beta[v],beta[u]+1)
        maxplus=max([lengths[-1],*(beta[v]+alpha[p.word_code(p.rc(words[v]))] for v in groups[0])])
        if maxplus!=L: raise AssertionError('max-plus interface mismatch')
        # Direct check of the line-graph normalization at k=4 (all three rows).
        upper_path=None
        if k==4:
            remaining=[i for i,x in enumerate(terminal) if not x];ids={x:i for i,x in enumerate(remaining)}
            upperadj=[[] for _ in remaining]
            for i,code in enumerate(remaining):
                for ch in range(4):
                    nxt=(code%n)*4+ch
                    if nxt in ids: upperadj[i].append(ids[nxt])
            upper_path=dag_longest(upperadj)[0]
            if upper_path!=L-1: raise AssertionError('line-graph normalization')
        rows.append({'row_id':spec['row_id'],'k':k,'weight':spec['weight'],'matches_original':True,
                     'selected_count':len(selected),'lower_states':n,'global_structural_lambda_states':struct,
                     'symmetry_reduced_form_count':2*len(reps),'lower_graph_longest_edges':L,
                     'positive_block_longest_edges':lengths[-1],'nonpositive_core_longest_edges':lengths[0],
                     'block_upper_bound':2*lengths[-1]+lengths[0]+2,'maxplus_exact':maxplus,
                     'upper_graph_longest_edges_direct_k4':upper_path})
        checks+=10+(k==4);decisions+=4*n
    # Minimal root locations: keep a point whenever one of its three local outputs differs.
    essential={}
    for order in original['orders']:
        if not order['complete_declared_range']: continue
        ids=order['ordered_row_ids'];rootdata=[]
        for j in range(1,len(ids),2):
            triple=[byid[x]['terminal']['selected_word_sha256'] for x in ids[j-1:j+2]]
            point=next(x for o in skeleton['orders'] if o['k']==order['k'] for x in o['critical_points'] if 'ATLAS-'+x['cell_id']==ids[j])
            rootdata.append({'ratio':point['ratio'],'essential':len(set(triple))>1,
                             'point_only_anomaly':triple[0]==triple[2]!=triple[1],
                             'three_different':len(set(triple))==3})
        E=sum(x['essential'] for x in rootdata)
        essential[str(order['k'])]={'roots':rootdata,'essential_count':E,'minimum_root_partition_cells':2*E+1,
                                    'distinct_terminal_classes':order['unique_terminal_geometry_count']}
    output={'status':'PASS','scope':'new formulas, source exact sign oracle reused','rows_checked':len(rows),
            'selected_word_decisions_checked':decisions,'checks':checks,'seconds':time.time()-start,
            'rows':rows,'minimal_one_parameter_root_locations':essential}
    (ROOT/'results/two_family_and_residual.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(rows),'rows;',decisions,'word decisions;',checks,'checks;',round(time.time()-start,2),'seconds')

if __name__=='__main__': main()
