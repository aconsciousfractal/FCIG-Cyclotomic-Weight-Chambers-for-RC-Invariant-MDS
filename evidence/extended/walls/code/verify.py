"""Certificate checker. It does not import walls.py or its sparse producer.
With --parent-replay it also invokes the inherited independent original-seed
verifier on all 199 strata, without changing the frozen predecessor.
"""
from core import *
from components import replay_pair,verify_separator
import argparse,sys,time,importlib.util
from itertools import combinations

def verify_sparse_certificate(m1,m2,c):
    delta={e:b-a for e,(a,b) in enumerate(zip(m1,m2)) if b!=a}
    stored=dict(c['patch'])
    check(len(stored)==len(c['patch']),'PATCH_DUPLICATES')
    check(stored==delta,'PATCH_LITERAL_IDENTITY')
    check(all(x in (-1,1) for x in stored.values()),'PATCH_DELTA_RANGE')
    check(c['external_companion_pairs']==len(c['external_exchanges']),'EXTERNAL_COUNT')
    check(c['changed_null_pcrs']==len(c['changed_cycle_indices']),'CYCLE_COUNT')
    check(c['symmetric_difference']==len(delta)==8*c['regular_midpoint_states']+2*c['external_companion_pairs']+2*c['changed_null_pcrs'],'MIXED_COUNT_IDENTITY')

def verify_mixed(parent,out,k,A,ts,mod,models,cs):
    checked=0;words=0
    for rr in A['roots']:
        if 'lower_lambda' not in rr['families']:continue
        rec=json.loads((out/f'mixed_k{k}_r{rr["id"]:03d}.json').read_text());j=2*rr['id'];n=len(models['lw']);point=cs[j+1]
        check(rec['ratio']==rr['value'],'MIXED_RATIO_IDENTITY')
        Z={u for u,x in enumerate(point['sl']) if x==0}
        S={u for u in range(n) if point['sl'][u]>0 and point['sm'][u]==0 and models['mr'][u]!=(0,0,0,0)}
        check(S==set(rec['regular_midpoint_states']),'MIXED_MIDPOINT_LIST')
        check(len(Z)==rec['singular_vertices'],'MIXED_ZERO_VERTEX_COUNT')
        Zcross={u for u in Z if models['ll'][u]!=(0,0,0,0)}
        check(Zcross==set(rec['transverse_lambda_states']),'MIXED_TRANSVERSE_LIST')
        # Verify the partition of internal zero edges, independently of cycle discovery.
        zeroedges={e for e in range(4*n) if point['sl'][e//4]==point['sl'][e%n]==0}
        seen=set()
        for cr in rec['cycle_records']:
            cy=cr['edges'];check(bool(cy) and len(cy)<=k,'MIXED_CYCLE_LENGTH')
            check(not (seen&set(cy)),'MIXED_CYCLE_DUPLICATE');seen.update(cy)
            check(all(models['rot'][e]==cy[(i+1)%len(cy)] for i,e in enumerate(cy)),'MIXED_ROTATION')
            check(len({e//4 for e in cy})==len(cy),'MIXED_CYCLE_VERTEX_REPEAT')
            for t in range(3):
                selected=[e for e in cy if ts[j+t][e]]
                check(selected==[cr['representatives'][t]],'MIXED_REPRESENTATIVE')
            broken=any(e//4 in Zcross for e in cy)
            check(broken==cr['transverse'] and (not broken or len(cy)==k),'MIXED_PRIMITIVITY')
        check(seen==zeroedges,'MIXED_INTERNAL_COVERAGE')
        check(rec['transverse_null_pcr_count']==sum(x['transverse'] for x in rec['cycle_records']),'MIXED_TRANSVERSE_PCR_COUNT')
        for comp,(s,t) in zip(rec['comparisons'],((0,1),(1,2),(0,2))):
            verify_sparse_certificate(ts[j+s],ts[j+t],comp)
            check(comp['source']==A['rows'][j+s]['id'] and comp['target']==A['rows'][j+t]['id'],'MIXED_PAIR_ID')
            delta={e:ts[j+t][e]-ts[j+s][e] for e in range(4*n) if ts[j+t][e]!=ts[j+s][e]}
            stored=dict(comp['patch']);check(len(stored)==len(comp['patch']),'PATCH_DUPLICATES');check(stored==delta,'PATCH_LITERAL_IDENTITY')
            da=[cs[j+t]['a'][u]-cs[j+s]['a'][u] for u in range(n)]
            dh=[int(cs[j+t]['sl'][u]>0)-int(cs[j+s]['sl'][u]>0) for u in range(n)]
            regular=sum(da[u]!=0 for u in S);external=[]
            for u in sorted(Z):
                if not da[u] and not dh[u]:continue
                for ch in range(4):
                    e=4*u+ch
                    if e in zeroedges:continue
                    x=mod.at(models['re'][e],rr['value'])
                    check(x!=0,'EXTERNAL_REAL_ZERO')
                    c=da[u]-dh[u]*int(x<0)
                    if c:external.append([u,ch,c])
            check(external==comp['external_exchanges'],'EXTERNAL_COEFFICIENT_CERT')
            changed=[i for i,c in enumerate(rec['cycle_records']) if c['representatives'][s]!=c['representatives'][t]]
            check(changed==comp['changed_cycle_indices'],'CHANGED_CYCLE_LIST')
            check(regular==comp['regular_midpoint_states'],'REGULAR_COUNT')
            check(len(external)==comp['external_companion_pairs'],'EXTERNAL_COUNT')
            check(len(changed)==comp['changed_null_pcrs'],'CYCLE_COUNT')
            check(len(delta)==comp['symmetric_difference']==8*regular+2*len(external)+2*len(changed),'MIXED_COUNT_IDENTITY')
            checked+=1;words+=4*n
        # Check frozen-midpoint residuals with a fresh full DAG computation.
        for t,fr in enumerate(rec['frozen_midpoint_residuals']):
            d=[0]*n
            for u in S:d[u]=cs[j+1]['a'][u]-cs[j+t]['a'][u]
            b=bytearray(4*n)
            for e in range(4*n):
                x=ts[j+t][e]+d[e//4]-d[e%n];check(x in (0,1),'VERIFY_FROZEN_BOOLEAN');b[e]=x
            check(sha(pack(b))==fr['bits_sha256'],'FROZEN_BITS')
            check(all(b[e]==b[models['rcc'][e]] for e in range(4*n)),'FROZEN_RC')
            check(sum(b)==A['one_per_pcr'],'FROZEN_SIZE')
            _,L,_,_=dag(b,n);check(L==fr['length'],'FROZEN_LONGEST_PATH')
            check(abs(L-A['rows'][j+t]['lower_longest_edges'])<=2,'FROZEN_VARIATION')
            pp=fr['longest_path'];check(len(pp)==L+1,'FROZEN_PATH_SIZE')
            check(all(u% (n//4)==v//4 and not b[u*4+v%4] for u,v in zip(pp,pp[1:])),'FROZEN_PATH_EDGES')
    return checked,words

def verify_binary(out):
    r=json.loads((out/'binary_parity_obstruction.json').read_text());q=2;k=4;n=8
    exact=[]
    # Independent carrier: all 65536 masks, not one-per-PCR products.
    for mask in range(1<<16):
        if mask.bit_count()!=6 or not (mask&1 and mask&(1<<15)):continue
        m=bytearray((mask>>e)&1 for e in range(16))
        if any(m[e]!=m[rc(e,4,2)] for e in range(16)):continue
        try:dag(m,n,q)
        except CheckError:continue
        exact.append(tuple(e for e,v in enumerate(m) if v))
    check(set(exact)==set(map(tuple,r['selected_edge_codes'])) and len(exact)==r['rc_mds_count']==6,'BINARY_INDEPENDENT_CENSUS')
    ms=[bytearray(int(e in row) for e in range(16)) for row in r['selected_edge_codes']]
    for c in r['odd_parity_pairs']:
        m1,m2=ms[c['source']],ms[c['target']];f=c['potential'];C=c['constant']
        check(len(f)==8 and all(m2[e]-m1[e]==f[e//2]-f[e%8] for e in range(16)),'BINARY_CUT_CERTIFICATE')
        check(C%2==1 and all(f[u]+f[rc(u,3,2)]==C for u in range(8)),'BINARY_ODD_PARITY')
    for c in r['negative_cycles']:verify_separator(ms[c['source']],ms[c['target']],k,c,q)
    # Direct paired move graph: there are no admissible nontrivial steps.
    count=0
    for m in ms:
        for u in range(8):
            v=rc(u,3,2);iu=set(inedges(u,n,q));ou=set(outedges(u,n,q));iv=set(inedges(v,n,q));ov=set(outedges(v,n,q))
            if (iu|ou)&(iv|ov):continue
            if all(m[e] for e in iu) and all(m[e] for e in ov):count+=1
    check(count==0,'BINARY_PAIRED_EDGES')
    return {'all_masks':65536,'rc_mds':6,'ordinary_F_components_met':3,'paired_edges':0,'paired_components':6}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path);ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results');ap.add_argument('--parent-replay',action='store_true');args=ap.parse_args()
    parent,tmp=locate_parent(args.parent);start=time.time();receipt=[]
    spec=importlib.util.spec_from_file_location('inherited_independent_verifier',parent/'code/verify_atlas.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    try:
        for k in (4,6,8):
            A,ts=parent_data(parent,k);models=mod.get_models(k);mod.verify_catalog(A,models);cs=[]
            for row,t in zip(A['rows'],ts):
                c=mod.compute(models,row['ratio']);check(c['term']==t,'INDEPENDENT_PARENT_TERM');cs.append(c)
                if args.parent_replay:
                    mod.verify_certificate(row,models,c,(parent/f"results/certificates/{row['id']}.bits").read_bytes(),(parent/f"results/certificates/{row['id']}.levels").read_bytes())
            if args.parent_replay:mod.verify_global(A,models,cs)
            R=json.loads((args.out/f'components_k{k}.json').read_text());lookup={r['id']:i for i,r in enumerate(A['rows'])}
            flat=[rid for cls in R['classes'] for rid in cls]
            check(len(flat)==len(set(flat))==len(lookup) and set(flat)==set(lookup),'COMPONENT_COVERAGE')
            check(len(R['classes'])==R['component_count'],'COMPONENT_COUNT')
            for cls in R['classes']:
                check(all(R['row_component'][r]==R['row_component'][cls[0]] for r in cls),'COMPONENT_LABEL')
                for r in cls:
                    rec=next((x for x in R['positive_geodesics'] if x['source']==cls[0] and x['target']==r),None)
                    check(rec is not None,'MISSING_POSITIVE_GEODESIC')
                    replay_pair(ts[lookup[cls[0]]],ts[lookup[r]],k,rec)
            lambda_roots=[r['value'] for r in A['roots'] if 'lower_lambda' in r['families']]
            for row in A['rows']:
                expected=2*sum(mod.cmp(r,row['ratio'])<0 for r in lambda_roots)+int(any(mod.cmp(r,row['ratio'])==0 for r in lambda_roots))
                check(R['row_component'][row['id']]==expected,'LAMBDA_COMPONENT_DIAGRAM')
            negpairs=set()
            for rec in R['negative_cycles']:
                i,j=rec['source_class'],rec['target_class'];check(i<j,'NEGATIVE_ORDER')
                check(rec['source']==R['classes'][i][0] and rec['target']==R['classes'][j][0],'NEGATIVE_REP_ID')
                check((i,j) not in negpairs,'NEGATIVE_DUPLICATE');negpairs.add((i,j))
                verify_separator(ts[lookup[rec['source']]],ts[lookup[rec['target']]],k,rec)
            check(negpairs==set(combinations(range(len(R['classes'])),2)),'NEGATIVE_PAIR_COVERAGE')
            wc,words=verify_mixed(parent,args.out,k,A,ts,mod,models,cs)
            row={'k':k,'parent_rows':len(ts),'parent_word_decisions':len(ts)*4**k,'new_components':len(R['classes']),
                 'positive_geodesics':len(R['positive_geodesics']),'paired_moves':sum(len(x['moves']) for x in R['positive_geodesics']),
                 'negative_cycle_pairs':len(negpairs),'mixed_pairs':wc,'mixed_word_decisions':words}
            receipt.append(row);print('VERIFIED',row,flush=True)
        binary=verify_binary(args.out)
        result={'status':'PASS','binary_control':binary,'parent_full_certificate_replay':args.parent_replay,'orders':receipt,
                'independence':'Original-seed reconstruction and Fraction-pair signs imported from inherited verifier, not new wall producer; positive paths checked by set operations, negatives by explicit cycle hits. No wall-discovery import.'}
        dump(args.out/'verification_round3.json',result)
        print('VERIFY_PASS seconds',round(time.time()-start,3),flush=True)
    finally:
        if tmp:tmp.cleanup()
if __name__=='__main__':main()
