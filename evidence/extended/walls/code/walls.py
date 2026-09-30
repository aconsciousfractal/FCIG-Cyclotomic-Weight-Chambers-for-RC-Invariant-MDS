"""Exact mixed-wall decomposition, including the central stratum.
The baseline terms are inherited evidence. Only local one-jet forms,
external companion exchanges, and null PCR representatives enter the update.
"""
from core import *
from collections import Counter
from array import array

def wall_record(k,rec,rows,terms,ff):
    n=4**(k-1);j=2*rec['id'];r=rec['value'];lf,rf,cf=ff
    idx=[j,j+1,j+2];tr=[terms[i] for i in idx]
    state3=[state(k,rows[i]['ratio'],ff) for i in idx]
    sl0,sr0,a0=state3[1]
    Z={u for u in range(n) if sl0[u]==0}
    S={u for u in range(n) if sl0[u]>0 and sr0[u]==0 and der(rf[u])!=0}
    cross={u for u in Z if der(lf[u])!=0}
    # One-jet signs, not numerical epsilon.
    for t,st in zip((-1,0,1),state3):
        sl,sr,a=st
        for u in Z:check(sl[u]==t*der(lf[u]),'LAMBDA_GERM')
        for u in range(n):
            if sr0[u]==0:check(sr[u]==t*der(rf[u]),'MIDPOINT_GERM')
    cycles=singular_cycles(k,Z)
    edges0={e for c in cycles for e in c};reps=[];cycle_records=[]
    for cy in cycles:
        rr=[]
        for ii,(sl,sr,a) in enumerate(state3):
            on=[]
            for e in cy:
                u,v=e//4,e%n
                seed=tr[1][e] if sl[u]==sl[v]==0 else int(sl[u]<=0<sl[v])
                m=seed+a[u]-a[v]
                check(m in (0,1),'LOCAL_CYCLE_NONBOOLEAN')
                if m:on.append(e)
                check(m==tr[ii][e],'LOCAL_CYCLE_REPLAY')
            check(len(on)==1,'CYCLE_SINGLE_SELECTION');rr.append(on[0])
        broken=any(e//4 in cross for e in cy)
        if broken:check(len(cy)==k,'TRANSVERSE_NONPRIMITIVE')
        rvals={value(rf[e//4],r) for e in cy}
        check(len(rvals)==1,'NULL_PCR_POTENTIAL_NOT_CONSTANT')
        cycle_records.append({'edges':cy,'representatives':rr,'transverse':broken,'twice_midpoint_numerator':list(next(iter(rvals)))})
        reps.append(rr)
    isolated=sorted(Z-{e//4 for e in edges0})
    external={}
    for z in sorted(cross|{u for u in Z if any(st[2][u] for st in state3)}):
        a,b,c,d=cf[z];le=[]
        for ch,(h,l) in enumerate(V):
            ex=z*4+ch;en=ch*n+z;x=ev((a+2*h,b+2*l,c,d),r)
            if x==0:
                check(ex in edges0 and en in edges0,'ZERO_COMPANION_NOT_INTERNAL');continue
            check(ex not in edges0 and en not in edges0,'EXTERNAL_COMPANION_INTERNAL')
            check(sl0[ex%n]==-x and sl0[en//4]==x,'COMPANION_OPPOSITE_SIGNS')
            le.append((ch,x))
        external[z]=le
    comparisons=[]
    for s,t in ((0,1),(1,2),(0,2)):
        da=[state3[t][2][u]-state3[s][2][u] for u in range(n)]
        dh=[int(state3[t][0][u]>0)-int(state3[s][0][u]>0) for u in range(n)]
        patch={};support_regular=set();middle=0;ext=0;cycnt=0;extrecords=[]
        def add(e,v):
            if v:patch[e]=patch.get(e,0)+v
        for u in sorted(S):
            if da[u]==0:continue
            middle+=1
            for e in outedges(u,n):add(e,da[u]);support_regular.add(e)
            for e in inedges(u,n):add(e,-da[u]);support_regular.add(e)
        check(len(support_regular)==8*middle,'REGULAR_STARS_OVERLAP')
        for z,vals in external.items():
            for ch,x in vals:
                c=da[z]-dh[z]*int(x<0)
                check(c in (-1,0,1),'EXTERNAL_COEFFICIENT')
                if c==0:continue
                e=z*4+ch;f=ch*n+z
                check(e not in support_regular and f not in support_regular,'REGULAR_SINGULAR_INTERFERENCE')
                check(e not in patch and f not in patch,'EXTERNAL_STARS_OVERLAP')
                add(e,c);add(f,-c);ext+=abs(c);extrecords.append([z,ch,c])
        changed_cycles=[]
        for ci,rr in enumerate(reps):
            if rr[s]!=rr[t]:
                check(rr[s] not in patch and rr[t] not in patch,'NULL_CYCLE_SUPPORT_OVERLAP')
                add(rr[s],-1);add(rr[t],1);cycnt+=1;changed_cycles.append(ci)
        patch={e:x for e,x in patch.items() if x}
        check(all(x in (-1,1) for x in patch.values()),'PATCH_NONBOOLEAN_DELTA')
        for e in range(4*n):check(patch.get(e,0)==tr[t][e]-tr[s][e],'SPARSE_WALL_IDENTITY')
        count=8*middle+2*ext+2*cycnt
        check(count==len(patch),'EXACT_WALL_COUNT')
        comps={'source':rows[idx[s]]['id'],'target':rows[idx[t]]['id'],
               'regular_midpoint_states':middle,'external_companion_pairs':ext,
               'changed_null_pcrs':cycnt,'symmetric_difference':count,
               'external_exchanges':extrecords,'changed_cycle_indices':changed_cycles,
               'patch':[[e,patch[e]] for e in sorted(patch)]}
        comparisons.append(comps)
    # The deletion of the exceptional states leaves one literal common induced DAG.
    C=[u not in Z and u not in S for u in range(n)]
    for e in range(4*n):
        if C[e//4] and C[e%n]:check(tr[0][e]==tr[1][e]==tr[2][e],'COMMON_CORE_MISMATCH')
    _,LC,pc,dc=dag(tr[1],n,allowed=C)
    paths=[]
    for ii,row in enumerate(rows[x] for x in idx):
        pp=row['longest_path_lower_vertex_codes'];runs=[];a=0
        while a<len(pp):
            if pp[a] not in Z:a+=1;continue
            b=a+1
            while b<len(pp) and pp[b] in Z:b+=1
            rr=pp[a:b];check(len(rr)<=k,'SINGULAR_PATH_EXCEEDS_ORDER')
            runs.append({'start_index':a,'end_index':b-1,'vertices':rr});a=b
        paths.append({'row':row['id'],'length':row['lower_longest_edges'],
                      'singular_vertices':sum(len(x['vertices']) for x in runs),'singular_runs':runs})
    frozen=[]
    if cross:
        for ii in range(3):
            d=[0]*n
            for u in S:d[u]=state3[1][2][u]-state3[ii][2][u]
            b=bytearray(4*n)
            for e in range(4*n):
                v=tr[ii][e]+d[e//4]-d[e%n]
                check(v in (0,1),'FROZEN_MIDPOINT_NONBOOLEAN');b[e]=v
            check(sum(b)==sum(tr[ii]),'FROZEN_MIDPOINT_CARDINALITY')
            check(all(b[e]==b[rc(e,k)] for e in range(4*n)),'FROZEN_MIDPOINT_RC')
            _,L,p,_=dag(b,n)
            check(abs(L-rows[idx[ii]]['lower_longest_edges'])<=2,'FROZEN_LENGTH_BOUND')
            frozen.append({'length':L,'longest_path':p,'bits_sha256':sha(pack(b))})
    return {'k':k,'root_id':rec['id'],'ratio':r,'lambda_wall':bool(cross),
            'singular_vertices':len(Z),'transverse_lambda_states':sorted(cross),
            'isolated_singular_count':len(isolated),'regular_midpoint_states':sorted(S),
            'singular_cycle_profile':dict(sorted(Counter(len(c) for c in cycles).items())),
            'transverse_null_pcr_count':sum(c['transverse'] for c in cycle_records),
            'cycle_records':cycle_records,'comparisons':comparisons,
            'lengths':[rows[x]['lower_longest_edges'] for x in idx],
            'common_induced_core_length':LC,'common_core_path':pc,'longest_path_decompositions':paths,'frozen_midpoint_residuals':frozen}

def run(parent,out):
    result=[];summary=[]
    for k in (4,6,8):
        a,terms=parent_data(parent,k);ff=forms(k)
        for rr in a['roots']:
            r=wall_record(k,rr,a['rows'],terms,ff)
            if r['lambda_wall']:
                dump(out/f'mixed_k{k}_r{rr["id"]:03d}.json',r)
                summary.append({key:r[key] for key in ('k','root_id','ratio','singular_vertices','isolated_singular_count','transverse_null_pcr_count','lengths','common_induced_core_length') }|
                    {'transverse_lambda_count':len(r['transverse_lambda_states']),
                     'regular_midpoint_count':len(r['regular_midpoint_states']),
                     'jump_parts':[{key:c[key] for key in ('regular_midpoint_states','external_companion_pairs','changed_null_pcrs','symmetric_difference')} for c in r['comparisons']]})
            result.append({'k':k,'root_id':rr['id'],'mixed':r['lambda_wall'],'comparison_count':3,
                           'symmetric_differences':[x['symmetric_difference'] for x in r['comparisons']]})
        print('WALLS_OK',k,len(a['roots']),flush=True)
    dump(out/'wall_summary.json',{'status':'PASS','walls':len(result),'comparisons':3*len(result),'mixed_walls':summary,'all_walls':result})
    return summary
