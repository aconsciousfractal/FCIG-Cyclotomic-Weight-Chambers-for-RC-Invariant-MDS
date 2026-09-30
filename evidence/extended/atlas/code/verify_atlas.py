"""Independent exact verifier for the full k=4,6,8 DNA atlas.
No imports from the producer or the inherited sign oracle.
Uses Fraction pairs and Q[X]/Phi_k for original seed reconstruction,
compressed symmetric-pair enumeration for walls, and sink elimination.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
from collections import Counter,deque
from functools import lru_cache,cmp_to_key
from pathlib import Path
import json,hashlib,math,struct,sys,time,copy

class VerificationError(RuntimeError):pass
def chk(b,code):
    if not b:raise VerificationError(code)
def sh(data):return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def bh(data):return hashlib.sha256(data).hexdigest()
def plus(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def times(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def invert(x):
    d=x[0]**2-2*x[1]**2;chk(d!=0,'DIVISION_ZERO');return x[0]/d,-x[1]/d
def scalar(n):return F(n),F(0)
def qvalue(t):return F(t[0],t[2]),F(t[1],t[2])
def triple(x):
    d=math.lcm(x[0].denominator,x[1].denominator);a=int(x[0]*d);b=int(x[1]*d);g=math.gcd(math.gcd(abs(a),abs(b)),d)
    return a//g,b//g,d//g
def sg(x):
    a,b=x
    if b==0:return (a>0)-(a<0)
    if a==0:return (b>0)-(b<0)
    if a>0 and b>0:return 1
    if a<0 and b<0:return -1
    dif=a*a-2*b*b
    return ((dif>0)-(dif<0)) * ((a>0)-(a<0))
def cmp(t,u):return sg(plus(qvalue(t),neg(qvalue(u))))
def at(f,r):
    a,b,c,d=f
    return sg(plus(times((F(a),F(c)),qvalue(r)),(F(b),F(d))))
def zero_ratio(f):
    a,b,c,d=f
    if a==c==0:return None
    return triple(times(neg((F(b),F(d))),invert((F(a),F(c)))))
def vadd(x,y):return x[0]+y[0],x[1]+y[1]
V=((-1,0),(0,-1),(0,1),(1,0))
D=Counter(vadd(x,y) for x in V for y in V)

def compressed_roots(k):
    walls=set();lr=set()
    for x in product(D,repeat=(k-2)//2):
        for z in V:
            if k==4:
                h=-2*x[0][0]-4*z[0];l=-2*x[0][1]-4*z[1];f=(h,l,0,0)
            elif k==6:
                h=-x[0][0]-3*x[1][0]-4*z[0];l=-x[0][1]-3*x[1][1]-4*z[1];f=(h,l,0,0)
            else:
                h=-2*sum(y[0] for y in x)-4*z[0];l=-2*sum(y[1] for y in x)-4*z[1]
                f=(h,l,x[0][0]-x[2][0],x[0][1]-x[2][1])
            r=zero_ratio(f)
            if r and cmp(r,(1,0,1))>0:walls.add(r)
        if k==4:f=(x[0][0],x[0][1],0,0)
        elif k==6:f=(x[0][0]+x[1][0],x[0][1]+x[1][1],0,0)
        else:f=(x[0][0]+x[2][0],x[0][1]+x[2][1],x[1][0],x[1][1])
        r=zero_ratio(f)
        if r and cmp(r,(1,0,1))>0:lr.add(r)
    return sorted(walls,key=cmp_to_key(cmp)),lr

def embedding8():
    # Four antipodal differences with their exact ordered-letter multiplicities.
    counts=Counter();structural=0
    for ds in product(D,repeat=4):
        mult=math.prod(D[d] for d in ds)
        if all(d==(0,0) for d in ds):structural+=mult;continue
        h=[d[0] for d in ds];l=[d[1] for d in ds]
        re=(2*h[0],2*l[0],h[1]-h[3],l[1]-l[3]);im=(2*h[2],2*l[2],h[1]+h[3],l[1]+l[3])
        r=zero_ratio(re) or zero_ratio(im)
        if r and cmp(r,(1,0,1))>0 and at(re,r)==at(im,r)==0:counts[r]+=mult
    return structural,counts

def encode(w):
    n=0
    for ch in w:n=4*n+ch
    return n
def rc(w):return tuple(3-ch for ch in w[::-1])
@lru_cache(maxsize=None)
def bzero(w):
    own={w[i:]+w[:i] for i in range(len(w))};mate=rc(w)
    if mate in own:
        fixed=[v for v in own if v==rc(v)];chk(len(fixed)==2,'TIE_FIXED_COUNT')
        return int(w==min(fixed))
    other={mate[i:]+mate[:i] for i in range(len(w))};small=min(own|other)
    return int(w==small or w==rc(small))

def powers(k):
    phi={4:(1,0,1),6:(1,-1,1),8:(1,0,0,0,1)}[k];d=len(phi)-1
    p=[1]+[0]*(d-1);ans=[]
    for _ in range(k):
        z=[0]+p;c=z.pop()
        p=[z[i]-c*phi[i] for i in range(d)];ans.append(tuple(p))
    chk(ans[-1]==tuple([1]+[0]*(d-1)),'ROOT_POWER');return ans

def scalar_reim(k,x):
    if k==4:return (2*x[0],0),(x[1],0)
    if k==6:return (2*x[0]+x[1],0),(x[1],0)
    return (2*x[0],x[1]-x[3]),(2*x[2],x[1]+x[3])

def make_forms(k,n):
    pp=powers(k);d=len(pp[0]);re=[];im=[];rm=[];ww=list(product(range(4),repeat=n))
    for w in ww:
        h=[0]*d;l=[0]*d;hs=ls=0
        for i,ch in enumerate(w):
            a,b=V[ch];hs+=a;ls+=b
            for j in range(d):h[j]+=a*pp[i][j];l[j]+=b*pp[i][j]
        hr,hi=scalar_reim(k,h);lr,li=scalar_reim(k,l)
        re.append((hr[0],lr[0],hr[1],lr[1]));im.append((hi[0],li[0],hi[1],li[1]))
        rm.append((hr[0]-2*hs,lr[0]-2*ls,hr[1],lr[1]))
    return ww,re,im,rm

def get_models(k):
    ww,re,im,_=make_forms(k,k);lw,_,ll,mr=make_forms(k,k-1)
    rcc=[encode(rc(w)) for w in ww];lrc=[encode(rc(w)) for w in lw]
    rot=[encode(w[1:]+w[:1]) for w in ww];orbits={}
    for i,w in enumerate(ww):
        if i not in orbits:
            cur=i;cyc=[]
            while not cyc or cur!=i:cyc.append(cur);cur=rot[cur]
            for c in cyc:orbits[c]=min(cyc)
    canon=[orbits[i] for i in range(len(ww))]
    return {'k':k,'ww':ww,'lw':lw,'re':re,'im':im,'ll':ll,'mr':mr,'rcc':rcc,'lrc':lrc,'rot':rot,'canon':canon,'pcrcount':len(set(canon))}

def compute(m,r):
    forms=set(m['re'])|set(m['im'])|set(m['ll'])|set(m['mr']);ev={f:at(f,r) for f in forms}
    sr=[ev[f] for f in m['re']];si=[ev[f] for f in m['im']];sl=[ev[f] for f in m['ll']];sm=[ev[f] for f in m['mr']]
    a=[int(sl[i]>0 and (sm[i]<0 or (sm[i]==0 and i<m['lrc'][i]))) for i in range(len(sl))]
    n=len(sl);term=[];nzero=0
    for i,w in enumerate(m['ww']):
        if sr[i]==si[i]==0:b=bzero(w);nzero+=1
        elif si[i]==0:b=int(sr[i]<0)
        else:b=int(si[i]<0 and si[m['rot'][i]]>0)
        val=b+a[encode(w[:-1])]-a[encode(w[1:])]
        chk(val in (0,1),'ORIGINAL_RULE_NONBOOLEAN');term.append(val)
    return {'term':bytes(term),'a':a,'sl':sl,'sm':sm,'zero':nzero,
      'col':[0 if sl[i]<=0 else (1 if a[i] else -1) for i in range(n)]}

def verify_certificate(row,m,c,bits,levelbytes):
    term=c['term'];n=len(c['a']);size=len(term)
    chk(len(bits)==(size+7)//8,'BITSET_SIZE')
    actual=bytes((bits[i//8]>>(i%8))&1 for i in range(size))
    chk(actual==term,'SELECTOR_MISMATCH')
    chk(bh(bits)==row['bitset_sha256'],'BITSET_DIGEST')
    chk(sh([i for i,v in enumerate(term) if v])==row['terminal_sha256'],'TERMINAL_DIGEST')
    chk(sum(term)==row['selected_count']==m['pcrcount'],'CARDINALITY')
    counts=Counter(m['canon'][i] for i,v in enumerate(term) if v)
    chk(len(counts)==m['pcrcount'] and all(v==1 for v in counts.values()),'ONE_PER_PCR')
    chk(all(term[i]==term[m['rcc'][i]] for i in range(size)),'LITERAL_RC')
    chk(row['zero_embedding_words']==c['zero'],'ZERO_EMBEDDING_COUNT')
    chk(row['zero_lower_words']==c['sl'].count(0),'ZERO_LOWER_COUNT')
    chk(row['positive_lower_words']==sum(x>0 for x in c['sl']),'POSITIVE_LOWER_COUNT')
    chk(row['ideal_count']==sum(c['a']),'IDEAL_COUNT')
    chk(row['midpoint_active_ties']==sum(c['sl'][i]>0 and c['sm'][i]==0 for i in range(n)),'ACTIVE_TIE_COUNT')
    chk(row['midpoint_active_accidental_ties']==sum(c['sl'][i]>0 and c['sm'][i]==0 and m['mr'][i]!=(0,0,0,0) for i in range(n)),'ACCIDENTAL_TIE_COUNT')
    chk(len(levelbytes)==2*n,'LEVEL_SIZE');levels=struct.unpack('<'+'H'*n,levelbytes)
    chk(bh(levelbytes)==row['level_sha256'],'LEVEL_DIGEST')
    L=row['lower_longest_edges'];chk(min(levels)>=0 and max(levels)==L,'LEVEL_BOUNDS')
    pred=[[] for _ in range(n)];outdeg=[0]*n;ecount=0
    for i,v in enumerate(term):
        if v:continue
        # Lower-state endpoints read from tuple words, independent of the producer formula.
        u=encode(m['ww'][i][:-1]);v=encode(m['ww'][i][1:]);ecount+=1
        chk(levels[v]>=levels[u]+1,'LEVEL_EDGE_INEQUALITY')
        chk(c['col'][u]<=c['col'][v],'BLOCK_ORDER')
        pred[v].append(u);outdeg[u]+=1
    chk(ecount==row['residual_arcs'],'RESIDUAL_ARC_COUNT')
    path=row['longest_path_lower_vertex_codes']
    chk(len(path)==L+1 and len(set(path))==len(path),'PATH_SIZE_OR_REPEAT')
    chk(all(0<=v<n for v in path),'PATH_VERTEX')
    for u,v in zip(path,path[1:]):
        x=m['lw'][u];y=m['lw'][v]
        chk(x[1:]==y[:-1],'PATH_OVERLAP')
        chk(term[encode(x+(y[-1],))]==0,'PATH_DELETED_EDGE')
    # Second longest-path computation works backwards from sinks, no saved labels.
    q=deque(i for i in range(n) if outdeg[i]==0);dp=[0]*n;local=[0]*n;seen=0
    while q:
        v=q.popleft();seen+=1
        for u in pred[v]:
            dp[u]=max(dp[u],dp[v]+1)
            if c['col'][u]==c['col'][v]:local[u]=max(local[u],local[v]+1)
            outdeg[u]-=1
            if outdeg[u]==0:q.append(u)
    chk(seen==n and max(dp)==L,'SINK_DP')
    block={str(col):max((local[i] for i in range(n) if c['col'][i]==col),default=0) for col in (-1,0,1)}
    chk(block==row['block_lengths'],'BLOCK_LENGTHS')
    chk(row['block_upper_bound']==2*block['1']+block['0']+2,'BLOCK_BOUND')
    chk(row['block_bound_saturated']==(L==row['block_upper_bound']),'BLOCK_SATURATION')
    chk(row['upper_longest_edges']==L-1,'LINE_GRAPH_UNITS')

def verify_catalog(A,m):
    k=A['k'];roots,lam=compressed_roots(k)
    chk([tuple(r['value']) for r in A['roots']]==roots,'ROOT_COMPLETENESS_ORDER')
    chk(len(A['rows'])==2*len(roots)+1,'CELL_COVERAGE')
    chk(A['open_intervals']==len(roots)+1 and A['point_strata']==len(roots),'STRATA_COUNTS')
    for i,r in enumerate(A['roots']):
        chk(r['id']==i,'ROOT_ID')
        flags=r['families'];chk(('lower_lambda' in flags)==(roots[i] in lam),'LAMBDA_ROOT_TAG')
        chk('midpoint' in flags,'MIDPOINT_ROOT_TAG')
        u=r['active_midpoint_word_code'];chk(u is not None and 0<=u<len(m['lw']),'ACTIVE_WITNESS_RANGE')
        chk(m['mr'][u]!=(0,0,0,0) and at(m['mr'][u],roots[i])==0 and at(m['ll'][u],roots[i])>0,'ACTIVE_WITNESS')
        chk(r['active_midpoint_word']==''.join('ACGT'[c] for c in m['lw'][u]),'ACTIVE_WITNESS_WORD')
    # Also check every raw seed root is among the retained roots.
    fs={'seed_real':m['re'],'seed_imag':m['im'],'lower_lambda':m['ll'],'midpoint':m['mr']}
    for name,forms in fs.items():
        rs={zero_ratio(f) for f in set(forms)};rs={r for r in rs if r and cmp(r,(1,0,1))>0}
        chk(rs<=set(roots),'MISSING_SEED_ROOT')
        chk(A['family_counts'][name]['positive_roots']==len(rs),'FAMILY_ROOT_COUNT')
        chk(A['family_counts'][name]['structural_occurrences']==forms.count((0,0,0,0)),'FAMILY_STRUCTURAL_COUNT')
        chk(A['family_counts'][name]['distinct_forms']==len(set(forms)),'FAMILY_DISTINCT_COUNT')
        chk(A['family_counts'][name]['word_occurrences']==len(forms),'FAMILY_OCCURRENCES')
        for r in A['roots']:chk((name in r['families'])==(tuple(r['value']) in rs),'FAMILY_TAG')
    for j,row in enumerate(A['rows']):
        chk(row['id']==f'k{k}_c{j:03d}','CELL_ID');x=row['ratio']
        if j%2:
            chk(row['kind']=='point' and tuple(x)==roots[j//2],'POINT_VALUE')
        else:
            lo=roots[j//2-1] if j else (1,0,1);hi=roots[j//2] if j//2<len(roots) else None
            chk(row['kind']=='interval' and x[1]==0 and cmp(lo,x)<0 and (hi is None or cmp(x,hi)<0),'INTERVAL_SAMPLE')
    if k==8:
        structural,counts=embedding8();chk(structural==A['structural_embedding_words'],'STRUCTURAL_EMBEDDING')
        chk({tuple(r['value']):r['accidental_words'] for r in A['embedding_spectrum']}==dict(counts),'EMBEDDING_SPECTRUM')
    for rec in A['embedding_spectrum']:
        i=rec['witness_code'];r=rec['value']
        chk(at(m['re'][i],r)==at(m['im'][i],r)==0 and (m['re'][i]!=(0,0,0,0) or m['im'][i]!=(0,0,0,0)),'EMBEDDING_WITNESS')
        chk(rec['witness_word']==''.join('ACGT'[c] for c in m['ww'][i]),'EMBEDDING_WITNESS_WORD')
    return roots,lam

def verify_global(A,models,computed):
    terms=[c['term'] for c in computed];classes={}
    for row,t in zip(A['rows'],terms):classes.setdefault(t,[]).append(row['id'])
    chk(A['terminal_classes']==list(classes.values()) and A['distinct_terminals']==len(classes),'TERMINAL_CLASSES')
    essential=0;wallstats=[]
    for i,rec in enumerate(A['roots']):
        tri=terms[2*i:2*i+3];trirows=A['rows'][2*i:2*i+3];tricomputed=computed[2*i:2*i+3]
        for nm,x,y in [('left_to_right',0,2),('left_to_point',0,1),('point_to_right',1,2)]:
            ds=[j for j in range(len(tri[x])) if tri[x][j]!=tri[y][j]];r=rec[nm]
            chk(r['symmetric_difference']==len(ds),'WALL_DIFFERENCE_COUNT')
            chk(r['first_witness_code']==(ds[0] if ds else None),'WALL_DIFFERENCE_WITNESS')
            chk(r['first_witness_word']==(''.join('ACGT'[c] for c in models['ww'][ds[0]]) if ds else None),'WALL_DIFFERENCE_WORD')
        e=len(set(tri))>1;essential+=e
        chk(rec['essential']==e,'ESSENTIAL_WALL')
        chk(rec['point_only_anomaly']==(tri[0]==tri[2]!=tri[1]),'POINT_ANOMALY')
        chk(rec['three_distinct']==(len(set(tri))==3),'THIRD_TERMINAL')
        if 'lower_lambda' not in rec['families']:
            active=[j for j in range(len(models['lw'])) if tricomputed[1]['sl'][j]>0 and tricomputed[1]['sm'][j]==0 and models['mr'][j]!=(0,0,0,0)]
            S=set(active)
            chk(all(encode(models['lw'][u][1:]+(ch,)) not in S for u in S for ch in range(4)),'ACTIVE_STATES_NOT_INDEPENDENT')
            chk(rec['left_to_right']['symmetric_difference']==8*len(S),'WALL_EXACT_JUMP')
            aa=[c['a'] for c in tricomputed]
            for nm,x,y in [('left_to_point',0,1),('point_to_right',1,2)]:
                changed=sum(u!=v for u,v in zip(aa[x],aa[y]))
                chk(rec[nm]['symmetric_difference']==8*changed,'WALL_PARTIAL_JUMP')
            L=[r['lower_longest_edges'] for r in trirows];chk(max(L)-min(L)<=2,'WALL_LENGTH_LIPSCHITZ')
            wallstats.append({'root_id':i,'active_states':len(S),'symmetric_difference':rec['left_to_right']['symmetric_difference'],'length_span':max(L)-min(L)})
    chk(A['essential_roots']==essential,'ESSENTIAL_COUNT')
    values=[r['lower_longest_edges'] for r in A['rows']];mn=min(values);mx=max(values)
    chk(A['optimum_lower_longest_edges']==mn and A['maximum_lower_longest_edges']==mx,'OPTIMUM_VALUE')
    chk(A['optimum_cells']==[r['id'] for r in A['rows'] if r['lower_longest_edges']==mn],'OPTIMUM_LOCUS')
    chk(A['length_histogram']=={str(k):v for k,v in sorted(Counter(values).items())},'LENGTH_HISTOGRAM')
    chk(A['block_bound_nonsaturated_cells']==[r['id'] for r in A['rows'] if not r['block_bound_saturated']],'SATURATION_LOCUS')
    return wallstats

def main():
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'results';t=time.time();receipt=[]
    for k in (4,6,8):
        A=json.loads((root/f'atlas_k{k}.json').read_text());m=get_models(k);verify_catalog(A,m);computed=[]
        for row in A['rows']:
            c=compute(m,row['ratio']);cid=row['id'];verify_certificate(row,m,c,(root/'certificates'/f'{cid}.bits').read_bytes(),(root/'certificates'/f'{cid}.levels').read_bytes());computed.append(c)
        walls=verify_global(A,m,computed)
        receipt.append({'k':k,'rows':len(computed),'word_decisions':len(computed)*4**k,'critical_roots':len(A['roots']),'distinct_terminals':A['distinct_terminals'],'regular_wall_checks':walls,'status':'PASS'})
        print('VERIFIED',k,'cells',len(computed),'all original-seed decisions, path certificates and exact optimum; regular walls',len(walls),flush=True)
    result={'status':'PASS','independence':'No producer or inherited oracle imported; cyclotomic polynomial seed, Fraction-pair signs, symmetric-pair completeness, antipodal embedding count, sink DP','orders':receipt,'word_decisions':sum(r['word_decisions'] for r in receipt)}
    (root/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print('PASS all;',result['word_decisions'],'word decisions; seconds',round(time.time()-t,3),flush=True)
if __name__=='__main__':main()
