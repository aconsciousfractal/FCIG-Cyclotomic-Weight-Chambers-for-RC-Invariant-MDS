"""Exact DNA weight chambers at k=4,6,8. No external dependencies.
Q(sqrt(2)) arithmetic uses integer triples (a,b,d) for (a+b*sqrt(2))/d.
This producer uses the two-family selector; the verifier uses cyclotomic
polynomial coordinates and the original negative-ray seed independently.
"""
from __future__ import annotations
from collections import deque, Counter
from functools import cmp_to_key
from itertools import product
from pathlib import Path
from array import array
import json, hashlib, math, time, sys

V=((-1,0),(0,-1),(0,1),(1,0))
# 2*cos(e*2pi/k), and sin(e*2pi/k) divided by its positive common factor.
# Common factors: 1 (k4), sqrt3/2 (k6), sqrt2/2 (k8).
TRIG={
 4:(((2,0),(0,0),(-2,0),(0,0),(2,0)),((0,0),(1,0),(0,0),(-1,0),(0,0))),
 6:(((2,0),(1,0),(-1,0),(-2,0),(-1,0),(1,0),(2,0)),((0,0),(1,0),(1,0),(0,0),(-1,0),(-1,0),(0,0))),
 8:(((2,0),(0,1),(0,0),(0,-1),(-2,0),(0,-1),(0,0),(0,1),(2,0)),((0,0),(1,0),(0,1),(1,0),(0,0),(-1,0),(0,-1),(-1,0),(0,0)))
}

def require(b,code):
    if not b:raise RuntimeError(code)

def canonical(data):return json.dumps(data,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
def digest(data):return hashlib.sha256(canonical(data)).hexdigest()
def filehash(data):return hashlib.sha256(data).hexdigest()
def signpair(a:int,b:int)->int:
    if not a:return (b>0)-(b<0)
    if not b:return (a>0)-(a<0)
    sa=(a>0)-(a<0);sb=(b>0)-(b<0)
    if sa==sb:return sa
    return sa*((a*a>2*b*b)-(a*a<2*b*b))
def normalize(a,b,d):
    require(d!=0,'ZERO_DENOMINATOR')
    if d<0:a,b,d=-a,-b,-d
    g=math.gcd(math.gcd(abs(a),abs(b)),d)
    return a//g,b//g,d//g
def compare(x,y):return signpair(x[0]*y[2]-y[0]*x[2],x[1]*y[2]-y[1]*x[2])
def root(f):
    a,b,c,d=f
    if a==c==0:return None
    return normalize(-a*b+2*c*d,c*b-a*d,a*a-2*c*c)
def evaluate(f,r):
    a,b,c,d=f;x,y,z=r
    return signpair(a*x+b*z+2*c*y,a*y+c*x+d*z)
def rational_between(left,right):
    # dyadic bisection with integer sign comparisons only
    lo=(1,0,1);hi=(2,0,1)
    if right is None:
        while compare(hi,left)<=0:hi=(hi[0]*2,0,1)
        return hi
    while compare(hi,right)<=0:hi=(hi[0]*2,0,1)
    for _ in range(256):
        mid=normalize(lo[0]*hi[2]+hi[0]*lo[2],0,2*lo[2]*hi[2])
        if compare(left,mid)<0 and compare(mid,right)<0:return mid
        if compare(mid,left)<=0:lo=mid
        else:hi=mid
    raise RuntimeError('DYADIC_SAMPLE_NOT_FOUND')
def word(code,n):
    out=[0]*n
    for i in range(n-1,-1,-1):out[i]=code%4;code//=4
    return tuple(out)
def wordtext(code,n):return ''.join('ACGT'[c] for c in word(code,n))
def revcomp(code,n):
    z=0
    for _ in range(n):z=4*z+3-code%4;code//=4
    return z
def form(w,coeff):
    a=b=c=d=0
    for ch,(x,y) in zip(w,coeff):
        h,l=V[ch];a+=h*x;b+=l*x;c+=h*y;d+=l*y
    return a,b,c,d

def rotation_data(k):
    n=4**k;mask=4**(k-1);owner=[-1]*n;pcrs=[]
    for code in range(n):
        if owner[code]>=0:continue
        cur=code;cyc=[]
        while not cyc or cur!=code:
            cyc.append(cur);cur=(cur%mask)*4+cur//mask
        idx=len(pcrs)
        for c in cyc:owner[c]=idx
        pcrs.append(cyc)
    b=bytearray(n);visited=set();rc=[revcomp(c,k) for c in range(n)]
    for idx,cyc in enumerate(pcrs):
        if idx in visited:continue
        mate=owner[rc[cyc[0]]];visited.add(idx);visited.add(mate)
        if mate==idx:
            fixed=[c for c in cyc if rc[c]==c]
            require(len(fixed)==2,'ZERO_PCR_FIXED_RULE');b[min(fixed)]=1
        else:
            x=min(min(cyc),min(pcrs[mate]));b[x]=b[rc[x]]=1
    return owner,pcrs,b,rc

def build_catalog(k):
    C,S=TRIG[k];lw=list(product(range(4),repeat=k-1))
    lf=[form(w,S[1:k]) for w in lw]
    rf=[form(w,[(x-2,y) for x,y in C[1:k]]) for w in lw]
    roots={};families={};counts={}
    def collect(name,forms):
        rset=set()
        for code,f in enumerate(forms):
            r=root(f)
            if r is not None and compare(r,(1,0,1))>0:
                rset.add(r);roots.setdefault(r,{}).setdefault(name,code)
        families[name]=rset;counts[name]={'word_occurrences':len(forms),'distinct_forms':len(set(forms)),'structural_occurrences':forms.count((0,0,0,0)),'positive_roots':len(rset)}
    collect('lower_lambda',lf);collect('midpoint',rf)
    wf_re=[];wf_im=[];emb={};struct=0
    for w in product(range(4),repeat=k):
        re=form(w,C[1:]);im=form(w,S[1:]);wf_re.append(re);wf_im.append(im)
        if re==im==(0,0,0,0):struct+=1;continue
        r=root(re) or root(im)
        if r and compare(r,(1,0,1))>0 and evaluate(re,r)==evaluate(im,r)==0:
            if r not in emb:emb[r]={'accidental_words':0,'witness_code':len(wf_re)-1}
            emb[r]['accidental_words']+=1
    collect('seed_real',wf_re);collect('seed_imag',wf_im)
    ordered=sorted(roots,key=cmp_to_key(compare));root_records=[]
    for i,r in enumerate(ordered):
        act=next((c for c in range(len(lw)) if rf[c]!=(0,0,0,0) and evaluate(rf[c],r)==0 and evaluate(lf[c],r)>0),None)
        root_records.append({'id':i,'value':r,'families':sorted(roots[r]),'active_midpoint_word_code':act,'active_midpoint_word':wordtext(act,k-1) if act is not None else None})
    return lw,lf,rf,ordered,{'k':k,'roots':root_records,'family_counts':counts,'structural_embedding_words':struct,'embedding_spectrum':[{'value':r,**emb[r],'witness_word':wordtext(emb[r]['witness_code'],k)} for r in sorted(emb,key=cmp_to_key(compare))]}

def dag_data(n,terminal,colours):
    indeg=[0]*n;adj=[[] for _ in range(n)]
    for c,m in enumerate(terminal):
        if not m:
            u,v=c//4,c%n
            require(colours[u]<=colours[v],'RESIDUAL_COLOUR_ORDER')
            adj[u].append(v);indeg[v]+=1
    q=deque(i for i in range(n) if not indeg[i]);dist=[0]*n;pred=[-1]*n;local=[0]*n;order=[]
    while q:
        u=q.popleft();order.append(u)
        for v in adj[u]:
            if dist[v]<dist[u]+1:dist[v]=dist[u]+1;pred[v]=u
            if colours[u]==colours[v]:local[v]=max(local[v],local[u]+1)
            indeg[v]-=1
            if not indeg[v]:q.append(v)
    require(len(order)==n,'RESIDUAL_CYCLE')
    best=max(range(n),key=lambda v:dist[v]);path=[best]
    while pred[path[-1]]>=0:path.append(pred[path[-1]])
    path.reverse();L=dist[best]
    require(len(path)==L+1,'PATH_LENGTH')
    lengths={str(c):max((local[v] for v in range(n) if colours[v]==c),default=0) for c in (-1,0,1)}
    require(lengths['-1']==lengths['1'],'RC_BLOCK_LENGTH')
    return dist,path,lengths,len(adj),sum(len(a) for a in adj)

def packbits(data):
    out=bytearray((len(data)+7)//8)
    for i,x in enumerate(data):
        if x:out[i//8]|=1<<(i%8)
    return bytes(out)

def compute_row(k,r,lw,lf,rf,bd,cell_id,kind,outdir):
    owner,pcrs,b,rc=bd;n=len(lw)
    sl=[evaluate(f,r) for f in lf];sr=[evaluate(f,r) for f in rf]
    lrc=[revcomp(c,k-1) for c in range(n)]
    a=[int(sl[i]>0 and (sr[i]<0 or (sr[i]==0 and i<lrc[i]))) for i in range(n)]
    colours=[0 if sl[i]<=0 else (1 if a[i] else -1) for i in range(n)]
    term=bytearray(4*n);zero=0;pcrsum=[0]*len(pcrs)
    for c in range(4*n):
        u,v=c//4,c%n
        if sl[u]==sl[v]==0:m0=b[c];zero+=1
        else:m0=int(sl[u]<=0<sl[v])
        m=m0+a[u]-a[v];require(m in (0,1),'NONBOOLEAN')
        term[c]=m;pcrsum[owner[c]]+=m
    require(all(s==1 for s in pcrsum),'ONE_PER_PCR')
    require(all(term[c]==term[rc[c]] for c in range(4*n)),'RC_TERMINAL')
    dist,path,bl,_,arc_count=dag_data(n,term,colours)
    L=max(dist);upper=2*bl['1']+bl['0']+2
    require(L<=upper,'THREE_BLOCK_UPPER_BOUND')
    bits=packbits(term);levels=array('H',dist)
    if sys.byteorder!='little':levels.byteswap()
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/(cell_id+'.bits')).write_bytes(bits)
    (outdir/(cell_id+'.levels')).write_bytes(levels.tobytes())
    selected=[i for i,x in enumerate(term) if x]
    row={'id':cell_id,'kind':kind,'ratio':r,'terminal_sha256':digest(selected),'bitset_sha256':filehash(bits),
      'level_sha256':filehash(levels.tobytes()),'selected_count':sum(term),'zero_embedding_words':zero,
      'zero_lower_words':sl.count(0),'positive_lower_words':sum(x>0 for x in sl),'ideal_count':sum(a),
      'midpoint_active_ties':sum(sl[i]>0 and sr[i]==0 for i in range(n)),
      'midpoint_active_accidental_ties':sum(sl[i]>0 and sr[i]==0 and rf[i]!=(0,0,0,0) for i in range(n)),
      'lower_longest_edges':L,'upper_longest_edges':L-1,'block_lengths':bl,'block_upper_bound':upper,
      'block_bound_saturated':L==upper,'residual_arcs':arc_count,'longest_path_lower_vertex_codes':path}
    return row,term

def build(k,out):
    lw,lf,rf,roots,cat=build_catalog(k);bd=rotation_data(k);rows=[];terms=[]
    for j in range(2*len(roots)+1):
        kind='point' if j%2 else 'interval'
        if j%2:r=roots[j//2]
        else:r=rational_between(roots[j//2-1] if j else (1,0,1),roots[j//2] if j//2<len(roots) else None)
        cid=f'k{k}_c{j:03d}';row,term=compute_row(k,r,lw,lf,rf,bd,cid,kind,out/'certificates')
        rows.append(row);terms.append(term)
    # Full equality classes are compared as literal byte vectors, not via hash alone.
    classes={}
    for row,t in zip(rows,terms):
        key=bytes(t)
        if key not in classes:classes[key]=[]
        classes[key].append(row['id'])
    for i,rr in enumerate(cat['roots']):
        left,point,right=terms[2*i:2*i+3]
        def delta(a,b):
            ds=[j for j in range(len(a)) if a[j]!=b[j]]
            return {'symmetric_difference':len(ds),'first_witness_code':ds[0] if ds else None,
              'first_witness_word':wordtext(ds[0],k) if ds else None}
        rr['left_to_right']=delta(left,right);rr['left_to_point']=delta(left,point);rr['point_to_right']=delta(point,right)
        rr['essential']=left!=point or point!=right
        rr['point_only_anomaly']=left==right and left!=point
        rr['three_distinct']=left!=point and point!=right and left!=right
    best=min(r['lower_longest_edges'] for r in rows);worst=max(r['lower_longest_edges'] for r in rows)
    cat.update({'status':'PRODUCED_EXACT','rows':rows,'terminal_classes':list(classes.values()),'distinct_terminals':len(classes),
      'open_intervals':len(roots)+1,'point_strata':len(roots),'essential_roots':sum(r['essential'] for r in cat['roots']),
      'optimum_lower_longest_edges':best,'optimum_cells':[r['id'] for r in rows if r['lower_longest_edges']==best],
      'maximum_lower_longest_edges':worst,'length_histogram':dict(sorted(Counter(r['lower_longest_edges'] for r in rows).items())),
      'block_bound_nonsaturated_cells':[r['id'] for r in rows if not r['block_bound_saturated']],
      'one_per_pcr':len(bd[1]),'scope':'DNA s=(-rho,-1,1,rho), rho>1, fixed zero-PCR and lex midpoint ties; exact lower-edge path units'})
    (out/f'atlas_k{k}.json').write_text(json.dumps(cat,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PRODUCED',k,'roots',len(roots),'cells',len(rows),'terminals',len(classes),'min',best,'max',worst,'nonsaturated',cat['block_bound_nonsaturated_cells'],flush=True)
    return cat

def main():
    out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'results'
    out.mkdir(parents=True,exist_ok=True);t=time.time()
    for k in (4,6,8):build(k,out)
    print('seconds',round(time.time()-t,3),flush=True)
if __name__=='__main__':main()
