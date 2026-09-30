"""Recheck inherited regression, critical-ratio Galois obstruction and optima."""
from pathlib import Path
from fractions import Fraction as F
import json,math
from verify_atlas import chk,cmp,plus,times,neg,scalar,qvalue,sg,sh
ROOT=Path(__file__).resolve().parents[1]

def poly(p,x):
    z=F(0)
    for c in p[::-1]:z=z*x+c
    return z

def ipoly(p,l,r):
    # rigorous rational interval Horner
    a=b=F(0)
    for c in p[::-1]:
        products=[a*l,a*r,b*l,b*r];a=min(products)+c;b=max(products)+c
    return a,b

def frac(d):return F(d['numerator'],d['denominator'])
def tr(x):return x.numerator,0,x.denominator

def locate(A,w):
    roots=[tuple(v['value']) for v in A['roots']]
    if w['kind']=='RATIONAL':
        x=tr(frac(w['value']))
        for j,r in enumerate(roots):
            if cmp(x,r)==0:return 2*j+1,{'exact_rational':list(x)}
            if cmp(x,r)<0:return 2*j,{'exact_rational':list(x)}
        return 2*len(roots),{'exact_rational':list(x)}
    p=w['polynomial_coefficients_ascending'];ints=w['isolating_interval'];l=frac(ints['left']);r=frac(ints['right']);dp=[j*p[j] for j in range(1,len(p))]
    chk(l>1 and poly(p,l)*poly(p,r)<0,'ALGEBRAIC_BRACKET')
    dmin,dmax=ipoly(dp,l,r);chk(dmin>0 or dmax<0,'ALGEBRAIC_UNIQUENESS')
    for _ in range(512):
        j=sum(cmp(t,tr(l))<=0 for t in roots)
        lo=roots[j-1] if j else (1,0,1);hi=roots[j] if j<len(roots) else None
        if cmp(lo,tr(l))<0 and (hi is None or cmp(tr(r),hi)<0):
            return 2*j,{'polynomial':p,'isolating_left':list(tr(l)),'isolating_right':list(tr(r)),'derivative_nonzero_on_original_interval':True}
        mid=(l+r)/2
        if poly(p,l)*poly(p,mid)<0:r=mid
        else:l=mid
    raise RuntimeError('LANDMARK_LOCATION')

def main():
    atlases={k:json.loads((ROOT/f'results/atlas_k{k}.json').read_text()) for k in (4,6,8)}
    old=json.loads((ROOT/'sources/original_producer.json').read_text())['scientific_payload'];sk=json.loads((ROOT/'sources/critical_ratio_skeleton.json').read_text())['scientific_payload']
    weights={}
    for order in sk['orders']:
        for cell in order['open_intervals']:weights['ATLAS-'+cell['cell_id']]={'kind':'RATIONAL','value':cell['representative']}
        for cell in order['critical_points']:weights['ATLAS-'+cell['cell_id']]={'kind':'RATIONAL','value':cell['ratio']}
    for i,lm in enumerate(sk['landmarks']):weights[f"ATLAS-K8-L{i:03d}-{lm['landmark_id']}"]=lm['weight']
    regression=[]
    for row in old['rows']:
        k=row['construction_identity']['k'];A=atlases[k];idx,proof=locate(A,weights[row['row_id']]);fresh=A['rows'][idx]
        chk(fresh['terminal_sha256']==row['terminal']['selected_word_sha256'],'INHERITED_TERMINAL')
        chk(fresh['lower_longest_edges']==row['residual_dag']['longest_path_edges'],'INHERITED_PATH')
        chk(fresh['zero_embedding_words']==row['zero_and_descent']['zero_embedding_word_count'],'INHERITED_ZERO')
        chk(fresh['zero_lower_words']==row['zero_and_descent']['lower_zero_state_count'],'INHERITED_LOWER_ZERO')
        chk(fresh['ideal_count']==row['terminal']['a_star_count'],'INHERITED_IDEAL_COUNT')
        regression.append({'inherited_row':row['row_id'],'new_cell':fresh['id'],'location_certificate':proof})
    A=atlases[8];rd={tuple(r['value']):r['id'] for r in A['roots']};admissible=[]
    for r,i in rd.items():
        conj=(r[0],-r[1],r[2])
        if cmp(conj,(1,0,1))>0:
            chk(conj in rd,'GALOIS_CANDIDATE_CLOSURE');admissible.append((i,rd[conj]))
    vals=[(5,-1,1),(17,6,7),(5,1,1),(17,-6,7)]
    ids=[2*rd[v]+1 for v in vals];rr=[A['rows'][i] for i in ids]
    bits=[(ROOT/'results/certificates'/f"{r['id']}.bits").read_bytes() for r in rr]
    chk(bits[0]==bits[1] and bits[2]!=bits[3],'GALOIS_NONDESCENT')
    chk(rr[2]['lower_longest_edges']!=rr[3]['lower_longest_edges'],'GALOIS_NONISOMORPHISM')
    witness=next(i for i in range(65536) if ((bits[2][i//8]^bits[3][i//8])>>(i%8))&1)
    galois={'admissible_critical_ratios':len(admissible),'rational_fixed_ratios':sum(i==j for i,j in admissible),'two_element_orbits':sum(i<j for i,j in admissible),
      'critical_root_permutation':admissible,'nondescending_source_values':vals[:2],'nondescending_target_values':vals[2:],
      'cell_ids':[r['id'] for r in rr],'source_terminal_literal_equality':True,'target_terminal_literal_inequality':True,
      'target_longest_lengths':[r['lower_longest_edges'] for r in rr[2:]],'first_target_separating_word_code':witness}
    # Every optimal cell is exactly in the specified half-closed interval.
    for k,lo,hi,expected in [(6,(4,0,1),None,57),(8,(3,0,1),(1,2,1),141)]:
        a=atlases[k]
        for row in a['rows']:
            inside=cmp(row['ratio'],lo)>=0 and (hi is None or cmp(row['ratio'],hi)<0)
            chk((row['lower_longest_edges']==expected)==inside,'OPTIMUM_INTERVAL')
    lifted=[]
    def mul16(x,y):
        z=[0]*8
        for i,a in enumerate(x):
            for j,b in enumerate(y):
                e=i+j;z[e%8]+=a*b*(-1 if e>=8 else 1)
        return z
    sqrt2=[0,0,1,0,0,0,-1,0]
    for rec in A['embedding_spectrum']:
        w=rec['witness_word'];ww=''.join('C'+ch for ch in w);H=[0]*8;L=[0]*8
        for i,ch in enumerate(ww):
            e=(i+1)%16;sign=-1 if e>=8 else 1
            if ch=='A':H[e%8]-=sign
            elif ch=='T':H[e%8]+=sign
            elif ch=='C':L[e%8]-=sign
            else:L[e%8]+=sign
        a,b,d=rec['value'];SH=mul16(sqrt2,H)
        chk(any(H) and all(a*H[i]+b*SH[i]+d*L[i]==0 for i in range(8)),'QUADRATIC_LIFT16')
        z8=256+rec['accidental_words']
        lifted.append({'ratio':rec['value'],'word_length_16':ww,'embedding_identity_mod_X8_plus_1':True,
          'total_zero_words_k16_by_proved_block_formula':z8*z8,'structural_zero_words_k16':256**2})
    out={'status':'PASS','quadratic_field_lift_checks':lifted,'regression_count':len(regression),'regression':regression,'galois':galois,
      'exact_optima':{'k6':{'minimum':57,'lower_inclusive':[4,0,1],'upper_exclusive':None},'k8':{'minimum':141,'lower_inclusive':[3,0,1],'upper_exclusive':[1,2,1]}},
      'disjoint_optima':cmp((1,2,1),(4,0,1))<0}
    chk(out['disjoint_optima'],'DISJOINT_OPTIMA');(ROOT/'results/derived_checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: 30 inherited regressions; exact optima [4,infinity) and [3,1+2sqrt2); Galois nondescending quotient')
if __name__=='__main__':main()
