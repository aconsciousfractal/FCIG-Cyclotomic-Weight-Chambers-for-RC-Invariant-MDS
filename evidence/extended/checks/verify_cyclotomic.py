"""Reviewer checks of bounded positional witnesses and k=30 counts.
No import from the external producer. SymPy supplies cyclotomic polynomials;
monic integer division and ordered-letter counting are implemented here.
"""
from pathlib import Path
from collections import Counter
from itertools import product
from fractions import Fraction
from functools import lru_cache
import argparse,json
import sympy as sp

LETTERS={'A':(-1,0),'C':(0,-1),'G':(0,1),'T':(1,0)}
def need(x,msg):
    if not x:raise ValueError(msg)
@lru_cache(None)
def phi(k):
    x=sp.Symbol('x');return list(reversed([int(v) for v in sp.Poly(sp.cyclotomic_poly(k,x),x).all_coeffs()]))
def rem(a,k):
    f=phi(k);d=len(f)-1;a=list(a)
    for i in range(len(a)-1,d-1,-1):
        c=a[i]
        if c:
            for j,v in enumerate(f):a[i-d+j]-=c*v
    return (a+[0]*d)[:d]
def parts(word):return [[LETTERS[c][i] for c in word] for i in (0,1)]
def count_rectangle(rho=None):
    units=list(LETTERS.values());d=Counter()
    for a,b in product(units,repeat=2):
        v=(a[0]-b[0],a[1]-b[1])
        d[v if rho is None else (rho.numerator*v[0]+rho.denominator*v[1],)]+=1
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    add=lambda a,b:tuple(x+y for x,y in zip(a,b))
    shifts={sub(a,b) for a,b in product(d,repeat=2)};total=0
    # Transpose the 3 by 5 CRT rectangle: explicit first row of width 3,
    # four further independent rows; no histogram producer imported.
    for row in product(d,repeat=3):
        first=1
        for v in row:first*=d[v]
        translate=0
        for shift in shifts:
            ways=1
            for v in row:ways*=d.get(add(shift,v),0)
            translate+=ways
        total+=first*translate**4
    return total
def main():
    ap=argparse.ArgumentParser();ap.add_argument('round1',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    src=json.loads((args.round1/'results/new_results.json').read_text());rational=[];algebraic=[]
    for rec in src['rational_witnesses']:
        word=rec['word_exponent_order'];k=len(word);H,L=parts(word);rho=Fraction(rec['rho'])
        need(any(rem(H,k)) and not any(rem([rho.numerator*h+rho.denominator*l for h,l in zip(H,L)],k)),'rational witness')
        rational.append({'k':k,'rho':str(rho)})
    for rec in src['algebraic_witnesses']:
        word=rec['word_exponent_order'];k=len(word);H,L=parts(word);n=rec['verification']['base_order']
        U={int(i)*(k//n):v for i,v in rec['U'].items()};a=[0]*(k+max(U))
        for i,h in enumerate(H):
            for j,u in U.items():a[i+j]+=h*u
        for i,l in enumerate(L):a[i]+=rec['denominator']*l
        need(any(rem(H,k)) and not any(rem(a,k)),'algebraic witness')
        mutated=a[:];mutated[0]+=1;need(any(rem(mutated,k)),'negative polynomial control')
        algebraic.append({'label':rec['label'],'k':k,'identity_verified':True})
    expected={'structural':21566097316,'5/3':21573498016,'2':27518030956,'3':288112407796,'5':22778390836,'7':21566157256}
    counts={key:count_rectangle(None if key=='structural' else Fraction(key)) for key in expected}
    need(counts==expected,'k30 counts')
    def template(word,k,role):
        a=[0]*k
        for e,ch in enumerate(word,1):
            h=LETTERS[ch][0];a[e%k]+=h
            a[(-e)%k]+=h if role=='midpoint' else -h
            if role=='midpoint':a[0]-=2*h
        return a
    coefficient=template('TTT',4,'midpoint')[0]
    need(coefficient==-6,'role height counterexample')
    need(not any(rem(template('TTTAATA',8,'midpoint'),8)),'structural midpoint')
    need(rem(template('TTTAATA',8,'lambda'),8)==[0,4,0,4],'positive lambda identity')
    result={'status':'PASS','rational_positional_witnesses':rational,'algebraic_positional_witnesses':algebraic,'k30_independent_transposed_counts':counts,
      'counterexamples_checked':{'midpoint_TTT_k4_constant_coefficient':coefficient,'TTTAATA_k8_midpoint_structural':True,'TTTAATA_k8_lambda':'2*sqrt(2)*rho'},
      'scope':'Finite exact positional identities and counts; uniform necessity and sufficiency checked in written mathematical review, not proved by these tests.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':'PASS','rational_witnesses':len(rational),'algebraic_witnesses':algebraic,'counts':counts}))
if __name__=='__main__':main()
