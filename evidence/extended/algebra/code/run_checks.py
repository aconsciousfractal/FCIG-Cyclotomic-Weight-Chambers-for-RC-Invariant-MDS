"""Run deterministic checks and save a machine-readable scientific receipt."""
from __future__ import annotations
from fractions import Fraction
from itertools import product, permutations
from math import gcd
from pathlib import Path
import hashlib,json,time
import sympy as sp
import strengthening as s

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
checks=[]
def check(name,condition,**details):
    if not condition: raise AssertionError(name)
    checks.append({'name':name,'status':'PASS',**details})

def main():
    start=time.time()
    spectra={}
    for t in range(6):
        S=s.embedding_spectrum(t)
        Farey=sorted({Fraction(v+u,v-u) for v in range(2,2**t+1) for u in range(1,v) if gcd(u,v)==1})
        check(f'Farey_spectrum_t{t}',S==Farey,count=len(S))
        check(f'Farey_count_t{t}',len(S)==sum(int(sp.totient(i)) for i in range(2,2**t+1)))
        spectra[str(t)]=[str(x) for x in S]
    witnesses=[]
    for primes in ([3],[5],[3,5],[3,7],[5,7],[3,5,7]):
        for rho in s.embedding_spectrum(len(primes)):
            w=s.embedding_witness(list(primes),rho)
            v=s.verify_word(w,rho)
            check(f'word_{len(w)}_{rho}',v['scalar_zero'] and v['nonstructural'])
            witnesses.append({'primes':primes,'rho':str(rho),'word_exponent_order':w,'verification':v})
    for r in [2,3,5,6,10,15]:
        w=s.embedding_witness([3,5],Fraction(2),r)
        v=s.verify_word(w,Fraction(2))
        check(f'power_lift_{len(w)}',v['scalar_zero'] and v['nonstructural'])
    for bad in [Fraction(4),Fraction(3,2),Fraction(9),Fraction(7,5)]:
        try: s.embedding_witness([3,5],bad)
        except ValueError: check(f'reject_unrealizable_k30_{bad}',True)
        else: check(f'reject_unrealizable_k30_{bad}',False)
    dp={}
    for k in [4,6,10]:
        for role in ['embedding','real','imag','lower','midpoint']:
            data=s.bounded_role_dp(k,role); dp[f'{k}_{role}']=data
            if role=='embedding':
                t=len([p for p in sp.factorint(k) if p!=2])
                check(f'dp_embedding_spectrum_{k}',set(data['accidental'])==set(spectra[str(t)]))
            if k in (6,10):
                p=k//2; h=(p-1)//2
                if role=='real':
                    st=s.rational_seed_real_count(p)
                    check(f'seed_real_struct_{k}',data['structural']==st)
                    check(f'seed_real_spectrum_{k}',set(data['accidental'])==set(spectra['2']))
                    for rho in s.embedding_spectrum(2):
                        check(f'seed_real_count_{k}_{rho}',data['accidental'][str(rho)]+st==s.rational_seed_real_count(p,rho))
                if role=='embedding':
                    M=4**p+4+4*2**p; M3=4**p+2+2*2**p+2*3**p
                    check(f'prime_embedding_counts_{k}',data['structural']==M and data['accidental']=={'3':M3-M})
                if role=='midpoint':
                    check(f'midpoint_struct_{k}',data['structural']==4)
                    check(f'midpoint_spectrum_{k}',set(data['accidental'])==set(map(str,s.midpoint_spectrum(p))))
                    for rho in s.midpoint_spectrum(p):
                        check(f'midpoint_count_{k}_{rho}',data['accidental'][str(rho)]+4==s.midpoint_scalar_count(p,rho))
    # Independent scalar-transfer search of every ratio allowed by the coefficient norm bound.
    mid={}
    for p in [3,5,7,11,13]:
        candidates={Fraction(a,b) for a in range(2,4*p+5) for b in range(1,a) if a+b<=4*p+4}
        observed={r:s.midpoint_scalar_count(p,r) for r in candidates}
        roots=sorted(r for r,n in observed.items() if n>4)
        check(f'midpoint_complete_bound_scan_p{p}',roots==s.midpoint_spectrum(p),ratios_tested=len(candidates),resonances=len(roots))
        h=(p-1)//2
        check(f'midpoint_max_p{p}',max(roots)==2*p+4-(-1)**h)
        mid[str(p)]={'rational_resonances':{str(r):observed[r] for r in roots},'coarse_candidates_checked':len(candidates),'maximum':str(max(roots))}
    counts={}
    for p,q in [(3,5),(3,7),(5,7)]:
        st=s.two_prime_count(p,q)
        table={'structural':st}
        for rho in [*s.embedding_spectrum(2),Fraction(4)]:
            total=s.two_prime_count(p,q,rho)
            check(f'matrix_count_{2*p*q}_{rho}',total>=st and ((total>st)==(rho in s.embedding_spectrum(2))))
            table[str(rho)]={'total':total,'accidental':total-st}
        k=2*p*q
        tail=sum(int(sp.totient(k//d))*4**d for d in sp.divisors(k) if d<k)
        check(f'Burnside_integrality_{k}',all((v['total']+tail)%k==0 for key,v in table.items() if key!='structural') and (st+tail)%k==0)
        table['structural_zero_PCR']=(st+tail)//k
        counts[str(k)]=table
    algebraic=[]
    for label,n,U,d,primes in [('sqrt2',8,{1:1,7:1},1,[3]),('gold',10,{0:2,1:3,9:3},1,[3,7,11])]:
        w=s.algebraic_lift(n,U,d,primes); data=s.verify_algebraic_lift(w,n,U,d)
        check(f'algebraic_lift_{label}',data['scalar_zero'] and data['component_H_nonzero'])
        algebraic.append({'label':label,'U':U,'denominator':d,'primes':primes,'word_exponent_order':w,'verification':data})
        # A single letter mutation must destroy the weighted zero: nonzero scalar
        # monomial, since distinct nonzero weights were used.
        mutated=('A' if w[0]!='A' else 'T')+w[1:]
        v=s.verify_algebraic_lift(mutated,n,U,d)
        check(f'reject_mutated_algebraic_{label}',not v['scalar_zero'])
    # Position permutations: brute check all permutations, binary words, k=2,3,4.
    for k in [2,3,4]:
        allowed=[]
        for pi in permutations(range(k)):
            ok=True
            for w in product(range(2),repeat=k):
                for x in range(2):
                    v=w[1:]+(x,); a=tuple(w[j] for j in pi); b=tuple(v[j] for j in pi)
                    if a[1:]!=b[:-1]: ok=False; break
                if not ok: break
            if ok: allowed.append(pi)
        check(f'directed_coordinate_permutations_k{k}',allowed==[tuple(range(k))])
    payload={'status':'PASS','scope':'Exact finite checks, not a formal proof assistant certificate',
             'check_count':len(checks),'checks':checks,'embedding_spectra':spectra,'positional_DP':dp,
             'midpoint_transfer':mid,'two_prime_counts':counts,'rational_witnesses':witnesses,'algebraic_witnesses':algebraic}
    raw=json.dumps(payload,indent=2,sort_keys=True)
    (OUT/'new_results.json').write_text(raw+'\n',encoding='utf-8',newline='\n')
    (OUT/'run_summary.json').write_text(json.dumps({'status':'PASS','check_count':len(checks),'seconds':time.time()-start,'new_results_sha256':hashlib.sha256((raw+'\n').encode()).hexdigest()},indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(checks),'checks',round(time.time()-start,2),'seconds',flush=True)

if __name__=='__main__': main()
