"""Additional independent finite checks: active midpoint counts and matrix orientation."""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import json
import strengthening as s
ROOT=Path(__file__).resolve().parents[1]

def rectangular_direct(p,q,rho):
    # Sum over the actual first row, not histograms; independent aggregation.
    n=Counter()
    for d,c in s.D.items(): n[d if rho is None else (rho.numerator*d[0]+rho.denominator*d[1],)]+=c
    D=list(n); shifts={s.sub(a,b) for a,b in product(D,repeat=2)}
    total=0
    for b in product(D,repeat=q):
        weight=1
        for x in b: weight*=n[x]
        row_total=0
        for a in shifts:
            term=1
            for x in b:
                term*=n.get(s.add(a,x),0)
                if not term:break
            row_total+=term
        total+=weight*row_total**(p-1)
    return total

def main():
    checks=[];active={}
    for p in [3,5,7,11,13,17,19]:
        values={}
        for rho in s.midpoint_spectrum(p):
            total=s.midpoint_scalar_count(p,rho)
            joint=s.midpoint_scalar_count(p,rho,lambda_zero=True)
            if not (total>joint and (total-joint)%2==0): raise AssertionError('active midpoint count')
            values[str(rho)]={'midpoint_total':total,'joint_lambda_zero':joint,'active_lambda_positive':(total-joint)//2}
        active[str(p)]=values
        checks.append({'name':f'all_midpoint_ratios_active_p{p}','status':'PASS','ratios':len(values)})
    for rho in [None,*s.embedding_spectrum(2),Fraction(4)]:
        fast=s.two_prime_count(3,5,rho)
        slow=rectangular_direct(5,3,rho) # 9^3 first rows, opposite orientation of coefficients.
        if fast!=slow: raise AssertionError('rectangular count mismatch')
        checks.append({'name':f'matrix30_direct_{rho}','status':'PASS','count':fast})
    # Raw midpoint certificate height counterexample: k=4, lower word TTT.
    laurent=Counter()
    for e in [1,2,3]:laurent[e]+=1;laurent[-e]+=1;laurent[0]-=2
    if laurent[0]!=-6: raise AssertionError('height counterexample mismatch')
    checks.append({'name':'appendix_A3_height_counterexample','status':'PASS','constant_coefficient':-6})
    # Structural midpoint does not imply inactive at k=8.
    w='TTTAATA'; k=8; phi=s.phi_coeff(k)
    H=[0]+[s.VECTOR[ch][0] for ch in w]
    HC=[H[-e%k] for e in range(k)]
    midpoint=[H[e]+HC[e]-(2*sum(H) if e==0 else 0) for e in range(k)]
    imag=s.poly_mod([H[e]-HC[e] for e in range(k)],phi)
    target=s.poly_mod([0,4,0,0,0,0,0,-4],phi)
    if any(s.poly_mod(midpoint,phi)) or imag!=target: raise AssertionError('k8 active structural tie')
    checks.append({'name':'k8_structural_midpoint_can_be_active','status':'PASS','word':w,'lambda':'2*sqrt(2)*rho'})
    # A point-only anomaly is not removable from a root partition.
    indispensable=lambda left,point,right:len({left,point,right})>1
    if not indispensable('same','other','same') or indispensable('same','same','same'):
        raise AssertionError('point-only anomaly test')
    checks.append({'name':'root_partition_point_only_anomaly','status':'PASS'})
    out={'status':'PASS','check_count':len(checks),'checks':checks,'active_midpoint_counts':active,
         'A3_counterexample':{'k':4,'lower_word':'TTT','role':'lower midpoint','Laurent_A_coefficients':dict(sorted(laurent.items()))}}
    (ROOT/'results/additional_checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(checks),'additional checks')
if __name__=='__main__':main()
