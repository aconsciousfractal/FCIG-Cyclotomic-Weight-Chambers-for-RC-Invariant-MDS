"""Exact finite checks for the accompanying mathematical supplement.

All DNA words returned here are indexed by exponent e=0,...,k-1.
To use the source paper's e=j+1 convention, left-rotate once: word[1:]+word[:1].
No numerical root approximations or floating-point decisions are used.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import factorial, gcd, prod
from functools import lru_cache
import sympy as sp

V = ((1,0),(-1,0),(0,1),(0,-1))
LETTER = {(1,0):'T',(-1,0):'A',(0,1):'G',(0,-1):'C'}
VECTOR = {v:k for k,v in LETTER.items()}
D = Counter((a[0]-b[0], a[1]-b[1]) for a,b in product(V,repeat=2))
PAIR = {(a[0]-b[0],a[1]-b[1]):(a,b) for a,b in product(V,repeat=2)}

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def scale(s,a): return tuple(s*x for x in a)

def embedding_spectrum(t: int) -> list[Fraction]:
    """Complete rational accidental embedding spectrum, t distinct odd primes."""
    if t < 0: raise ValueError('t must be nonnegative')
    M=2**(t+1)
    return sorted({Fraction(a,b) for b in range(1,M) for a in range(b+1,M-b+1)
                   if gcd(a,b)==1 and ((a+b)%2==0 or 2*(a+b)<=M)})

def unit_decomposition(v: tuple[int,int], length: int) -> list[tuple[int,int]]:
    used=abs(v[0])+abs(v[1])
    if used>length or (length-used)%2: raise ValueError('Infeasible unit-vector length')
    out=[]
    for i,c in enumerate(v):
        u=tuple((1 if c>0 else -1) if j==i else 0 for j in range(2))
        out.extend([u]*abs(c))
    out.extend([(1,0),(-1,0)]*((length-used)//2))
    return out

def embedding_witness(primes: list[int], rho: Fraction, multiplicity: int=1) -> str:
    """Position-aware witness for k=2*product(primes), optionally power-lifted.

    Distinct odd primes are required. The lift parameter must have no prime
    factors outside 2 and primes; thus Phi_{k*r}(X)=Phi_k(X^r).
    """
    if len(set(primes)) != len(primes) or any(p==2 or not sp.isprime(p) for p in primes):
        raise ValueError('distinct odd primes required')
    if rho<=1 or rho not in embedding_spectrum(len(primes)):
        raise ValueError('ratio outside proved spectrum')
    if multiplicity<1 or any(p not in [2,*primes] for p in sp.factorint(multiplicity)):
        raise ValueError('invalid power lift')
    t=len(primes); M=2**(t+1); a,b=rho.numerator,rho.denominator
    g=1 if (a+b)%2==0 else 2
    units=unit_decomposition((g*b,-g*a),M)
    corner={}
    for i,bits in enumerate(product((0,1),repeat=t)):
        corner[bits]=scale((-1)**sum(bits),add(units[2*i],units[2*i+1]))
    N=prod(primes); word=[None]*(2*N)
    for indices in product(*(range(p) for p in primes)):
        j=sum((N//p)*u for p,u in zip(primes,indices))%N
        bits=tuple(int(u==p-1) for u,p in zip(indices,primes))
        d=scale((-1)**j,corner[bits])
        word[j],word[j+N]=PAIR[d]
    if any(v is None for v in word): raise AssertionError('unfilled position')
    if multiplicity==1: return ''.join(LETTER[v] for v in word)
    out=[(1,0)]*(len(word)*multiplicity)
    for j,v in enumerate(word): out[j*multiplicity]=v
    return ''.join(LETTER[v] for v in out)

@lru_cache(None)
def phi_coeff(k: int) -> tuple[int,...]:
    z=sp.Symbol('z'); f=sp.Poly(sp.cyclotomic_poly(k,z),z)
    return tuple(int(f.nth(i)) for i in range(f.degree()+1))

def poly_mod(values: list[int] | tuple[int,...], phi: tuple[int,...]) -> tuple[int,...]:
    out=list(values); degree=len(phi)-1
    if len(out)<degree: out.extend([0]*(degree-len(out)))
    for j in range(len(out)-1,degree-1,-1):
        c=out[j]
        if c:
            for i,v in enumerate(phi): out[j-degree+i]-=c*v
    return tuple(out[:degree])

def verify_word(word: str, rho: Fraction) -> dict:
    k=len(word); phi=phi_coeff(k)
    hv=[VECTOR[ch][0] for ch in word]; lv=[VECTOR[ch][1] for ch in word]
    H=poly_mod(hv,phi); L=poly_mod(lv,phi)
    zero=not any(rho.numerator*a+rho.denominator*b for a,b in zip(H,L))
    return {'k':k,'ratio':str(rho),'scalar_zero':zero,'nonstructural':any(H) or any(L),
            'component_H_nonzero':any(H)}

def bounded_role_dp(k: int,role: str,state_cap: int=4_000_000) -> dict:
    """Independent positional cyclotomic remainder DP (compressed words)."""
    if role not in ('embedding','real','imag','lower','midpoint'): raise ValueError(role)
    phi=phi_coeff(k); dim=len(phi)-1
    powers=[poly_mod([0]*e+[1],phi) for e in range(k)]
    exponents=range(1,k) if role in ('lower','midpoint') else range(k)
    states={(0,)*(2*dim):1}
    for e in exponents:
        a,b=powers[e],powers[-e%k]
        if role=='embedding': v=a
        elif role=='real': v=add(a,b)
        elif role in ('imag','lower'): v=sub(a,b)
        else: v=sub(add(a,b),scale(2,powers[0]))
        changes=(v+(0,)*dim,scale(-1,v)+(0,)*dim,(0,)*dim+v,(0,)*dim+scale(-1,v))
        nxt=defaultdict(int)
        for s,n in states.items():
            for change in changes: nxt[add(s,change)]+=n
        states=nxt
        if len(states)>state_cap: raise RuntimeError('state cap exceeded')
    structural=0; accidental=Counter()
    for state,n in states.items():
        a,b=state[:dim],state[dim:]
        if not any(a):
            if not any(b): structural+=n
            continue
        i=next(i for i,x in enumerate(a) if x)
        rho=Fraction(-b[i],a[i])
        if rho>1 and not any(rho*x+y for x,y in zip(a,b)): accidental[rho]+=n
    return {'structural':structural,'accidental':{str(r):accidental[r] for r in sorted(accidental)},'states':len(states)}

def two_prime_count(p: int,q: int,rho: Fraction|None=None) -> int:
    """Weighted rectangular-kernel sum, by first-row histograms.

    None = vector/structural, rational rho = scalar total embedding zeros.
    Complexity polynomial in q for a fixed difference alphabet; use the
    smaller of the two primes as the histogram length.
    """
    if p==q or min(p,q)<3 or not sp.isprime(p) or not sp.isprime(q):
        raise ValueError('two distinct odd primes required')
    if p<q: p,q=q,p
    counts=Counter()
    for d,n in D.items():
        key=d if rho is None else (rho.numerator*d[0]+rho.denominator*d[1],)
        counts[key]+=n
    keys=list(counts); shifts={sub(x,y) for x,y in product(keys,repeat=2)}
    dlen=len(keys)
    # Precompute weights of all allowed translated entries.
    translations=[[counts.get(add(a,b),0) for b in keys] for a in shifts]
    fact=[factorial(i) for i in range(q+1)]
    total=0
    def visit(i,left,hist,den,weight):
        nonlocal total
        if i==dlen-1:
            h=hist+[left]; w=weight*counts[keys[i]]**left
            rows=0
            nz=[(j,c) for j,c in enumerate(h) if c]
            for weights in translations:
                pr=1
                for j,c in nz:
                    if not weights[j]: pr=0; break
                    pr*=weights[j]**c
                rows+=pr
            total+=fact[q]//(den*fact[left])*w*rows**(p-1)
            return
        for n in range(left+1): visit(i+1,left-n,hist+[n],den*fact[n],weight*counts[keys[i]]**n)
    visit(0,q,[],1,1)
    return total

def rational_seed_real_count(p: int,rho: Fraction|None=None) -> int:
    if p<3 or not sp.isprime(p): raise ValueError('odd prime required')
    h=(p-1)//2; base=4*36**h+8*6**h+4
    if rho==Fraction(5,3): return 4*36**h+4*10**h+4*6**h+2*5**h+2
    if rho==2: return 4*36**h+4*22**h+4*6**h+2*17**h+2
    if rho==3: return 4*44**h+6*31**h+4*10**h+2
    if rho==5: return 4*36**h+4*10**h+4*6**h+2*25**h+2
    if rho==7: return 4*36**h+8*6**h+2*5**h+2
    return base

def midpoint_spectrum(p: int) -> list[Fraction]:
    if p<3 or not sp.isprime(p): raise ValueError('odd prime required')
    h=(p-1)//2; out=set()
    def put(x,y):
        a,b=sorted((abs(x),abs(y)))
        if a and b>a: out.add(Fraction(b,a))
    for a in range(-p,p+1):
        for b in range(-p,p+1):
            if abs(a)+abs(b)<=p and (abs(a)+abs(b))%2: put(a,b)
    for r in range(-h-1,h+2):
        for s in range(-h-1,h+2):
            if abs(r)+abs(s)<=h+1 and (abs(r)+abs(s)-h-1)%2==0: put(1+2*r,2*s)
    for u in range(-h,h+1,2):
        for v in range(-h,h+1):
            for e in (-1,1): put(1+2*u+4*v+4*e,1+2*u-4*v)
    for t in range(-h,h+1,2):
        put(7+2*t,1-2*t); put(3+2*t,5-2*t)
    for t in range(-h,h+1): put(3+2*t,1-2*t)
    return sorted(out)

def midpoint_scalar_count(p: int,rho: Fraction,lambda_zero: bool=False) -> int:
    """Independent scalar transfer DP on common C and sum of D."""
    h=(p-1)//2; a,b=rho.numerator,rho.denominator
    S=(a,-a,b,-b); table=defaultdict(Counter)
    for x1,x2,y1,y2 in product(S,repeat=4):
        if lambda_zero and x1-x2+y1-y2 != 0: continue
        c=x1+x2-y1-y2; d=x1+x2+y1+y2
        table[c][2*d]+=1
    count=0
    for c,T in table.items():
        sums={0:1}
        for _ in range(h):
            nxt=defaultdict(int)
            for u,n in sums.items():
                for v,m in T.items(): nxt[u+v]+=n*m
            sums=nxt
        count+=sum(sums.get(-c-4*x,0) for x in S)
    return count

def algebraic_lift(n: int,U: dict[int,int],denominator: int,primes: list[int]) -> str:
    """One-hot CRT lift for rho=U(zeta_n)/denominator, when rho is real>1.

    Reality and size of rho are assumptions checked independently by caller.
    Coefficients U use representatives 0<=e<n; all positions are filled.
    """
    if n<4 or n%2 or denominator<=0: raise ValueError('even n>=4, positive denominator required')
    if any(e<0 or e>=n for e in U): raise ValueError('exponent outside base block')
    if len(set(primes))!=len(primes) or any(not sp.isprime(p) or gcd(p,n)!=1 for p in primes):
        raise ValueError('fresh distinct primes required')
    t=len(primes); slots=2**t
    corner={}
    for e in range(n):
        target=(2*denominator if e==0 else 0,-2*U.get(e,0))
        units=unit_decomposition(target,slots)
        corner[e]={bits:scale((-1)**sum(bits),u) for bits,u in zip(product((0,1),repeat=t),units)}
    k=n*prod(primes); w=[None]*k
    for e in range(n):
        for indices in product(*(range(p) for p in primes)):
            E=(k//n*e+sum(k//p*u for p,u in zip(primes,indices)))%k
            bits=tuple(int(u==p-1) for u,p in zip(indices,primes))
            w[E]=corner[e][bits]
    if any(v is None for v in w): raise AssertionError('unfilled position')
    return ''.join(LETTER[v] for v in w)

def verify_algebraic_lift(word: str,n: int,U: dict[int,int],denominator: int) -> dict:
    k=len(word); phi=phi_coeff(k); d=len(phi)-1
    H=poly_mod([VECTOR[ch][0] for ch in word],phi)
    L=poly_mod([VECTOR[ch][1] for ch in word],phi)
    # Verify d*L + U(zeta_k^(k/n))*H == 0 in Q[z]/Phi_k.
    values=[0]*(k+d)
    for j,x in enumerate(L): values[j]+=denominator*x
    for e,c in U.items():
        for j,x in enumerate(H): values[e*(k//n)+j]+=c*x
    return {'k':k,'base_order':n,'degree':d,'scalar_zero':not any(poly_mod(values,phi)),
            'component_H_nonzero':any(H)}
