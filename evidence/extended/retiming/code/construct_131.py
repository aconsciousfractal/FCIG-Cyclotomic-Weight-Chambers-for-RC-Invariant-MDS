"""Construct the 8230-word RC MDS with lower residual length 131.
Only the frozen base input and four explicit patches are used. Each endpoint
is checked for literal RC, one word per PCR, and residual acyclicity.
"""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json
import retime as rt
PATCHES=((3494,(3,)),(3870,(2,3)),(5718,(1,)),(12158,(3,)))

def construct(root,output):
    q,k,n=4,8,16384
    m=rt.read_bits(root/'inputs/k8_c094.bits',65536)
    history=[rt.dag(m,q,k)['L']]
    for u,letters in PATCHES:
        delta=Counter()
        for a in letters:
            for e,sgn in ((a*n+u,-1),(q*u+a,1)):
                delta[e]+=sgn;delta[rt.theta(e,k,q)]+=sgn
        for e,sgn in delta.items():m[e]+=sgn
        rt.require(all(x in (0,1) for x in m),'PATCH_BOOLEAN')
        rt.require(all(m[e]==m[rt.theta(e,k,q)] for e in range(len(m))),'PATCH_RC')
        rt.require(rt.pcr_check(m,q,k)==8230,'PATCH_SIZE')
        history.append(rt.dag(m,q,k)['L'])
    rt.require(history==[141,135,135,133,131],'PATCH_LENGTHS')
    rt.write_bits(output,m)
    rt.require(output.read_bytes()==(root/'certificates/escape_131.bits').read_bytes(),'FINAL_BITS')
    return {'status':'PASS','lengths':history,'size':sum(m),'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    print(json.dumps(construct(Path(__file__).resolve().parents[1],a.out),sort_keys=True,indent=2))
