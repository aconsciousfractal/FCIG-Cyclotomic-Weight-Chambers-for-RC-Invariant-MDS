"""Read-only selected-case verifier. Python standard library only.
Usage: python -B -X utf8 verify_case.py --capsule .
The kernel is a declared shared dependency, not an independent second source.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse, hashlib, importlib.util, json

def verify(root, certificate=None):
    root=Path(root).resolve()
    kernel=root/'code/exact_kernel.py'
    expected='0a221c5088a08bf3f762ef418ae1442a85d1839100dc044466bc93ca497c72d0'
    if hashlib.sha256(kernel.read_bytes()).hexdigest()!=expected:
        raise ValueError('kernel identity mismatch')
    spec=importlib.util.spec_from_file_location('selected_case_kernel',kernel)
    K=importlib.util.module_from_spec(spec);spec.loader.exec_module(K)
    n=K.Net();need=K.need
    need((n.n,n.m,len(n.orbits),len(n.stable_p),len(n.S),len(n.groups),len(n.options),len(n.greps))==(16384,65536,8230,128,988,4051,32274,622),'carrier census')
    c=K.load(root/'results/joint_core106.json.gz') if certificate is None else certificate
    need((c['q'],c['k'],c['T'],c['L'])==(4,8,107,106),'certificate parameters')
    need(n.wordsfile(root/'results/joint_core106.words.txt')==set(c['selected']),'word-list identity')
    _,length=n.core(c['selected'],c['rank'],c['path'],106)
    need(len(c['selected'])==8102,'selected count')
    walks=K.pool(n,root);rows,mapping=K.compile_rows(n,walks)
    old=K.load(root/'inputs/compact_module.json.gz')
    need(rows==old['rows'] and mapping==old['source_walk_to_row'] and n.greps==old['stable_representatives'],'literal pool reconstruction')
    need(len(walks)==1804 and len(rows)==1031,'pool/module counts')
    need(K.lift_check(n,rows,c['selected'],c['gamma'])==128,'signed lift mass')
    selected=set(c['selected']);values=[];seen=set();ext=Counter()
    need(isinstance(c['families'],list),'return families')
    for f in c['families']:
        need(type(f['stable_index']) is int and 0<=f['stable_index']<len(n.stable_p),'PCR index')
        p=n.stable_p[f['stable_index']]
        need(f['stable_index'] not in seen,'distinct PCR');seen.add(f['stable_index'])
        need(f['representative']==n.word(min(p)),'PCR representative')
        need(sorted(t['stable_edge'] for t in f['walks'])==p,'complete return family')
        a=Counter();rhs=0
        for t in f['walks']:
            w=[t['stable_edge']]+t['return'];n.walk(w,True)
            need(all(e not in n.S for e in t['return']),'E0 return')
            v=(len(w)+106)//107;need(v==t['rhs'],'ceiling');rhs+=v
            a.update(t['return'])
        balance=[0]*n.n
        for e,v in a.items():balance[e%n.n]+=v;balance[e//n.q]-=v
        need(not any(balance),'external circulation')
        value=sum(v for e,v in a.items() if e in selected)
        need(rhs-value==2,'stable lower two');values.append([f['representative'],rhs,value])
        ext.update(a)
    need(len(seen)==5 and values==[['ATCCCGGG',10,8],['ATCGCGCG',9,7],['CCCGGGCG',9,7],['CCGCGCGG',9,7],['CCGCGGCG',9,7]],'family values')
    need(sum(v for e,v in ext.items() if e in selected)==36,'weighted cut value')
    need(sum(len(f['walks']) for f in c['families'])==40,'return count')
    return {'status':'PASS','selected':8102,'core_length':length,'old_pool_walks':len(walks),'old_module_rows':len(rows),'stable_coordinates':len(n.greps),'signed_lift_mass':128,'return_walks':40,'stable_lower':133,'full_class_lower':8235,'family_values':values,'scope':'One exact core and its E0 real-gradient class; no endpoint, attained minimum, or global lower.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--capsule',type=Path,default=Path(__file__).resolve().parent)
    print(json.dumps(verify(ap.parse_args().capsule),sort_keys=True))
if __name__=='__main__':main()
