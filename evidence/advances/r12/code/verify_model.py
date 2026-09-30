"""Independent literal-word reconstruction of every model column and row."""
from pathlib import Path
from collections import Counter
import argparse,gzip,json
from literal import ROOT, WORDS, VERTICES, rc, ER, VR, pcrs, read_set, graph, need, save, load
from verify_domains import universe, verify_trace, verify_nonstable_path
import hashlib

def rebuild(B,cuts):
    gs,own,opts,_=universe()
    names=['x_'+w for w in ER]+['r_'+v for v in VR];idx={v:i for i,v in enumerate(names)}
    variables=[{'name':'x_'+w,'type':'binary','lb':0,'ub':1,'cost':1 if rc(w)==w else 2} for w in ER]
    variables +=[{'name':'r_'+v,'type':'real','lb':0,'ub':B,'cost':0} for v in VR]
    def col(w):return idx['x_'+min(w,rc(w))]
    period={w:len(p) for p in pcrs() for w in p}
    def row(name,d,low=None,high=None):
        return {'name':name,'terms':[[idx[k],v] for k,v in sorted(d.items(),key=lambda kv:idx[kv[0]]) if v], 'lb':low,'ub':high}
    rows=[]
    for w in ER:
        d={};constant=0
        for v,coefficient in ((w[1:],1),(w[:-1],-1)):
            key='r_'+min(v,rc(v))
            if v>rc(v):constant+=B*coefficient;coefficient=-coefficient
            d[key]=d.get(key,0)+coefficient
        q=d.copy();q['x_'+w]=B+1
        rows.append(row('edge_'+w,q,1-constant))
        if not gs[own[w]]['stable']:
            q=d.copy();q['x_'+w]=B+period[w]-1
            rows.append(row('drop_'+w,q,high=B-constant))
    for g in gs.values():
        rep=g['representative']; p=next(o for o in pcrs() if o[0]==rep)
        d=Counter('x_'+min(w,rc(w)) for w in p)
        rows.append(row('group_'+rep,d,1,2 if g['stable'] else 1))
    rows.append(row('cardinality',{v['name']:v['cost'] for v in variables if v['cost']},8230,8231))
    if cuts:
        need(B==106,'MODEL_CUT_BOUND')
        base=read_set(ROOT/'inputs/S8232.words.txt',8232)
        dig=hashlib.sha256((ROOT/'inputs/S8232.words.txt').read_bytes()).hexdigest()
        for f in sorted((ROOT/'certificates/regions').glob('*.json')):
            c=load(f);verify_trace(c,base,dig)
            rows.append(row('escape_'+f.stem,{'x_'+w:1 for w in c['boundary_words']},1))
        pc=load(ROOT/'certificates/nonstable_path.json')
        verify_nonstable_path(pc,read_set(ROOT/'inputs/M107.words.txt',8230))
        rows.append(row('path_nonstable_M107',Counter('x_'+min(w,rc(w)) for w in pc['walk']),1))
    return variables,rows

def verify_model(m):
    need(m.get('schema')=='R12-GLOBAL-NEAR-1','MODEL_SCHEMA')
    need(m.get('all_groups_free') is True,'MODEL_NOT_GLOBAL')
    need(m.get('length_bound')==106 and m.get('cardinality_upper')==8231,'MODEL_TARGET')
    need(m.get('escape_cuts') is True,'MODEL_CUTS')
    v,r=rebuild(106,True)
    need(m.get('variables')==v,'MODEL_VARIABLES')
    need(m.get('rows')==r,'MODEL_ROWS')
    need((len(v),len(r))==(41088,69359),'MODEL_COUNTS')
    return {'status':'PASS','variables':len(v),'Boolean_RC_indicators':len(ER),
            'real_symmetric_ranks':len(VR),'rows':len(r),'escape_clauses':8,'nonstable_path_cuts':1}

def positive_control():
    # This checks the uncut core at B107, NOT the B106 target. The source M107
    # cannot pass the target's new path/escape requirements by being renamed.
    v,rows=rebuild(107,False);s=read_set(ROOT/'inputs/M107.words.txt',8230);d,_,p=graph(s)
    z=[2*int(w in s) for w in ER]+[107+d[u]-d[rc(u)] for u in VR]
    for vv,zz in zip(v,z):need(2*vv['lb']<=zz<=2*vv['ub'],'POSITIVE_CONTROL_DOMAIN')
    for r in rows:
        val=sum(a*z[j] for j,a in r['terms'])
        need(r['lb'] is None or val>=2*r['lb'],'POSITIVE_CONTROL_LOWER')
        need(r['ub'] is None or val<=2*r['ub'],'POSITIVE_CONTROL_UPPER')
    return {'status':'PASS','cardinality':len(s),'L':len(p),'test_bound':107,'escape_cuts_included':False,'target_witness':False}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    m=json.loads(gzip.decompress((ROOT/'model/global_8231_L106.json.gz').read_bytes()))
    result=verify_model(m);result['positive_control']=positive_control()
    save(a.out/'model_verification.json',result);print(json.dumps(result,indent=2))
