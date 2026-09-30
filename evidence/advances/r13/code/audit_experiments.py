"""Reconstruct numerical outputs and signed cycle rejections with exact integers.
No solver result is used as a proof of infeasibility or optimality.
"""
from literal import *
from verify_certificates import check_endpoint
from collections import Counter,deque
import argparse,struct
REPS=[v for v in V if v<rc(v)];NZ=3*len(REPS)

def signed_generators(m):
    G=[]
    for u in REPS:
        for a in 'ACG':
            d=Counter()
            for w,c in [(u+a,1),(rc(u+a),1),(a+u,-1),(rc(a+u),-1)]:d[w]+=c
            G.append({w:c for w,c in d.items() if c})
    for g in groups():
        if not g['stable']:continue
        selected=sorted(g['words']&m);fixed=sorted(w for w in g['words'] if w==rc(w))
        need(len(selected)==1 and len(fixed)==2 and selected[0] in fixed,'SIGNED_STABLE_SOURCE')
        other=next(w for w in fixed if w!=selected[0]);G.append({other:1,selected[0]:-1})
    need(len(G)==24704,'SIGNED_GENERATOR_COUNT');return G

def fixed_tree():
    seen={'AAAAAAA'};q=deque(['AAAAAAA']);tree=[]
    while q:
        u=q.popleft()
        for a in ABC:
            w=u+a;v=w[1:]
            if v not in seen:seen.add(v);q.append(v);tree.append(w)
    need(seen==set(V) and len(tree)==len(V)-1,'NORMALIZATION_TREE');return tree
TREE=fixed_tree()
CANONICAL_BYTES={}
def normalize(x):
    f={'AAAAAAA':0}
    for w in TREE:f[w[1:]]=f[w[:-1]]-x[w]
    y=[x[w]-f[w[:-1]]+f[w[1:]] for w in W]
    data=struct.pack('<'+str(len(W))+'i',*y)
    key=hashlib.sha256(data).hexdigest()
    if key in CANONICAL_BYTES:need(CANONICAL_BYTES[key]==data,'CANONICAL_HASH_COLLISION')
    else:CANONICAL_BYTES[key]=data
    return key

def check_signed(rec):
    m=read_set(safe_path(rec['source']));need(digest(m)==rec['source_sha256'],'SIGNED_SOURCE_HASH');G=signed_generators(m)
    need(rec['force_gamma'] is False,'SIGNED_FORCE_SCOPE')
    cuts=list(rec['initial_cuts']);need(cuts,'SIGNED_INITIAL_CUTS')
    for p in cuts:walk(p,closed=True)
    seen=set();classes=[];active=Counter();flip=Counter();out=[]
    need(rec['status']=='ITERATION_LIMIT' and len(rec['iterations'])==100,'SIGNED_RUN_STATUS')
    for i,r in enumerate(rec['iterations']):
        need(r['iteration']==i and r['incumbent'] is True,'SIGNED_ITERATION')
        z=r.get('coordinates');need(type(z) is dict,'SIGNED_COORDINATES');co={}
        for k,v in z.items():
            need(str(int(k))==k and 0<=int(k)<len(G) and type(v) is int and v!=0,'SIGNED_INDEX_TYPE')
            j=int(k);need((-2<=v<=2) if j<NZ else v==1,'SIGNED_DOMAIN');co[j]=v
        x={w:int(w in m) for w in W}
        for j,v in co.items():
            for w,c in G[j].items():x[w]+=v*c
        need(all(x[w]==x[rc(w)] for w in W),'SIGNED_RC')
        need(all(sum(x[w] for w in p)==1 for p in pcrs()),'SIGNED_PCR_SUM')
        raw=hashlib.sha256(struct.pack('<'+str(len(W))+'i',*(x[w] for w in W))).hexdigest()
        need(raw==r.get('cochain_sha256_i32'),'SIGNED_COCHAIN_HASH');need(r.get('distinct')==(raw not in seen),'SIGNED_DISTINCT_FLAG');seen.add(raw)
        k=len({j//3 for j in co if j<NZ});st=sum(v for j,v in co.items() if j>=NZ)
        need(k==r.get('active_sites') and st==r.get('stable_changes'),'SIGNED_ACTIVE_CHANGES');active[k]+=1;flip[st]+=1
        for p in cuts:need(107*sum(x[w] for w in p)>=len(p),'SIGNED_PREVIOUS_CUT')
        p=r.get('obstruction');walk(p,closed=True)
        val=len(p)-107*sum(x[w] for w in p)
        need(r.get('oracle_status')==1 and val>0 and val==r.get('violation'),'SIGNED_REJECTION')
        need(Counter(p) not in [Counter(q) for q in cuts],'SIGNED_DUPLICATE_OBSTRUCTION');cuts.append(p)
        key=normalize(x);classes.append(key);out.append({'iteration':i,'cochain_hash':raw,'gradient_class_fingerprint':key,'violation':val})
    return {'status':'PASS','candidates':len(out),'distinct_raw_cochains':len(seen),'distinct_gradient_classes':len(set(classes)),'active_sites_histogram':dict(active),'stable_changes_histogram':dict(flip),'final_constraints':len(cuts),'rows':out}

def stats(s):
    d,b,p=graph(s);gs=groups();cen=next(g for g in gs if g['representative']=='ACACACAG');B=load(ROOT/'upstream/round12/certificates/regions/adaptive_1.json')['boundary_words']
    return {'cardinality':len(s),'L':len(p),'sha256':digest(s),'center_count':len(set(cen['pcr'])&s),'B21_hits':[w for w in B if w in s],'overselected_groups':[{'group':g['representative'],'stable':g['stable'],'selected_count_full_group':len(g['words']&s)} for g in gs if len(g['words']&s)>(1 if g['stable'] else 2)]}

def check_milp(name,source):
    logs=load(ROOT/'experiments'/name/'log.json');m=read_set(safe_path(source));gs={g['representative']:g for g in groups()};out=[]
    need(len(logs)==4,'MILP_CALL_COUNT')
    for i,r in enumerate(logs):
        st=stats(m);need(r['iteration']==i and r['input']==st,'MILP_INPUT_METADATA')
        free=r['free_groups'];need(type(free) is list and free==sorted(set(free)) and all(g in gs for g in free),'MILP_FREE_GROUPS')
        need(all(g in free for g in gs if gs[g]['stable']),'MILP_STABLE_AVAILABLE')
        need(r['target_L']==106 and r['center_single'] and r['border_hit'] and r['lifted'],'MILP_GATES')
        if r['incumbent']:
            p=ROOT/'experiments'/name/f'output{i:03}.words.txt';n=read_set(p);ns=stats(n)
            need(r['output']==ns and ns['cardinality']<=r['cap'] and ns['L']<=106,'MILP_OUTPUT_METADATA')
            need(ns['center_count']==1 and ns['B21_hits'],'MILP_TRUE_ESCAPE')
            allowed=frozenset().union(*(gs[g]['words'] for g in free));need(m^n<=allowed,'MILP_FROZEN_INDICATORS')
            need(r['delta']==len(m^n),'MILP_DELTA')
            if len(n)<=8232:m=n
        else:need('output' not in r,'MILP_NO_INCUMBENT_NOT_A_WITNESS')
        out.append({'iteration':i,'status':r['solver_status'],'incumbent':r['incumbent'],'cap':r['cap'],'nonstable_groups':sum(not gs[g]['stable'] for g in free),'word_cost':r.get('output',{}).get('cardinality')})
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    report={'status':'RUNNING','signed':[],'milp':{}}
    allraw=set();allclasses=set()
    for name in ['signed_alt2_seed13220','signed_alt3_seed13230']:
        r=check_signed(load(ROOT/'experiments'/name/'log.json'));r['name']=name;report['signed'].append(r)
        allraw.update(x['cochain_hash'] for x in r['rows']);allclasses.update(x['gradient_class_fingerprint'] for x in r['rows'])
    report['signed_total_raw']=len(allraw);report['signed_total_gradient_classes']=len(allclasses);report['equal_fingerprints_checked_by_exact_bytes']=True
    jobs=[('forced_migration','inputs/S8232.words.txt'),('repair_migrated','experiments/forced_migration/output000.words.txt'),('alternative_2_cap8231','experiments/reuse_107/extension_002_000.words.txt'),('alternative_3_cap8231','experiments/reuse_107/extension_003_000.words.txt')]
    for name,src in jobs:report['milp'][name]=check_milp(name,src)
    old=load(ROOT/'experiments/coalesce_budget_original/log.json');need(len(old)==24 and all(r['status']==2 and r['operations']==8000001 and 'cycle' not in r for r in old),'BUDGET_FAILURE_PRESERVED')
    report['failed_budget_cases']=len(old);report['status']='PASS';save(a.out/'experiments.json',report)
    brief={k:v for k,v in report.items() if k!='signed'};brief['signed']=[{k:v for k,v in r.items() if k!='rows'} for r in report['signed']];print(json.dumps(brief,indent=2))
if __name__=='__main__':main()
