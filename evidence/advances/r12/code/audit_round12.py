"""Exact audit of saved endpoints and experiment scope; no solver imported.
A time limit and an infeasible numeric status are metadata, not proof premises.
"""
from pathlib import Path
from collections import Counter
import argparse, hashlib, json
from literal import ROOT, WORDS, rc, pcrs, groups, read_set, height, load, save, need
from verify_domains import universe, verify_nonstable_path


def packed_digest(selected):
    b=bytearray(len(WORDS)//8)
    for i,w in enumerate(WORDS):
        if w in selected:b[i//8] |= 1<<(7-i%8)
    return hashlib.sha256(b).hexdigest()


def model_counts(rec,base,drop):
    gs=list(groups()); ids=rec.get('free_groups')
    need(isinstance(ids,list) and ids==sorted(set(ids)) and all(type(i)is int and 0<=i<len(gs) for i in ids),'MODEL_GROUP_IDS')
    need({i for i,g in enumerate(gs) if g['stable']}<=set(ids),'ALL_STABLE_MUST_BE_FREE')
    free=set().union(*(gs[i]['words'] for i in ids)); ER=[w for w in WORDS if w<=rc(w)]
    nx=sum(w in free for w in ER)
    rows=sum(w in free or w not in base for w in ER)+sum(p[0] in free for p in pcrs())+1
    if drop:rows+=sum(w in free and not universe()[0][universe()[1][w]]['stable'] for w in ER)
    need((rec['selection_variables'],rec['variables'],rec['rows'])==(nx,nx+8192,rows),'MODEL_COUNTS')
    need(rec['stable_groups_free']==128 and rec['nonstable_groups_free']==len(ids)-128,'MODEL_FREE_COUNTS')
    return free


def verify_neutral(root):
    log=load(root/'log.json');need(len(log)==24,'NEUTRAL_COUNT')
    previous=read_set(ROOT/'inputs/S8232.words.txt',8232);base=previous;seen=set();rows=[]
    for j,r in enumerate(log):
        need(r['iteration']==j and r['seed']==12500+j,'NEUTRAL_INDEX')
        before,bi=height(load(root/f'input{j:02}.json'),previous,full_graph=False)
        after,ai=height(load(root/f'output{j:02}.json'),full_graph=True)
        free=model_counts(r,before,False)
        need(before^after<=free,'FROZEN_ENDPOINTS')
        need(r['before']==len(before) and r['after']==len(after)==r['selected_cardinality']==8232,'NEUTRAL_COST')
        need(r['L']==ai['L']==105,'NEUTRAL_L')
        need(r['selected_difference']==len(before^after),'NEUTRAL_PATCH')
        digest=packed_digest(after)
        need(r['selected_sha256_packbits']==digest,'NEUTRAL_HASH')
        need(r['distinct']==(digest not in seen),'NEUTRAL_DISTINCT')
        seen.add(digest);previous=after
        bad=[]
        for g in groups():
            c=len(g['words']&after)
            if c>(1 if g['stable'] else 2):bad.append(g['representative'])
        need(bad==['ACACACAG'],'NEUTRAL_EXCESS_LOCATION')
        rows.append({'iteration':j,'cardinality':len(after),'L':ai['L'],
                     'word_distance_from_input':len(base^after),'word_distance_from_previous':len(before^after),
                     'word_sha256':hashlib.sha256(('\n'.join(sorted(after))+'\n').encode()).hexdigest()})
    return {'status':'PASS','outputs':24,'distinct_selected_sets':len(seen),'base_seen_as_output':packed_digest(base) in seen,
            'best_cardinality':8232,'max_changed_words_from_base':max(r['word_distance_from_input'] for r in rows),
            'excess_group':'ACACACAG','rows':rows}


def verify_stable_cases():
    log=load(ROOT/'experiments/stable_defect/cases.json');gs=list(groups())
    expect=[None]
    for i,g in enumerate(gs):
        if not g['stable']:continue
        fixed=sorted(w for w in g['words'] if rc(w)==w)
        need(len(fixed)==2,'STABLE_FIXED_COUNT')
        expect.append({'group':i,'selected':[WORDS.index(w) for w in fixed]})
        pairs=sorted({tuple(sorted((w,rc(w)))) for w in g['words'] if w!=rc(w)})
        for pair in pairs:expect.append({'group':i,'selected':[WORDS.index(w) for w in pair]})
    need(len(log)==len(expect)==495,'STABLE_CASE_COUNT')
    # Enumeration order is not mathematical; compare the whole indexed universe.
    key=lambda c: None if c is None else (c['group'],tuple(c['selected']))
    need(len({key(r['case']) for r in log})==495,'STABLE_CASE_DUPLICATE')
    need({key(r['case']) for r in log}==set(map(key,expect)),'STABLE_CASE_COMPLETENESS')
    base=read_set(ROOT/'inputs/M107.words.txt',8230)
    p=load(ROOT/'certificates/nonstable_path.json');verify_nonstable_path(p,base)
    for j,r in enumerate(log):
        need(r['index']==j and r['status']=='BASE_PATH','STABLE_RUN_STATUS')
        need(all(type(i)is int and 0<=i<65536 for i in r['walk']),'STABLE_INDEX_RANGE')
        need([WORDS[i] for i in r['walk']]==p['walk'],'STABLE_SAME_PATH')
    return {'status':'PASS','cases':495,'distinct_refuting_paths':1,'SAT_branch_exercised':False,
            'complete_domain_nonfixed_excess_options':366,'fixed_pair_excess_options':128,
            'claim_boundary':'Prototype 2-SAT branch not promoted by these 495 early path rejections.'}


def audit_metadata():
    base=read_set(ROOT/'inputs/S8232.words.txt',8232);result=[]
    batches=[('near1',load(ROOT/'experiments/near1/log.json'),True),
             ('near_global',load(ROOT/'experiments/near_global/log.json'),True),
             ('adaptive',[load(ROOT/f'experiments/adaptive/solve{i}.json') for i in range(2)],True),
             ('migrate',load(ROOT/'experiments/migrate/log.json'),False)]
    for name,log,drop in batches:
        for i,r in enumerate(log):
            model_counts(r,base,drop)
            need(r['incumbent_available'] is False,'EXPECTED_NO_INCUMBENT')
            need(r['status'] in ('HighsModelStatus.kInfeasible','HighsModelStatus.kTimeLimit'),'SOLVER_STATUS')
            row={'batch':name,'index':i,'status_metadata':r['status'],'limit_seconds':r['limit'],
                 'target_cardinality':r['target'],'nonstable_groups_free':r['nonstable_groups_free'],
                 'variables':r['variables'],'rows':r['rows'],'incumbent_available':False}
            result.append(row)
            if name=='near1':
                n=r['nonstable_groups_free'];c=load(ROOT/f'certificates/regions/geometric_{n}.json')
                need({groups()[k]['representative'] for k in r['free_groups'] if not groups()[k]['stable']}==set(c['free_groups']),
                     'GEOMETRIC_MODEL_CERTIFICATE_JOIN')
    need(len(result)==10,'NUMERICAL_CALL_COUNT')
    last=load(ROOT/'experiments/global_final/solve.json')
    need(last['model_sha256']==hashlib.sha256((ROOT/'model/global_8231_L106.json.gz').read_bytes()).hexdigest(),'FINAL_SOLVE_MODEL_JOIN')
    need(last['variables']==41088 and last['rows']==69359,'FINAL_SOLVE_COUNTS')
    need(last['incumbent_available'] is False and last['status_metadata']=='HighsModelStatus.kTimeLimit','FINAL_SOLVE_STATUS')
    result.append({'batch':'global_final','status_metadata':last['status_metadata'],'limit_seconds':last['limit_seconds'],
                   'variables':last['variables'],'rows':last['rows'],'nonstable_groups_free':4051,'target_cardinality':8231,
                   'incumbent_available':False,'model_sha256':last['model_sha256']})
    return {'status':'PASS','calls':len(result),'records':result,
            'scope':'Reconstructed masks and model counts; no numeric negative is a theorem. No internal solver replay.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    result={'status':'PASS','neutral':verify_neutral(ROOT/'experiments/neutral'),
            'stable_cases':verify_stable_cases(),'negative_metadata':audit_metadata()}
    save(a.out/'experiment_audit.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
