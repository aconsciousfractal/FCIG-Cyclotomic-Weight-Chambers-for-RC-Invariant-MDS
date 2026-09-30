"""Deliberate corruptions of certificate data. No discovery or solver imports."""
from pathlib import Path
from copy import deepcopy
import argparse,gzip,hashlib,json
from literal import ROOT, read_set, load, save, height, need, rc
from verify_domains import verify_trace,verify_nonstable_path,universe
from verify_model import verify_model

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    base=read_set(ROOT/'inputs/S8232.words.txt',8232);dig=hashlib.sha256((ROOT/'inputs/S8232.words.txt').read_bytes()).hexdigest()
    small=load(ROOT/'certificates/regions/adaptive_1.json');large=load(ROOT/'certificates/regions/adaptive_139.json')
    partial=json.loads(gzip.decompress((ROOT/'certificates/adaptive_686_pruning.json.gz').read_bytes()))
    path=load(ROOT/'certificates/nonstable_path.json');m=read_set(ROOT/'inputs/M107.words.txt',8230)
    rows=[]
    def reject(name,obj,mutate,check):
        c=deepcopy(obj);mutate(c)
        try:check(c)
        except ValueError as e:rows.append({'name':name,'status':'REJECTED','reason':str(e)})
        else:raise ValueError('MUTATION_ACCEPTED: '+name)
    tr=lambda c:verify_trace(c,base,dig)
    rej=lambda name,mut:reject(name,small,mut,tr)
    rej('schema',lambda c:c.update(schema='other'))
    rej('length_105',lambda c:c.update(target_L=105))
    rej('length_float',lambda c:c.update(target_L=106.0))
    rej('base_digest',lambda c:c.update(base_sha256='0'*64))
    rej('duplicate_free_group',lambda c:c['free_groups'].append(c['free_groups'][0]))
    st=next(g for g in universe()[0] if universe()[0][g]['stable'])
    rej('stable_group_claimed_nonstable',lambda c:c.update(free_groups=sorted(c['free_groups']+[st])))
    rej('empty_steps',lambda c:c.update(steps=[]))
    rej('missing_final_elimination',lambda c:c['steps'].pop())
    rej('duplicate_elimination',lambda c:c['steps'].insert(1,deepcopy(c['steps'][0])))
    rej('bad_test_option',lambda c:c['steps'][0]['option'].__setitem__(0,'AAAAAAAA'))
    rej('bad_step_group',lambda c:c['steps'][0].update(group='AAAAAAAA'))
    rej('unknown_kind',lambda c:c['steps'][0].update(kind='maybe'))
    rej('broken_continuity',lambda c:c['steps'][0]['walk'].__setitem__(0,'AAAAAAAA'))
    rej('open_cycle',lambda c:c['steps'][0]['walk'].pop())
    rej('cycle_called_long_path',lambda c:c['steps'][0].update(kind='path'))
    w=st;cy=[]
    while w not in cy:cy.append(w);w=w[1:]+w[0]
    rej('stable_PCR_cycle_cannot_be_mandatory',lambda c:c['steps'][0].update(walk=cy))
    rej('root_not_empty',lambda c:c.update(empty_group='AAAAAAAA'))
    rej('wrong_end_status',lambda c:c.update(result='SAT'))
    rej('boundary_word_missing',lambda c:c['boundary_words'].pop())
    rej('boundary_selected_word',lambda c:c['boundary_words'].append('AAAAAAAA'))
    rej('boundary_group_missing',lambda c:c['boundary_groups'].pop())
    rej('single_assumptions_missing',lambda c:c.update(single_assumptions=[]))
    reject('path_length_106',large,lambda c:next(s for s in c['steps'] if s['kind']=='path')['walk'].pop(),tr)
    # Locate and reverse a genuine dependency, instead of assuming adjacent steps depend.
    _,own,_,selecting=universe();done={};pair=None
    for i,s in enumerate(large['steps']):
        for w in s['walk']:
            if own[w] in large['free_groups'] and own[w]!=s['group']:
                j=done.get((own[w],selecting[w]))
                if j is not None:pair=(j,i);break
        if pair:break
        done[(s['group'],tuple(s['option']))]=i
    need(pair is not None,'MUTATION_DEPENDENCY_LOCATED')
    reject('future_dependency',large,lambda c:c['steps'].insert(0,c['steps'].pop(pair[1])),tr)
    ptr=lambda c:verify_trace(c,base,dig,partial=True)
    first=partial['free_groups'][0]
    reject('partial_missing_option',partial,lambda c:c['remaining_options'][first].pop(),ptr)
    reject('partial_duplicate_option',partial,lambda c:c['remaining_options'][first].append(c['remaining_options'][first][0]),ptr)
    reject('partial_missing_group',partial,lambda c:c['remaining_options'].pop(first),ptr)
    reject('partial_false_unsat_status',partial,lambda c:c.update(result='EMPTY_DOMAIN'),ptr)
    reject('partial_omitted_step',partial,lambda c:c['steps'].pop(),ptr)
    reject('nonstable_path_short',path,lambda c:c['walk'].pop(),lambda c:verify_nonstable_path(c,m))
    reject('nonstable_path_invalid_word',path,lambda c:c['walk'].__setitem__(0,'AAAAAAAN'),lambda c:verify_nonstable_path(c,m))
    h=load(ROOT/'inputs/S8232.height.json')
    reject('height_cost',h,lambda c:c.update(cardinality=8231),lambda c:height(c))
    reject('height_bound',h,lambda c:c['heights'].__setitem__(0,54),lambda c:height(c))
    reject('height_nonintegral',h,lambda c:c['heights'].__setitem__(0,0.5),lambda c:height(c))
    mod=json.loads(gzip.decompress((ROOT/'model/global_8231_L106.json.gz').read_bytes()))
    reject('model_restricted',mod,lambda c:c.update(all_groups_free=False),verify_model)
    reject('model_wrong_target',mod,lambda c:c.update(length_bound=107),verify_model)
    reject('model_relaxed_indicator',mod,lambda c:c['variables'][0].update(type='real'),verify_model)
    reject('model_missing_row',mod,lambda c:c['rows'].pop(),verify_model)
    reject('model_false_escape_sign',mod,lambda c:next(r for r in c['rows'] if r['name'].startswith('escape_'))['terms'][0].__setitem__(1,-1),verify_model)
    # Positive order-invariance check for finite domains, which are mathematical sets.
    c=deepcopy(partial)
    for v in c['remaining_options'].values():v.reverse()
    ptr(c)
    result={'status':'PASS','rejected':len(rows),'tests':rows,'positive_control':'Remaining-domain option order is immaterial.'}
    save(a.out/'mutation_tests.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
