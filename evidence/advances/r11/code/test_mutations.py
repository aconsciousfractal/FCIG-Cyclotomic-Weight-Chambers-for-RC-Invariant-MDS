#!/usr/bin/env python3
"""Fault injection against finite certificates. Each alteration must be rejected."""
from exact import *
from audit_experiments import fusion_row
from copy import deepcopy
import argparse
import tempfile


def main(out):
    s=read_set(ROOT/'certificates/S8232.words.txt',8232)
    old=read_set(ROOT/'inputs/S8234.words.txt',8234)
    hc=load(ROOT/'certificates/S8232.height.json')
    rp=load(ROOT/'certificates/S8232.rank_path.json')
    ir=load(ROOT/'certificates/irredundancy.json')
    block=load(ROOT/'experiments/pcr_phase1/log.json')[2]
    rows=[]
    def bad(name, fn):
        try:fn()
        except (ValueError,KeyError,TypeError,IndexError) as exc:
            rows.append({'name':name,'status':'REJECTED','output':f'{type(exc).__name__}: {exc}'});return
        raise RuntimeError('ALTERATION_ACCEPTED: '+name)
    def obj(name,original,change,check):
        c=deepcopy(original);change(c);need(c!=original,'NO_OP_MUTATION: '+name)
        bad(name,lambda:check(c))
    with tempfile.TemporaryDirectory(prefix='round11_mutations_') as td:
        p=Path(td)/'words.txt';words=sorted(s)
        variants=[('duplicate_word',words+[words[0]],None),
                  ('missing_cardinality',words[1:],8232),
                  ('invalid_alphabet',['N'+words[0][1:]]+words[1:],None),
                  ('invalid_length',[words[0][:-1]]+words[1:],None)]
        for name,ws,card in variants:
            p.write_text('\n'.join(ws)+'\n');bad(name,lambda card=card:read_set(p,card))
        w=next(w for w in words if rc(w)!=w)
        p.write_text('\n'.join(x for x in words if x!=w)+'\n')
        bad('broken_RC',lambda:read_set(p))
    H=lambda c:height(c,s)
    obj('height_noninteger',hc,lambda c:c['heights'].__setitem__(0,0.5),H)
    obj('height_bool_not_integer',hc,lambda c:c['heights'].__setitem__(0,True),H)
    obj('height_outside_domain',hc,lambda c:c['heights'].__setitem__(0,54),H)
    obj('height_truncated',hc,lambda c:c['heights'].pop(),H)
    i=next(i for i,x in enumerate(hc['heights']) if x<53)
    obj('height_antisymmetry',hc,lambda c:c['heights'].__setitem__(i,c['heights'][i]+1),H)
    obj('height_wrong_target',hc,lambda c:c.__setitem__('target_L',105),H)
    obj('height_wrong_half_bound',hc,lambda c:c.__setitem__('half_height',52),H)
    obj('height_wrong_cardinality',hc,lambda c:c.__setitem__('cardinality',8231),H)
    bad('height_wrong_join',lambda:height(hc,s-{next(iter(s))}))
    R=lambda c:rank_path(c,s)
    obj('rank_wrong_hash',rp,lambda c:c.__setitem__('selected_sha256','0'*64),R)
    obj('rank_outside_bound',rp,lambda c:c['rank'].__setitem__(0,106),R)
    obj('rank_noninteger',rp,lambda c:c['rank'].__setitem__(0,0.5),R)
    obj('rank_truncated',rp,lambda c:c['rank'].pop(),R)
    obj('rank_edge_violation',rp,lambda c:c.__setitem__('rank',[0]*16384),R)
    obj('rank_wrong_length',rp,lambda c:c.__setitem__('L',106),R)
    obj('path_truncated',rp,lambda c:c['path'].pop(),R)
    obj('path_bad_label',rp,lambda c:c['path'].__setitem__(0,'NNNNNNNN'),R)
    obj('path_broken_continuity',rp,lambda c:c['path'].__setitem__(1,c['path'][3]),R)
    obj('avoiding_word_wrong_join',rp,lambda c:c.__setitem__('avoiding_word',c['avoiding_word']+'A'),R)
    I=lambda c:irredundancy(c,s)
    obj('irredundancy_missing_PCR',ir,lambda c:c['overselected_PCRs'].pop(),I)
    obj('irredundancy_wrong_PCR_selection',ir,lambda c:c['overselected_PCRs'][0]['selected'].pop(),I)
    obj('irredundancy_missing_cycle',ir,lambda c:c['cycles'].pop(),I)
    obj('irredundancy_duplicate_cycle',ir,lambda c:c['cycles'].append(deepcopy(c['cycles'][0])),I)
    obj('irredundancy_invalid_edge',ir,lambda c:c['cycles'][0]['cycle'].__setitem__(0,'NNNNNNNN'),I)
    obj('irredundancy_cycle_continuity',ir,lambda c:c['cycles'][0]['cycle'].__setitem__(1,'AAAAAAAA'),I)
    obj('irredundancy_repeated_cycle',ir,lambda c:c['cycles'][0].__setitem__('cycle',c['cycles'][0]['cycle']*2),I)
    def rotate_labels(c):
        c['cycles'][0]['selected_edge'],c['cycles'][1]['selected_edge']=c['cycles'][1]['selected_edge'],c['cycles'][0]['selected_edge']
    obj('irredundancy_wrong_unique_cut',ir,rotate_labels,I)
    B=lambda c:block_summary(c,old,s)
    stable=next(i for i in block['free_groups'] if groups()[i]['stable'])
    obj('block_stable_group_frozen',block,lambda c:c['free_groups'].remove(stable),B)
    changed=next(i for i in block['free_groups'] if not groups()[i]['stable'] and
                 (groups()[i]['words']&old)!=(groups()[i]['words']&s))
    obj('block_changed_group_frozen',block,lambda c:c['free_groups'].remove(changed),B)
    obj('block_duplicate_group',block,lambda c:c['free_groups'].append(c['free_groups'][-1]),B)
    obj('block_row_count',block,lambda c:c.__setitem__('rows',c['rows']+1),B)
    obj('block_variable_count',block,lambda c:c.__setitem__('variables',c['variables']-1),B)
    obj('block_false_improvement',block,lambda c:c.__setitem__('after',8230),B)
    obj('block_false_length',block,lambda c:c.__setitem__('L',104),B)
    f=ROOT/'experiments/fusion';fr=load(f/'log.json')[0]
    fi=load(f/'input_000.json');fg=load(f/'proposal_000.json');fo=load(f/'output_000.json')
    F=lambda c:fusion_row(c,fi,fg,fo)
    obj('fusion_false_hash',fr,lambda c:c.__setitem__('solution_sha256','f'*64),F)
    obj('fusion_false_change_count',fr,lambda c:c.__setitem__('changed_pairs',1),F)
    obj('fusion_false_cost',fr,lambda c:c.__setitem__('after',8232),F)
    save(out/'mutation_tests.json',{'status':'PASS','rejected':len(rows),'tests':rows})
    print(f'PASS: {len(rows)} genuine alterations rejected.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);main(a.out)
