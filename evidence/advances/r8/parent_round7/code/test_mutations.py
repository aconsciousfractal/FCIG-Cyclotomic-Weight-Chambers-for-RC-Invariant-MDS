#!/usr/bin/env python3
"""Targeted negative controls. Every mutation must raise the expected checker error.
Only the standard library and the independent checkers are used.
"""
from __future__ import annotations
import argparse,copy,json,tempfile,time
from pathlib import Path
import verify_round7 as v
import verify_model as model
import audit_search as audit
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args();start=time.time()
 base=v.read_words(ROOT/'inputs/M107.words.txt');basej=v.load(ROOT/'inputs/M107.json');stable=v.verify_base(base,basej)
 proof=v.load(ROOT/'certificates/two_sites_all_stable.json');data=v.cycle_data(proof['cycles'],base,stable)
 m=v.load(ROOT/'models/complete_T107_strengthened.json.gz');chart=v.load(ROOT/'inputs/chart_M107_T108.json.gz')
 search=v.load(ROOT/'discovery/results/three_sites.json');forced=v.load(ROOT/'discovery/results/three_sites_stable.json')
 v.verify_proof(base,stable,proof);model.check_model(base,stable,m);audit.audit_run(base,stable,search);audit.audit_run(base,stable,forced)
 results=[]
 def test(name,fn,expected):
  try:fn()
  except v.Invalid as e:
   v.require(str(e)==expected,'UNEXPECTED_MUTATION_ERROR:'+name+':'+str(e));results.append({'test':name,'status':'REJECTED','output':str(e)})
  else:raise RuntimeError('MUTATION_ACCEPTED:'+name)
 def mutate(obj,change,fn):
  x=copy.deepcopy(obj);change(x);return fn(x)
 def B(change):return mutate(basej,change,lambda x:v.verify_base(base,x))
 def P(change):return mutate(proof,change,lambda x:v.verify_proof(base,stable,x))
 def M(change):return mutate(m,change,lambda x:model.check_model(base,stable,x))
 def S(change,src=search):return mutate(src,change,lambda x:audit.audit_run(base,stable,x))
 def duplicate_word():
  with tempfile.TemporaryDirectory() as d:
   f=Path(d)/'words.txt';f.write_text('AAAAAAAA\nAAAAAAAA\n');v.read_words(f)
 test('duplicate_word_input',duplicate_word,'DUPLICATE_WORD')
 test('missing_base_word',lambda:v.verify_base(base-{'AAAAAAAA'},basej),'BASE_RC_CARDINALITY')
 test('base_digest',lambda:B(lambda x:x.__setitem__('sha256_words','0'*64)),'BASE_DIGEST')
 test('negative_base_rank',lambda:B(lambda x:x['rank'].__setitem__(0,-1)),'BASE_RANK_DOMAIN')
 test('constant_base_rank',lambda:B(lambda x:x.__setitem__('rank',[0]*16384)),'BASE_RANK_EDGE')
 test('shortened_base_path',lambda:B(lambda x:x['longest_path'].pop()),'BASE_LONGEST_PATH')
 test('broken_base_cycle',lambda:B(lambda x:x['critical_cycle'].__setitem__(0,'AAAAAAAA')),'WALK_ADJACENCY')
 test('wrong_base_metadata',lambda:B(lambda x:x.__setitem__('component_lower_bound',106)),'BASE_METADATA')
 test('wrong_target_period',lambda:P(lambda x:x.__setitem__('T',106)),'PROOF_SCOPE')
 test('wrong_reference_letter',lambda:P(lambda x:x.__setitem__('reference_letter','A')),'PROOF_SCOPE')
 test('weakened_stable_domain',lambda:P(lambda x:x.__setitem__('gamma_domain',[-1,1])),'COEFFICIENT_SCOPE')
 test('wrong_proof_base',lambda:P(lambda x:x.__setitem__('base_sha256','0'*64)),'PROOF_BASE')
 test('broken_obstruction_walk',lambda:P(lambda x:x['cycles'][0].__setitem__(0,'AAAAAAAA')),'WALK_ADJACENCY')
 test('duplicate_cover_walk',lambda:P(lambda x:x['cover_walks'].append(x['cover_walks'][0])),'COVER_INDICES')
 test('wrong_total_pair_count',lambda:P(lambda x:x.__setitem__('expected_total_pairs',33550335)),'TOTAL_PAIR_CENSUS')
 test('missing_pair_refutation',lambda:P(lambda x:x['pair_proofs'].pop()),'PAIR_COVERAGE_INCOMPLETE')
 test('duplicate_pair_refutation',lambda:P(lambda x:x['pair_proofs'].append(copy.deepcopy(x['pair_proofs'][0]))),'PAIR_COVERAGE_KEY')
 first=proof['pair_proofs'][0]
 def F(change):return mutate(first,change,lambda x:v.verify_pair_proof(x,data))
 test('negative_dual_multiplier',lambda:F(lambda x:x['multipliers'][0].__setitem__(1,-1)),'DUAL_COEFFICIENT')
 test('zero_dual_denominator',lambda:F(lambda x:x['multipliers'][0].__setitem__(2,0)),'DUAL_COEFFICIENT')
 test('out_of_range_dual_row',lambda:F(lambda x:x['multipliers'][0].__setitem__(0,100000)),'DUAL_COEFFICIENT')
 test('changed_dual_weight',lambda:F(lambda x:x['multipliers'][0].__setitem__(1,2)),'FARKAS_IDENTITY')
 test('reversed_pair',lambda:F(lambda x:x['sites'].reverse()),'PAIR_DOMAIN')
 strip=next(x for x in proof['pair_proofs'] if x['kind']=='integer_strip')
 def I(change):return mutate(strip,change,lambda x:v.verify_pair_proof(x,data))
 test('changed_strip_form',lambda:I(lambda x:x['form'][0].__setitem__(1,2)),'STRIP_DUAL_IDENTITY')
 test('changed_half_integer',lambda:I(lambda x:x['lower']['value'].__setitem__(0,0)),'STRIP_DUAL_IDENTITY')
 test('zero_strip_denominator',lambda:I(lambda x:x['lower']['value'].__setitem__(1,0)),'STRIP_VALUE')
 test('wrong_model_alphabet_size',lambda:M(lambda x:x.__setitem__('q',2)),'MODEL_SCOPE')
 test('expanded_model_domain',lambda:M(lambda x:x['variables'][0].__setitem__('upper',3)),'MODEL_VARIABLES')
 test('relaxed_model_base_row',lambda:M(lambda x:x['rows'][0].__setitem__('rhs',x['rows'][0]['rhs']-1)),'MODEL_ROWS')
 test('removed_unary_tightening',lambda:M(lambda x:x['rows'].pop(32896)),'MODEL_ROWS')
 test('altered_global_cycle_rhs',lambda:M(lambda x:x['rows'][-1].__setitem__('rhs',x['rows'][-1]['rhs']+1)),'MODEL_ROWS')
 test('wrong_control_period',lambda:mutate(chart,lambda x:x.__setitem__('T',107),lambda x:model.verify_control(base,stable,x,m['strengthening']['cycles'])),'CONTROL_BASE_PERIOD')
 test('wrong_initial_cut_count',lambda:S(lambda x:x.__setitem__('initial_cycles',192)),'INITIAL_CUT_CORPUS')
 test('wrong_search_period',lambda:S(lambda x:x.__setitem__('T',108)),'SEARCH_SCOPE')
 test('unbounded_search_metadata',lambda:S(lambda x:x.__setitem__('z_bounds',[-3,3])),'SEARCH_DOMAIN')
 test('wrong_actual_site_count',lambda:S(lambda x:x['audit'][0].__setitem__('actual_sites',2)),'ACTUAL_SUPPORT_METADATA')
 test('wrong_actual_stable_count',lambda:S(lambda x:x['audit'][0].__setitem__('stable_changes',1)),'ACTUAL_SUPPORT_METADATA')
 test('noninteger_search_coefficient',lambda:S(lambda x:x['audit'][0]['operations'][0].__setitem__(1,0.5)),'OPERATION_INDEX')
 test('wrong_separation_index',lambda:S(lambda x:x['audit'][0].__setitem__('cycle_index',0)),'SEPARATOR_INDEX')
 test('wrong_iteration_number',lambda:S(lambda x:x['audit'][0].__setitem__('iteration',1)),'SEARCH_ITERATION')
 def remove_stable(x):
  x['audit'][0]['operations']=[t for t in x['audit'][0]['operations'] if t[0]<285];x['audit'][0]['stable_changes']=0
 test('activation_without_actual_stable_flip',lambda:S(remove_stable,forced),'ACTUAL_STABLE_CHANGE_REQUIRED')
 out={'status':'PASS','pristine_controls':'PASS','mutations_rejected':len(results),'tests':results,'seconds':time.time()-start}
 if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
