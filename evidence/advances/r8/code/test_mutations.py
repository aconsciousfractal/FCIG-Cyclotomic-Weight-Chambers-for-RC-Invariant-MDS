#!/usr/bin/env python3
"""Adversarial finite tests. Every altered certificate must be rejected.
These tests are assurance, not a substitute for the mathematical arguments.
"""
from pathlib import Path
import argparse,copy,gzip,json,tempfile,time
import verify_round8 as v
import verify_model8 as vm
import audit_round8 as audit
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();start=time.time();base=v.p.read_words(ROOT/'inputs/M107.words.txt');stable=v.p.pcr_data(base)
 source=v.load(ROOT/'parent_round7/certificates/two_sites_all_stable.json');kernel=v.load(ROOT/'certificates/lifted_support_kernel.json');allp=v.load(ROOT/'certificates/triple_refutations.json.gz');data=v.rows_from_words(allp['cycles'],base,stable);results=[]
 def reject(name,call):
  try:call()
  except (v.Invalid,v.p.Invalid,ValueError,KeyError,IndexError,TypeError,ZeroDivisionError) as e:results.append({'test':name,'status':'REJECTED','output':f'{type(e).__name__}: {e}'});return
  raise AssertionError('MUTATION_ACCEPTED: '+name)
 def change(obj,edit):q=copy.deepcopy(obj);edit(q);return q
 def kt(name,edit):reject(name,lambda:v.verify_kernel(base,stable,change(kernel,edit),source))
 kt('kernel_wrong_period',lambda q:q.__setitem__('T',106))
 kt('kernel_wrong_base',lambda q:q.__setitem__('base_sha256','0'*64))
 kt('kernel_missing_parent_lift',lambda q:q['lifted_supports'].pop())
 kt('kernel_wrong_parent_index',lambda q:q['lifted_supports'][0].__setitem__('pair_proof',1))
 kt('kernel_deleted_bad_site',lambda q:q['lifted_supports'][0]['support'].pop())
 kt('kernel_missing_hyperedge',lambda q:q['minimal_edges'].pop())
 kt('kernel_missing_triple',lambda q:q['candidate_triples'].pop())
 kt('kernel_forged_small_support',lambda q:q['supports_below_three'].append(['AAAAAAA']))
 kt('kernel_wrong_total_population',lambda q:q.__setitem__('total_triples',1))
 examples={kind:next(p for p in allp['proofs'] if p['kind']==kind) for kind in allp['proof_kinds']}
 def pt(name,kind,edit):reject(name,lambda:v.verify_case(change(examples[kind],edit),data))
 F='gcd_farkas';S='integer_strip';D='dominated_gcd_farkas';M='dominated_modular_equalities';C='dominated_cg_farkas'
 pt('triple_repeated_site',F,lambda q:q['sites'].__setitem__(1,q['sites'][0]))
 pt('triple_wrong_site_orientation',F,lambda q:q['sites'].__setitem__(0,v.p.rc(q['sites'][0])))
 pt('triple_unavailable_cycle_prefix',F,lambda q:q.__setitem__('cycle_count',len(data)+1))
 pt('triple_noninteger_cycle_prefix',F,lambda q:q.__setitem__('cycle_count',1.5))
 pt('dual_changed_numerator',F,lambda q:q['multipliers'][0].__setitem__(1,q['multipliers'][0][1]+1))
 pt('dual_zero_denominator',F,lambda q:q['multipliers'][0].__setitem__(2,0))
 pt('dual_negative_multiplier',F,lambda q:q['multipliers'][0].__setitem__(1,-1))
 pt('dual_repeated_row',F,lambda q:q['multipliers'].append(q['multipliers'][0].copy()))
 pt('dual_unavailable_row',F,lambda q:q['multipliers'][0].__setitem__(0,q['cycle_count']+256))
 pt('strip_forged_bound',S,lambda q:q['lower']['value'].__setitem__(0,q['lower']['value'][0]+1))
 pt('strip_wrong_linear_form',S,lambda q:q['form'][0].__setitem__(1,q['form'][0][1]+1))
 pt('strip_zero_denominator',S,lambda q:q['upper_negative']['value'].__setitem__(1,0))
 def wrong_sign(q):
  for i,(j,b) in enumerate(q['fixed_gamma']):
   if any(r['gamma'][j] for r in data[:q['cycle_count']]):q['fixed_gamma'][i][1]=1-b;return
  raise RuntimeError('NO_NONZERO_FIXED_BIT')
 pt('dominance_wrong_sign',D,wrong_sign)
 pt('dominance_nonbinary_fix',D,lambda q:q['fixed_gamma'][0].__setitem__(1,2))
 pt('dominance_missing_bit',D,lambda q:q['fixed_gamma'].pop())
 pt('dominance_duplicate_free_bit',D,lambda q:q['free_gamma'].append(q['free_gamma'][0]))
 pt('modular_wrong_modulus',M,lambda q:q.__setitem__('modulus',1))
 pt('modular_zero_combination_term',M,lambda q:q['equalities'][0].__setitem__('coefficient',0))
 pt('modular_forged_reverse_proof',M,lambda q:q['equalities'][0]['reverse_multipliers'][0].__setitem__(1,q['equalities'][0]['reverse_multipliers'][0][1]+1))
 pt('cg_wrong_left_side',C,lambda q:q['cuts'][0]['coefficients'][0].__setitem__(1,q['cuts'][0]['coefficients'][0][1]+1))
 pt('cg_wrong_upward_rounding',C,lambda q:q['cuts'][0].__setitem__('rhs',q['cuts'][0]['rhs']+1))
 pt('cg_downward_rounding',C,lambda q:q['cuts'][0].__setitem__('rhs',q['cuts'][0]['rhs']-1))
 pt('cg_noninteger_rhs',C,lambda q:q['cuts'][0].__setitem__('rhs',.5))
 pt('cg_forward_reference',C,lambda q:q['cuts'][0]['multipliers'][0].__setitem__(0,q['cycle_count']+2*len(q['free_gamma'])))
 pt('cg_missing_cut',C,lambda q:q['cuts'].pop())
 pt('cg_negative_weight',C,lambda q:q['cuts'][0]['multipliers'][0].__setitem__(1,-1))
 def whole(name,edit):
  q=change(allp,edit)
  with tempfile.TemporaryDirectory() as d:
   path=Path(d)/'bad.json.gz'
   with gzip.open(path,'wt') as f:json.dump(q,f)
   reject(name,lambda:v.verify_cases(base,stable,kernel['candidate_triples'],path))
 whole('allcases_wrong_coefficient_scope',lambda q:q.__setitem__('site_coefficients','bounded [-2,2] only'))
 whole('allcases_omitted_triple',lambda q:q['proofs'].pop())
 def omitted_reported(q):r=q['proofs'].pop();q['unresolved_supports']=[r['sites']]
 whole('allcases_partial_certificate_cannot_close_gate',omitted_reported)
 whole('allcases_duplicate_triple',lambda q:q['proofs'].append(q['proofs'][0]))
 read=v.load(ROOT/'certificates/readable_refutation.json')
 reject('readable_wrong_rounding_fraction',lambda:v.verify_readable(base,stable,allp,change(read,lambda q:q['pre_rounding_rhs'].__setitem__(0,1))))
 model=v.load(ROOT/'models/complete_T107_round8.json.gz')
 reject('global_wrong_period',lambda:vm.verify(change(model,lambda q:q.__setitem__('T',108))))
 reject('global_missing_raw_walk',lambda:vm.verify(change(model,lambda q:q['rows'].pop())))
 reject('global_false_no_site_cap',lambda:vm.verify(change(model,lambda q:q['round8_metadata'].__setitem__('site_cap',3))))
 run=v.load(ROOT/'discovery/results/adaptive803_stable/result.json')
 reject('search_actual_activity_wrong',lambda:audit.audit_adaptive(base,stable,change(run,lambda q:q['audit'][0].__setitem__('actual_sites',1))))
 reject('search_false_stable_change',lambda:audit.audit_adaptive(base,stable,change(run,lambda q:q['audit'][0].__setitem__('stable_changes',0))))
 reject('search_forged_deficit',lambda:audit.audit_adaptive(base,stable,change(run,lambda q:q['audit'][0].__setitem__('deficit',0))))
 reject('search_forged_cochain_hash',lambda:audit.audit_adaptive(base,stable,change(run,lambda q:q['audit'][0].__setitem__('cochain_sha256_le_i64','0'*64))))
 reject('closed_walk_removed_edge',lambda:v.rows_from_words([allp['cycles'][0][1:]],base,stable))
 out={'status':'PASS','mutations':len(results),'all_rejected':True,'tests':results,'seconds':time.time()-start}
 if a.out:a.out.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
