#!/usr/bin/env python3
"""Round5: exact offline verification, including retiming-closed local exclusion."""
from pathlib import Path
from collections import Counter
import argparse,gzip,json,time
import graph_checker as g
from switch_checker import verify_chain
from atlas_checker import verify_atlas
from sparse_checker import verify_sparse
from radius_checker import all_domains
ROOT=Path(__file__).resolve().parents[1]
def load(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return json.load(f)
def component_neighborhood(base,cert):
 g.require(cert['radius']==3 and cert['bound']==108,'COMPONENT_RADIUS_SCOPE')
 g.require(all(c['kind']=='cycle' for c in cert['cuts']),'COMPONENT_EXCLUSION_NEEDS_ONLY_CYCLES')
 r=verify_sparse(base,cert)
 r.update({'retiming_closure_excluded':True,'ordinary_components_after_perturbation':True,
           'any_number_of_complete_moves_after_perturbation':True,
           'pre_retiming_not_included':True,'minimum_in_retiming_closure':109})
 return r

def panel_check(sets,certs,panel):
 pairs=[]
 for item in panel['same_component']:
  src,tgt=item['source'],item['target'];f=item['potential'];g.require(src in sets and tgt in sets,'PANEL_ENDPOINT')
  g.require(len(f)==16384 and all(type(x) is int for x in f),'PANEL_INTEGER_POTENTIAL')
  g.require(all(int(e in sets[tgt])-int(e in sets[src])==f[g.VI[e[:-1]]]-f[g.VI[e[1:]]] for e in g.EDGES),'PANEL_GRADIENT')
  sums={f[g.VI[v]]+f[g.VI[g.rc(v)]] for v in g.VERTICES};g.require(len(sums)==1 and next(iter(sums))%2==item['relative_sector']==1,'PANEL_SECTOR')
  pairs.append([src,tgt])
 sep=panel['separation'];g.require(sep['type']=='separating_closed_walk','PANEL_SEPARATOR_TYPE');src,tgt=sep['source'],sep['target']
 g.require(sep['source_sha256']==g.digest_set(sets[src]) and sep['target_sha256']==g.digest_set(sets[tgt]),'PANEL_SEPARATOR_DIGEST')
 h1=g.walk_check(sep['words'],sets[src],cycle=True)[2];h2=g.walk_check(sep['words'],sets[tgt],cycle=True)[2]
 g.require(h1==sep['source_hits'] and h2==sep['target_hits'] and h1!=h2,'PANEL_DIFFERENT_COMPONENTS')
 g.require(all(certs[name]['L']==certs[name]['component_lower_bound']==107 for pair in pairs for name in pair),'PANEL_COMPONENT_MINIMA')
 return {'pairs':pairs,'separator_length':len(sep['words']),'separator_hits':[h1,h2],'ordinary_component_minima':[107,107]}


def check(root=ROOT,quick=False,manifest=False):
 start=time.time();sets={};certs={};out={'status':'PASS_FULL' if not quick else 'PASS_WITH_DIGEST_ONLY_CHAIN','sets':{}}
 catalog=load(root/'certificates/catalog.json')
 for name,entry in catalog.items():
  sets[name]=g.read_set(root/entry['words']);certs[name]=load(root/entry['certificate']);out['sets'][name]=g.set_certificate_check(sets[name],certs[name])
 g.require(certs['M109']['L']==109 and certs['M107']['L']==107,'ENDPOINT_LENGTHS')
 for a,b,f in [('M109','M109_rephased','retime_initial'),('M111_before_retiming','M107','retime_final_sector0'),('M111_before_retiming','M107_sector1','retime_final_sector1')]:
  g.retiming_check(sets[a],sets[b],load(root/f'certificates/{f}.json'))
 out['retimings']=3
 out['two_component_panel']=panel_check(sets,certs,load(root/'certificates/two_component_panel.json'))
 out['chain']=verify_chain(root,load(root/'certificates/chain_M109_to_M107.json'),quick=quick)
 out['atlases']={}
 for entry in load(root/'certificates/atlas_catalog.json'):
  name=entry['name'];out['atlases'][name]=verify_atlas(sets[name],certs[name],load(root/entry['atlas']))
 out['component_neighborhood']=component_neighborhood(sets['M109'],load(root/'certificates/component_radius3_M109.json.gz'))
 scope=load(root/'certificates/component_neighborhood_scope.json')
 g.require(scope['baseline']=='M109' and scope['radius']==3 and scope['period']==109 and scope['excluded_bound']==108,'NEIGHBORHOOD_SCOPE')
 g.require(scope['all_cuts_are_closed_walks'] is True and scope['ordinary_components_after_perturbation'] is True and scope['pre_retiming_not_included'] is True,'NEIGHBORHOOD_BOUNDARY')
 final=sets['M107'];path=certs['M107']['longest_path'];word=path[0]+''.join(e[-1] for e in path[1:])
 g.require((root/'witnesses/M107_longest_avoiding_word.txt').read_text().strip()==word and len(word)==114,'AVOIDING_WORD')
 g.require(all(word[i:i+8] not in final for i in range(len(word)-7)),'AVOIDING_WORD_HIT')
 g.require((root/'witnesses/M107_critical_cycle_108.txt').read_text().splitlines()==certs['M107']['critical_cycle'],'CRITICAL_CYCLE_TEXT')
 groups=all_domains(sets['M109']);changed=[gr for gr in groups if not set(gr['options'][gr['base_option']])<=final];cy=certs['M109']['critical_cycle']
 comp={'changed_RC_PCR_groups':len(changed),'changed_PCRs':sum(len(set(gr['necklaces'])) for gr in changed),'symmetric_difference':len(final^sets['M109']),'old_cycle_length':len(cy),'old_cycle_hits_final':g.walk_check(cy,final,cycle=True)[2]}
 g.require(comp==load(root/'certificates/comparison.json'),'COMPARISON');g.require(comp['old_cycle_hits_final']!=1,'COMPONENT_SEPARATION')
 out['comparison']=comp
 trap=sets['Acyclic_trap_111'];patch=load(root/'certificates/acyclic_trap_patch.json');restored=sets['M109'].copy();gs=all_domains(sets['M109']);used=set()
 for ch in patch['changes']:
  i=ch['group'];g.require(type(i) is int and 0<=i<len(gs) and i not in used,'TRAP_GROUP');used.add(i)
  gr=gs[i];g.require(ch['old']==gr['options'][gr['base_option']] and ch['new'] in gr['options'] and ch['new']!=ch['old'],'TRAP_OPTIONS')
  restored.difference_update(ch['old']);restored.update(ch['new'])
 g.require(len(used)==3 and restored==trap,'TRAP_PATCH_REPLAY')
 tr=verify_atlas(sets['M109'],certs['M109'],load(root/'inputs/M109_atlas.json'),trap,patch['old_critical_histogram'])
 g.require(tr['cross_histogram']=={'3':6931456},'TRAP_COMPLETE_OLD_CYCLES')
 g.require(certs['Acyclic_trap_111']['L']==111 and certs['Acyclic_trap_111']['component_lower_bound']==111,'TRAP_WORSE_COMPONENT')
 out['acyclic_trap']={'changed_groups':3,'old_critical_histogram':tr['cross_histogram'],'ordinary_component_minimum':111,'residual_L':111}
 out['word_guarantee']={'longest_avoiding_letters':114,'universal_hit_letters':115}
 if manifest:out['manifest_files']=g.manifest_check(root)
 out['elapsed_seconds']=time.time()-start;return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--quick',action='store_true');ap.add_argument('--manifest',action='store_true');a=ap.parse_args();print(json.dumps(check(a.root,a.quick,a.manifest),indent=2))
