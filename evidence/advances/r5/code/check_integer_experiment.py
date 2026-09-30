from pathlib import Path
import sys,json
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'code'))
import graph_checker as g
from atlas_checker import component_counts
from radius_checker import all_domains
D=R/'discovery/runs/repair83';run=json.load(open(D/'result.json'));cuts=json.load(open(D/'cuts.json'));base=g.read_set(R/'inputs/M109.words.txt');at=json.load(open(R/'inputs/M109_atlas.json'))['components'][0];groups=all_domains(base);out=[]
for rec in run['iterations']:
 i=rec['iteration'];p=D/f'candidate_{i:03d}.words.txt'
 if not p.exists():
  g.require(rec['objective'] is None,'MISSING_INCUMBENT');continue
 S=g.read_set(p);g.require(len(S)==8230 and {g.rc(w) for w in S}==S,'INCUMBENT_DOMAIN');g.necklaces_check(S)
 changes=sum(not set(gr['options'][gr['base_option']])<=S for gr in groups);g.require(changes==rec['changes'],'CHANGES')
 hist=component_counts(at['vertices'],at['edges'],at['anchor'],at['T'],S)['histogram'];g.require(min(map(int,hist))>=2,'OLD_CRITICAL_CYCLE_NOT_REPAIRED')
 try:res=g.graph_check(S);ac=True
 except g.VerificationError as e:
  g.require(str(e)=='RESIDUAL_DIRECTED_CYCLE','UNEXPECTED_GRAPH_FAILURE');ac=False
 g.require(ac==rec['acyclic'],'GRAPH_RECEIPT')
 if ac:g.require(res['L']==rec['L'],'GRAPH_LENGTH_RECEIPT')
 first=rec['cycle_cuts'];new=rec['new_cycle_cuts'];cy=cuts[first]
 g.walk_check(cy['words'],S,cycle=True);hits=sum(w in S for w in cy['words']);g.require(len(cy['words'])>109*hits,'NEW_CYCLE_NOT_A_CLOCK_OBSTRUCTION')
 out.append({'iteration':i,'sha256_words':g.digest_set(S),'changed_groups':changes,'acyclic':ac,'L':res['L'] if ac else None,'old_critical_min_hits':min(map(int,hist)),'new_clock_obstruction':{'words':cy['words'],'length':len(cy['words']),'hits':hits}})
result={'status':'PASS_EXACT_INCUMBENT_AUDIT','T':109,'radius_budget':32,'solver_calls':len(run['iterations']),'incumbents_checked':len(out),'acyclic_incumbents':sum(r['acyclic'] for r in out),'cyclic_incumbents':sum(not r['acyclic'] for r in out),'termination':run['status'],'global_unsat_certified':False,'records':out}
print(json.dumps(result,indent=2))
