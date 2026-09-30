"""Deliberate corruptions; focused exact checks exclude repeated full DAG traversals."""
import copy,json,time
from verify_round9 import *

def main():
    t=time.time();base=read_words(ROOT/'inputs/M107.words.txt');stable=necklaces(base);c=json.loads((ROOT/'certificates/stable_cube.json').read_text());tests=[]
    def cube(label,edit):
        a=copy.deepcopy(c);edit(a)
        try:verify_cube(a,base,stable,direct=False)
        except (ValueError,KeyError,IndexError) as e:tests.append({'case':label,'status':'REJECTED','output':str(e)});return
        raise ValueError('MUTATION_ACCEPTED:'+label)
    cube('remove_unary',lambda a:a['acyclic_unary_forbidden'].pop())
    cube('remove_pair',lambda a:a['acyclic_pairs_forbidden'].pop())
    cube('invent_pair',lambda a:a['acyclic_pairs_forbidden'].append([0,0]))
    cube('alter_cost',lambda a:a['costs'].__setitem__(0,a['costs'][0]+2))
    cube('alter_fixed_L',lambda a:a.__setitem__('fixed_L',106))
    cube('remove_fixed_path_edge',lambda a:a['fixed_longest_path'].pop())
    cube('wrong_base',lambda a:a.__setitem__('base_sha256','0'*64))
    cube('duplicate_fixed_edge',lambda a:a['edges'].__setitem__(0,a['edges'][1]))
    cube('out_of_range_edge',lambda a:a['edges'].__setitem__(0,65536))
    cube('remove_forced',lambda a:a['acyclic_forced_literals'].pop())
    cube('invent_free_variable',lambda a:a['free_variables'].append(0))
    cube('alter_census_L',lambda a:a['complete_census'][0].__setitem__('L',106))
    cube('alter_census_bit',lambda a:a['complete_census'][0]['bits'].__setitem__(0,1))
    cube('remove_census_row',lambda a:a['complete_census'].pop())
    cube('alter_profile',lambda a:a['length_profile'].__setitem__('107',9))
    cube('false_optimum',lambda a:a.__setitem__('optimum',106))
    cube('component_cycle_break',lambda a:a['complete_census'][0]['component_cycle'].pop())
    cube('component_hit_error',lambda a:a['complete_census'][0].__setitem__('component_hits',2))
    p=json.loads((ROOT/'certificates/boundary_small.json').read_text())
    for label,edit in [('drop_boundary',lambda a:a['boundary'].pop()),('drop_path_constraint',lambda a:a['compressed_paths'].pop()),('alter_path_length',lambda a:a['compressed_paths'][0].__setitem__(2,99)),('alter_prefix_bound',lambda a:a['boundary_bounds'][0].__setitem__(1,99)),('alter_fixed_path_bound',lambda a:a.__setitem__('fixed_L',0))]:
        a=copy.deepcopy(p);edit(a)
        try:verify_boundary(a,base)
        except (ValueError,KeyError,IndexError) as e:tests.append({'case':label,'status':'REJECTED','output':str(e)});continue
        raise ValueError('MUTATION_ACCEPTED:'+label)
    from audit_half_depth import validate as vh
    from audit_master import validate_cut as vm
    pairs={frozenset(w for w in pp if rc(w)==w) for pp in stable};selectable=EDGESET-set().union(*stable)
    hc=json.loads((ROOT/'certificates/half_depth_cuts.json').read_text())[0]
    for label,edit in [('half_depth_short_arm',lambda a:a['paths'][0].pop(0)),('half_depth_wrong_target',lambda a:a.__setitem__('target_L',107)),('half_depth_wrong_pair',lambda a:a['fixed_arcs'].__setitem__(1,a['fixed_arcs'][0])),('half_depth_missing_support',lambda a:a['cut_words'].pop()),('half_depth_missing_selectable',lambda a:a['selectable_support'].pop())]:
        a=copy.deepcopy(hc);edit(a)
        try:vh(a,pairs,selectable)
        except (ValueError,KeyError,IndexError) as e:tests.append({'case':label,'status':'REJECTED','output':str(e)});continue
        raise ValueError('MUTATION_ACCEPTED:'+label)
    mc=json.loads((ROOT/'certificates/boolean_master_cuts.json').read_text())[0]
    for label,edit in [('master_short_path',lambda a:a['walk'].pop()),('master_missing_support',lambda a:a['cut_words'].pop()),('master_bad_type',lambda a:a.__setitem__('type','UNSAT_GLOBAL'))]:
        a=copy.deepcopy(mc);edit(a)
        try:vm(a,selectable,pairs)
        except (ValueError,KeyError,IndexError) as e:tests.append({'case':label,'status':'REJECTED','output':str(e)});continue
        raise ValueError('MUTATION_ACCEPTED:'+label)
    out={'status':'PASS','rejected':len(tests),'tests':tests,'seconds':time.time()-t};(OUT/'mutation_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
