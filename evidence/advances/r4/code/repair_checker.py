#!/usr/bin/env python3
"""Certify the critical-cycle repair distance and a non-MDS repair witness."""
from collections import Counter
import graph_checker as g
from radius_checker import all_domains
from atlas_checker import verify_atlas

def verify_repair(base, bc, atlas, candidate, c):
    g.require(c['schema']=='rc-mds-critical-repair-gap-v1','REPAIR_SCHEMA')
    g.require(c['baseline_sha256']==g.digest_set(base) and c['candidate_sha256']==g.digest_set(candidate),'REPAIR_DIGEST')
    g.require(c['target_hits']==2,'REPAIR_TARGET')
    g.require(len(atlas['components'])==1,'REPAIR_COMPONENT_SCOPE')
    at=atlas['components'][0];total=at['simple_cycle_count'];weights=at['edge_cycle_counts']
    # Reconstruct the whole tight graph and every exact cycle multiplicity.
    got=verify_atlas(base,bc,atlas,candidate,c['candidate_histogram'])
    g.require(c['cycles']==total and all(int(k)>=2 for k in got['cross_histogram']),'REPAIR_UPPER')
    g.require(len(candidate)==8230 and {g.rc(w) for w in candidate}==candidate,'REPAIR_PCR_RC')
    g.necklaces_check(candidate);groups=all_domains(base);exceptions={}
    for rec in c['exceptional_single_alternatives']:
        key=(rec['group'],rec['option']);g.require(key not in exceptions,'REPAIR_DUPLICATE_EXCEPTION');exceptions[key]=rec
    count=weak=0;wanted=set();chosen=[]
    for i,gr in enumerate(groups):
        old=gr['options'][gr['base_option']]
        selected=[j for j,op in enumerate(gr['options']) if set(op)<=candidate]
        g.require(len(selected)==1,'REPAIR_CANDIDATE_CHOICE')
        if selected[0]!=gr['base_option']:chosen.append([i,selected[0]])
        for j,new in enumerate(gr['options']):
            if j==gr['base_option']:continue
            count+=1;gain=sum(weights.get(w,0) for w in new)-sum(weights.get(w,0) for w in old)
            if gain<total:weak+=1;continue
            key=(i,j);wanted.add(key);g.require(key in exceptions,'REPAIR_EXCEPTION_COVERAGE');rec=exceptions[key]
            cyc=rec['cycle'];updated=(base-set(old))|set(new)
            g.require(len(cyc)==at['T'] and all(w in weights for w in cyc),'REPAIR_EXCEPTION_NOT_CRITICAL')
            _,_,hit=g.walk_check(cyc,updated,cycle=True)
            g.require(hit<2 and rec['hits']==hit and rec['weighted_gain']==gain,'REPAIR_EXCEPTION_HIT')
    g.require(set(exceptions)==wanted and c['aggregate_exclusions']==weak and c['alternative_count']==count,'REPAIR_LOWER_COVERAGE')
    g.require(len(chosen)==2 and chosen==c['groups_changed'],'REPAIR_TWO_CHANGES')
    _,_,hit=g.walk_check(c['uncut_cycle'],candidate,cycle=True)
    g.require(hit==0,'REPAIR_FALSE_CYCLE')
    return {'status':'PASS','critical_cycle_repair_distance':2,'RC_PCR_groups':len(groups),
      'single_alternatives':count,'aggregate_exclusions':weak,'explicit_exceptions':len(wanted),
      'critical_cycle_count':total,'candidate_cycle_hit_histogram':got['cross_histogram'],
      'residual_uncut_cycle_length':len(c['uncut_cycle']),'candidate_is_MDS':False}
