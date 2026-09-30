#!/usr/bin/env python3
"""Exact primal/dual certificate for the specified finite cut relaxation.
Not a relaxation result for every missing cycle or for the full rank model.
"""
from collections import Counter,defaultdict
from fractions import Fraction
import graph_checker as g
from radius_checker import all_domains

def verify_fractional(base, cuts, cert):
    g.require(cert['schema']=='rc-mds-finite-cut-lp-v1','LP_SCHEMA')
    g.require(cert['baseline_sha256']==g.digest_set(base),'LP_BASE')
    groups=all_domains(base);alph=[];inverse=defaultdict(list)
    for i,gr in enumerate(groups):
        old=set(gr['options'][gr['base_option']])
        for j,op in enumerate(gr['options']):
            if j==gr['base_option']:continue
            a=len(alph);alph.append((i,j));new=set(op)
            for w in old|new:
                d=int(w in new)-int(w in old)
                if d:inverse[w].append((a,d))
    ids={ij:a for a,ij in enumerate(alph)};x={};mass=defaultdict(Fraction)
    for row in cert['primal']:
        key=(row['group'],row['option']);g.require(key in ids and ids[key] not in x,'LP_PRIMAL_INDEX')
        value=Fraction(row['value']);g.require(value>=0,'LP_NEGATIVE_PRIMAL');x[ids[key]]=value;mass[key[0]]+=value
    g.require(all(t<=1 for t in mass.values()),'LP_GROUP_MASS')
    dual={}
    for row in cert['dual']:
        i=row['cut'];g.require(type(i) is int and 0<=i<len(cuts) and i not in dual,'LP_DUAL_INDEX')
        value=Fraction(row['weight']);g.require(value>=0,'LP_NEGATIVE_DUAL');dual[i]=value
    totals=defaultdict(Fraction);objective=Fraction(0)
    for i,cut in enumerate(cuts):
        kind=cut['kind'];g.require(kind in ('cycle','path'),'LP_CUT_KIND');words=cut['words']
        _,_,hit=g.walk_check(words,base,cycle=kind=='cycle')
        if kind=='path':g.require(len(words)==cert['bound']+1,'LP_PATH_LENGTH');rhs=1
        else:rhs=(len(words)+cert['bound'])//(cert['bound']+1)
        g.require(cut['rhs']==rhs,'LP_CUT_RHS');co=Counter()
        for w,k in Counter(words).items():
            for a,d in inverse[w]:co[a]+=k*d
        g.require(sum(Fraction(co[a])*v for a,v in x.items())>=rhs-hit,'LP_PRIMAL_CUT')
        if i in dual:
            for a,z in co.items():totals[a]+=dual[i]*z
            objective+=dual[i]*(rhs-hit)
    g.require(all(z<=1 for z in totals.values()),'LP_DUAL_COLUMN')
    optimum=Fraction(cert['optimum']);g.require(sum(x.values())==objective==optimum,'LP_OBJECTIVE_MATCH')
    return {'status':'PASS','specified_cut_count':len(cuts),'alternative_variables':len(alph),
      'exact_fractional_optimum':str(optimum),'primal_nonzeros':len(x),'dual_nonzeros':len(dual),
      'unrestricted_MDS_LP_claim':False}
