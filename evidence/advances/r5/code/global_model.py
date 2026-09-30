#!/usr/bin/env python3
"""Emit a complete global feasibility MILP in LP format (not a component model).

Only stdlib is needed. No solver is called. Infeasibility is NOT asserted.
Binary x variables index reverse-complement orbits of eight-letter words;
continuous h variables are residual ranks on seven-letter vertices.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse, gzip, json
import graph_checker as v


def linear_expression(terms: dict[str,int]) -> str:
    out=[]
    for name,coef in terms.items():
        if coef:
            out.append(('+' if coef>0 else '-')+' '+str(abs(coef))+' '+name)
    return ' '.join(out) if out else '0'


def write_model(bound: int, path: Path, *, cycle_cuts: bool = True) -> dict:
    v.require(type(bound) is int and bound>=0,'MODEL_BOUND')
    path.parent.mkdir(parents=True,exist_ok=True)
    canonical={w:min(w,v.rc(w)) for w in v.EDGES}
    reps=sorted(set(canonical.values()))
    variables={w:'x'+str(i) for i,w in enumerate(reps)}
    seen=set();pcrs=[]
    for w in v.EDGES:
        if w in seen:continue
        orbit={w[i:]+w[:i] for i in range(v.K)}
        seen.update(orbit)
        pcrs.append(Counter(variables[canonical[z]] for z in sorted(orbit)))
    opener=gzip.open if path.suffix=='.gz' else open
    cut_count=0
    aggregate_count=0
    with opener(path,'wt',encoding='ascii',newline='\n') as f:
        f.write('\\ DNA k=8; global L <= '+str(bound)+'; exact cardinality 8230\n')
        f.write('Minimize\n obj: 0 x0\nSubject To\n')
        for i,terms in enumerate(pcrs):
            f.write(f' pcr_{i}: '+linear_expression(terms)+' = 1\n')
        for e,w in enumerate(v.EDGES):
            terms=Counter()
            terms['h'+str(v.VI[w[1:]])]+=1
            terms['h'+str(v.VI[w[:-1]])]-=1
            terms[variables[canonical[w]]]+=bound+1
            f.write(f' clock_{e}: '+linear_expression(terms)+' >= 1\n')
        if cycle_cuts:
            catalog=json.loads((v.ROOT/'certificates/catalog.json').read_text())
            cut_seen=set()
            for name,entry in catalog.items():
                cert=json.loads((v.ROOT/entry['certificate']).read_text())
                if 'critical_cycle' not in cert:continue
                cyc=cert['critical_cycle']
                # Counts aggregate if the closed walk revisits an RC orbit.
                terms=Counter(variables[canonical[w]] for w in cyc)
                rhs=(len(cyc)+bound)//(bound+1)
                key=(tuple(sorted(terms.items())),rhs)
                if key in cut_seen:continue
                cut_seen.add(key)
                f.write(f' cyclecut_{cut_count}: '+linear_expression(terms)+f' >= {rhs}\n')
                cut_count+=1
        if cycle_cuts:
            # These are valid weighted SUMS of the exact simple-cycle families.
            # An aggregate inequality is weaker than its full constituent family.
            for entry in json.loads((v.ROOT/'certificates/atlas_catalog.json').read_text()):
                atlas=json.loads((v.ROOT/entry['atlas']).read_text())
                for component in atlas['components']:
                    terms=Counter()
                    for w,n in component['edge_cycle_counts'].items():
                        terms[variables[canonical[w]]]+=n
                    rhs=((component['T']+bound)//(bound+1))*component['simple_cycle_count']
                    f.write(f' aggregate_{aggregate_count}: '+linear_expression(terms)+f' >= {rhs}\n')
                    aggregate_count+=1
        # Complete current critical-cycle family via an anchored shortest-path model.
        at=json.loads((v.ROOT/'certificates/M107_atlas.json').read_text())['components'][0]
        anchor=at['anchor']; yy={z:'y'+str(i) for i,z in enumerate(at['vertices'])}; sink='ysink'
        threshold=(at['T']+bound)//(bound+1)
        for i,w in enumerate(at['edges']):
            terms=Counter({variables[canonical[w]]:1});terms[yy[w[:-1]]]+=1
            terms[sink if w[1:]==anchor else yy[w[1:]]]-=1
            f.write(f' tightpath_{i}: '+linear_expression(terms)+' >= 0\n')
        f.write('Bounds\n')
        for i in range(len(v.VERTICES)):
            f.write(f' 0 <= h{i} <= {bound}\n')
        for z in at['vertices']:
            f.write(f' 0 <= {yy[z]} <= {0 if z==anchor else threshold}\n')
        f.write(f' {threshold} <= ysink <= {threshold}\n')
        f.write('Binary\n')
        for w in reps:
            f.write(' '+variables[w]+'\n')
        f.write('End\n')
    metadata={'q':4,'k':8,'bound_L':bound,'binary_RC_orbit_variables':len(reps),
      'continuous_rank_variables':len(v.VERTICES),'necklace_equalities':len(pcrs),
      'clock_inequalities':len(v.EDGES),'additional_valid_cycle_cuts':cut_count,
      'additional_aggregate_cycle_cuts':aggregate_count,
      'global_model':True,'solver_run':False,
      'extended_critical_potentials':len(at['vertices'])+1,'extended_critical_rows':len(at['edges']),
      'individual_cycles_encoded':at['simple_cycle_count'],'individual_cycle_threshold':threshold,
      'indexing':'lexicographic A<C<G<T; x indexes sorted min(w,RC(w)); h indexes seven-letter words'}
    path.with_name(path.name+'.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


def assignment_check(selected: set[str], rank: list[int], bound: int) -> bool:
    """Exact check of one proposed model assignment; no infeasibility inference."""
    v.require(len(rank)==len(v.VERTICES),'MODEL_RANK_SIZE')
    v.require(all(type(h) is int and 0<=h<=bound for h in rank),'MODEL_RANK_RANGE')
    v.require({v.rc(w) for w in selected}==selected,'MODEL_RC')
    v.necklaces_check(selected)
    for w in v.EDGES:
        v.require(rank[v.VI[w[1:]]]-rank[v.VI[w[:-1]]]+(bound+1)*int(w in selected)>=1,
                  'MODEL_CLOCK')
    return True


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bound',type=int,default=106)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--no-cycle-cuts',action='store_true')
    a=p.parse_args()
    print(json.dumps(write_model(a.bound,a.out,cycle_cuts=not a.no_cycle_cuts),indent=2))

if __name__=='__main__':main()
