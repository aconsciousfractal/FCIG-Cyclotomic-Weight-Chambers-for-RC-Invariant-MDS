#!/usr/bin/env python3
"""Read every LP row; verify semantic coefficients independently of its writer."""
from pathlib import Path
from collections import Counter
import gzip,json,re
import graph_checker as g
from atlas_checker import verify_atlas
R=Path(__file__).resolve().parents[1]

def parsed(row):
    label,tail=row.split(':',1)
    expr,relation,rhs=re.split(r' (>=|=) ',tail.strip())
    tokens=expr.split();g.require(len(tokens)%3==0,'MODEL_TOKEN_COUNT');co=Counter()
    for k in range(0,len(tokens),3):
        sign,n,name=tokens[k:k+3];g.require(sign in ('+','-') and n.isdigit(),'MODEL_COEFFICIENT')
        co[name]+=int(n)*(1 if sign=='+' else -1)
    return label.strip(),{k:v for k,v in co.items() if v},relation,int(rhs)

def main(path=None):
    bound=108;reps=sorted({min(w,g.rc(w)) for w in g.EDGES});idx={w:i for i,w in enumerate(reps)}
    x=lambda w:'x'+str(idx[min(w,g.rc(w))]);pcrs=[];seen=set()
    for w in g.EDGES:
        if w in seen:continue
        orb={w[i:]+w[:i] for i in range(8)};seen|=orb;pcrs.append(dict(Counter(x(z) for z in orb)))
    cycles=set()
    for e in json.loads((R/'certificates/catalog.json').read_text()).values():
        c=json.loads((R/e['certificate']).read_text())
        if 'critical_cycle' not in c:continue
        cy=c['critical_cycle'];g.walk_check(cy,set(),cycle=True)
        cycles.add((tuple(sorted(Counter(x(z) for z in cy).items())),(len(cy)+bound)//(bound+1)))
    aggregates=set()
    for entry in json.loads((R/'certificates/atlas_catalog.json').read_text()):
        atlas=json.loads((R/entry['atlas']).read_text())
        aset=g.read_set(R/entry['words']);acert=json.loads((R/entry['certificate']).read_text())
        g.set_certificate_check(aset,acert);verify_atlas(aset,acert,atlas)
        for c in atlas['components']:
            weights=Counter()
            for w,n in c['edge_cycle_counts'].items():weights[x(w)]+=n
            aggregates.add((tuple(sorted(weights.items())),((c['T']+bound)//(bound+1))*c['simple_cycle_count']))
    at=json.loads((R/'certificates/M109_atlas.json').read_text())['components'][0]
    anchor=at['anchor'];yy={v:'y'+str(i) for i,v in enumerate(at['vertices'])}; threshold=(at['T']+bound)//(bound+1)
    nExtended=0
    nPCR=nClock=0;gotCycles=set();gotAggregates=set();bounds=[];binary=[];mode=None
    for row in gzip.open(path if path is not None else R/'models/global_L108.lp.gz','rt'):
        row=row.strip()
        if not row or row.startswith('\\'):continue
        if row in ('Minimize','Subject To','Bounds','Binary','End'):mode=row;continue
        if mode=='Minimize':g.require(row=='obj: 0 x0','MODEL_OBJECTIVE')
        elif mode=='Subject To':
            label,terms,relation,rhs=parsed(row)
            if label.startswith('pcr_'):
                g.require(label=='pcr_'+str(nPCR) and terms==pcrs[nPCR] and relation=='=' and rhs==1,'MODEL_PCR_ROW');nPCR+=1
            elif label.startswith('clock_'):
                w=g.EDGES[nClock];e=Counter({x(w):bound+1});e['h'+str(g.VI[w[1:]])]+=1;e['h'+str(g.VI[w[:-1]])]-=1
                expected={k:v for k,v in e.items() if v}
                g.require(label=='clock_'+str(nClock) and terms==expected and relation=='>=' and rhs==1,'MODEL_CLOCK_ROW');nClock+=1
            elif label.startswith('tightpath_'):
                w=at['edges'][nExtended];expected=Counter({x(w):1})
                expected[yy[w[:-1]]]+=1;expected['ysink' if w[1:]==anchor else yy[w[1:]]]-=1
                expected={k:v for k,v in expected.items() if v}
                g.require(label=='tightpath_'+str(nExtended) and terms==expected and relation=='>=' and rhs==0,'MODEL_EXTENDED_ROW');nExtended+=1
            elif label.startswith('aggregate_'):
                g.require(relation=='>=','MODEL_AGGREGATE_RELATION');gotAggregates.add((tuple(sorted(terms.items())),rhs))
            else:
                g.require(label.startswith('cyclecut_') and relation=='>=','MODEL_UNKNOWN_ROW');gotCycles.add((tuple(sorted(terms.items())),rhs))
        elif mode=='Bounds':bounds.append(row)
        elif mode=='Binary':binary.append(row)
        else:raise g.VerificationError('MODEL_SECTION')
    g.require(nPCR==8230 and nClock==65536,'MODEL_ROW_COVERAGE')
    g.require(cycles==gotCycles,'MODEL_CYCLE_CUTS')
    g.require(aggregates==gotAggregates,'MODEL_AGGREGATE_CUTS')
    g.require(nExtended==len(at['edges']),'MODEL_EXTENDED_COVERAGE')
    expected_bounds=[f'0 <= h{i} <= {bound}' for i in range(16384)]
    expected_bounds += [f'0 <= {yy[z]} <= {0 if z==anchor else threshold}' for z in at['vertices']]
    expected_bounds += [f'{threshold} <= ysink <= {threshold}']
    g.require(bounds==expected_bounds,'MODEL_BOUNDS')
    g.require(binary==[f'x{i}' for i in range(len(reps))],'MODEL_BINARY_COVERAGE')
    # Independently evaluate the base system at B=111 on the known witness.
    m=g.read_set(R/'witnesses/M109.words.txt');c=json.loads((R/'certificates/M109.json').read_text());h=c['rank']
    g.require(all(0<=z<=109 for z in h),'MODEL_POSITIVE_RANK')
    g.necklaces_check(m);g.require({g.rc(w) for w in m}==m,'MODEL_POSITIVE_RC')
    g.require(all(h[g.VI[w[1:]]]-h[g.VI[w[:-1]]]+110*int(w in m)>=1 for w in g.EDGES),'MODEL_POSITIVE_CLOCK')
    print(json.dumps({'status':'PASS','global_bound':bound,'PCR_rows':nPCR,'clock_rows':nClock,'RC_binary_variables':len(reps),'ranks':16384,'extended_potentials':len(at['vertices'])+1,'extended_rows':nExtended,'cycles_encoded':at['simple_cycle_count'],'valid_cycle_cuts':len(cycles),'aggregate_cycle_cuts':len(aggregates),'known_M109_assignment_at_109':'PASS','solver_run':False,'global_infeasibility_certified':False},indent=2))
if __name__=='__main__':main()
