"""Selected negative controls of actual certificate validators; not a proof of
exhaustive mutation resistance. No solver, no floats, no ignored assertions.
"""
from core import *
from components import replay_pair,verify_separator
from verify import verify_sparse_certificate
from copy import deepcopy
import argparse

def main():
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path);p.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results');args=p.parse_args()
    parent,tmp=locate_parent(args.parent);receipts=[]
    try:
        A,ts=parent_data(parent,8);ids={r['id']:i for i,r in enumerate(A['rows'])}
        c=json.loads((args.out/'components_k8.json').read_text())
        pos=next(x for x in c['positive_geodesics'] if len(x['moves'])>=3)
        neg=c['negative_cycles'][0]
        mix=json.loads((args.out/'mixed_k8_r064.json').read_text())['comparisons'][2]
        def trial(name,fn):
            try:fn()
            except CheckError as e:receipts.append({'name':name,'result':'REJECTED','output':str(e)});print(name,'REJECTED',e);return
            raise CheckError('MUTATION_ACCEPTED:'+name)
        def pv(z):return replay_pair(ts[ids[pos['source']]],ts[ids[pos['target']]],8,z)
        def nv(z):return verify_separator(ts[ids[neg['source']]],ts[ids[neg['target']]],8,z)
        def mv(z):return verify_sparse_certificate(ts[ids[mix['source']]],ts[ids[mix['target']]],z)
        z=deepcopy(pos);z['moves'].pop();trial('drop_pair',lambda:pv(z))
        z=deepcopy(pos);z['moves'][1]=z['moves'][0];trial('repeat_pair',lambda:pv(z))
        z=deepcopy(pos);z['cut_constant']+=1;trial('parity_flip',lambda:pv(z))
        z=deepcopy(pos);z['minimum_paired_moves']+=1;trial('false_geodesic_length',lambda:pv(z))
        z=deepcopy(pos);z['moves'][0]=4**7;trial('out_of_carrier_pair',lambda:pv(z))
        z=deepcopy(neg);z['hits_source']+=1;trial('false_cycle_hit_count',lambda:nv(z))
        z=deepcopy(neg);z['edge_codes'][0]=z['edge_codes'][1];trial('repeated_cycle_vertex',lambda:nv(z))
        z=deepcopy(neg);z['edge_codes'].reverse();trial('reverse_directed_cycle',lambda:nv(z))
        z=deepcopy(neg);z['length']=16;trial('oversized_separator',lambda:nv(z))
        z=deepcopy(mix);z['patch'][0][1]*=-1;trial('patch_wrong_sign',lambda:mv(z))
        z=deepcopy(mix);z['patch'].append(z['patch'][0]);trial('duplicate_patch_edge',lambda:mv(z))
        z=deepcopy(mix);z['patch'].pop();trial('missing_patch_edge',lambda:mv(z))
        z=deepcopy(mix);z['changed_null_pcrs']+=1;trial('false_null_cycle_count',lambda:mv(z))
        z=deepcopy(mix);z['regular_midpoint_states']-=1;z['external_companion_pairs']+=4;trial('coherent_total_wrong_breakdown',lambda:mv(z))
        z=deepcopy(mix);z['symmetric_difference']-=2;trial('wrong_total_jump',lambda:mv(z))
        check(len(receipts)==15,'MUTATION_PANEL_COUNT');dump(args.out/'mutations_round3.json',{'status':'PASS','rejected':len(receipts),'receipts':receipts})
    finally:
        if tmp:tmp.cleanup()
if __name__=='__main__':main()
