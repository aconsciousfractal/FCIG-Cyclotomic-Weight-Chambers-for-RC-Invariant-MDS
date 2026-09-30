"""Replay only the mathematical dependencies of the stated local exclusions."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('round',type=int,choices=[2,4,5]);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    root=ROOT/f'r{a.round}';sys.path.insert(0,str(root/'code'))
    import verify as v
    import graph_checker as g
    if a.round==2:
        result=v.local_exclusion(g.read_set(root/'witnesses/M121.words.txt'),v.load(root/'certificates/local_exclusion_99PCR.json'))
    elif a.round==4:
        from sparse_checker import verify_sparse
        from fractional_checker import verify_fractional
        from repair_checker import verify_repair
        base=g.read_set(root/'inputs/M111.words.txt');sparse=v.load(root/'certificates/radius4_M111.json.gz')
        result={'radius':verify_sparse(base,sparse),'fractional':verify_fractional(base,sparse['cuts'],v.load(root/'certificates/fractional_cut_relaxation.json')),
                'repair':verify_repair(base,v.load(root/'inputs/M111.json'),v.load(root/'inputs/M111_atlas.json'),g.read_set(root/'witnesses/critical_repair_not_MDS.words.txt'),v.load(root/'certificates/critical_repair_gap.json'))}
    else:
        from atlas_checker import verify_atlas
        base=g.read_set(root/'inputs/M109.words.txt');bc=v.load(root/'inputs/M109.json')
        result={'radius':v.component_neighborhood(base,v.load(root/'certificates/component_radius3_M109.json.gz'))}
        sets={};certs={}
        for name,row in v.load(root/'certificates/catalog.json').items():
            sets[name]=g.read_set(root/row['words']);certs[name]=v.load(root/row['certificate']);g.set_certificate_check(sets[name],certs[name])
        result['components']=v.panel_check(sets,certs,v.load(root/'certificates/two_component_panel.json'))
        patch=v.load(root/'certificates/acyclic_trap_patch.json')
        result['trap']=verify_atlas(base,bc,v.load(root/'inputs/M109_atlas.json'),sets['Acyclic_trap_111'],patch['old_critical_histogram'])
        result['atlases']={row['name']:verify_atlas(sets[row['name']],certs[row['name']],v.load(root/row['atlas'])) for row in v.load(root/'certificates/atlas_catalog.json')}
    output={'status':'PASS','round':a.round,'result':result}
    a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n');print(json.dumps(output))
if __name__=='__main__':main()
