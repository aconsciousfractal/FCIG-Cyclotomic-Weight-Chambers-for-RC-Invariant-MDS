"""Recompute all 42 RC integer retimings from exact ordinary potentials.
This is a producer replay, distinct from the no-import certificate checker.
"""
from pathlib import Path
import json,hashlib
import retime as rt

def run(root):
    checked=[]
    for path in sorted((root/'inputs/ordinary_potentials').glob('*.json')):
        rec=json.loads(path.read_text());i=rec['component']
        m=rt.read_bits(root/f"inputs/{rec['source']}.bits",65536)
        for ep in (0,1):
            result=rt.rounded_rc(m,4,8,rec['T'],rec['potential'],ep)
            rt.require(result['status']=='FEASIBLE_RC','ROUNDING_STATUS')
            name=f'k8_component_{i:02d}_sector_{ep}'
            expected=json.loads((root/f'certificates/{name}.json').read_text())
            rt.require(result['H']==expected['H'] and result['cut_potential']==expected['cut_potential'],'ROUNDING_COORDINATES')
            rt.require(result['indicator']==rt.read_bits(root/expected['bits'],65536),'ROUNDING_INDICATOR')
            checked.append(name)
    rt.require(len(checked)==42,'ROUNDING_COVERAGE')
    return {'status':'PASS','ordinary_potentials':21,'retimings':42,'names':checked}

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    print(json.dumps(run(root),sort_keys=True,indent=2))
