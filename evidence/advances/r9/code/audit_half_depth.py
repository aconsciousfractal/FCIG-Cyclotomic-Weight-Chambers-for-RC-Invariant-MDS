"""Exact global validity of six two-arm cuts, and their retrospective coverage."""
from verify_round9 import *

def validate(c,stable_pairs,selectable):
    req(c['type']=='two_arm_half_depth' and c['target_L']==106,'HALF_DEPTH_TYPE')
    fixed=[decode(e) for e in c['fixed_arcs']];req(frozenset(fixed) in stable_pairs and len(set(fixed))==2,'HALF_DEPTH_STABLE_PAIR')
    paths=[[decode(e) for e in pp] for pp in c['paths']];req(len(paths)==2,'HALF_DEPTH_ARMS')
    for e,p in zip(fixed,paths):
        req(rc(e)==e and len(p)==53,'HALF_DEPTH_53')
        walk_check(p);req(p[-1][1:]==e[:-1],'HALF_DEPTH_TAIL')
        walk_check(p+[e]+[rc(w) for w in reversed(p)])
    words=set(sum(paths,[]));req(words=={decode(e) for e in c['cut_words']},'HALF_DEPTH_SUPPORT')
    req(all(rc(w)!=w for w in words),'HALF_DEPTH_NONFIXED')
    req(words&selectable=={decode(e) for e in c['selectable_support']},'HALF_DEPTH_SELECTABLE')
    return words&selectable

def main():
    t=time.time();base=read_words(ROOT/'inputs/M107.words.txt');stable=necklaces(base);pairs={frozenset(w for w in p if rc(w)==w) for p in stable};selectable=EDGESET-set().union(*stable)
    cuts=json.loads((ROOT/'certificates/half_depth_cuts.json').read_text());supports=[validate(c,pairs,selectable) for c in cuts];old=json.loads((ROOT/'certificates/boolean_master_cuts.json').read_text());points={}
    for ob in old:points[ob['iteration']]=((base-set().union(*stable))-{decode(i) for i in ob['outer_removed_selected']})|{decode(i) for i in ob['outer_changed_selected']}
    result=json.loads((ROOT/'logs/half_depth_filter.json').read_text());viol={str(i):[j for j,s in enumerate(supports) if not (s&m)] for i,m in points.items()}
    req(viol==result['point_violations'] and sum(bool(x) for x in viol.values())==result['retrospectively_rejected'],'HALF_DEPTH_FILTER')
    req(result['global_cuts']==len(cuts)==6,'HALF_DEPTH_COUNT')
    out={'status':'PASS','globally_valid_cuts':len(cuts),'retrospective_points':len(points),'retrospectively_rejected':sum(bool(x) for x in viol.values()),'retrospective_only':True,'seconds':time.time()-t}
    print(json.dumps(out,indent=2));(OUT/'half_depth_audit.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
