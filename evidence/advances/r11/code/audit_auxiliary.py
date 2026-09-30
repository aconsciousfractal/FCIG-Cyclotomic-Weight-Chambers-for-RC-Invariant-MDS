"""Auxiliary finite controls preserved from exploratory branches."""
from exact import *
from itertools import combinations


def auxiliary():
    old=read_set(ROOT/'inputs/S8234.words.txt',8234)
    exceptional=set().union(*(set(p)&old for p in pcrs() if len(set(p)&old)>1))
    orbits=[];seen=set()
    for w in sorted(exceptional):
        if w in seen:continue
        o=frozenset((w,rc(w)));orbits.append(o);seen.update(o)
    expected=set()
    for mask in range(1,1<<len(orbits)):
        removed=set().union(*(o for j,o in enumerate(orbits) if mask>>j&1))
        if all(set(p)&(old-removed) for p in pcrs()):expected.add(tuple(sorted(removed)))
    rows=load(ROOT/'experiments/auxiliary/subset_scan.json')
    need({tuple(row['removed']) for row in rows}==expected and len(rows)==len(expected)==8,'OLD_SUBSET_COVERAGE')
    cyclic=too_long=0
    for row in rows:
        selected=old-set(row['removed'])
        need(len(selected)==row['cardinality'],'OLD_SUBSET_CARDINALITY')
        if 'error' in row:
            try:graph(selected)
            except ValueError as exc:need(str(exc)=='RESIDUAL_DIRECTED_CYCLE','OLD_SUBSET_ERROR');cyclic+=1
            else:raise ValueError('OLD_SUBSET_EXPECTED_CYCLE')
        else:
            d,b,p=graph(selected);walk(row['path'],simple=True)
            need(len(row['path'])==row['L']==len(p)==107,'OLD_SUBSET_PATH_LENGTH')
            need(not set(row['path'])&selected,'OLD_SUBSET_PATH_RESIDUAL');too_long+=1
    base=read_set(ROOT/'inputs/M107.words.txt',8230)
    fixed={w for w in WORDS if rc(w)==w};D_deleted=base|fixed
    data=load(ROOT/'experiments/auxiliary/stable_extension1.json')
    path=[WORDS[e] for e in data['input_path']];walk(path,simple=True)
    need(len(path)==107 and not set(path)&D_deleted,'AUX_STABLE_BASE_PATH')
    candidates={min(w,rc(w)) for w in path}
    records=data['candidates'];need(len(records)==len(candidates)==105,'AUX_STABLE_CANDIDATE_COUNT')
    need({WORDS[x['orbit']] for x in records}==candidates,'AUX_STABLE_CANDIDATE_COVERAGE')
    for row in records:
        e=WORDS[row['orbit']];p=[WORDS[x] for x in row['path']];walk(p,simple=True)
        need(len(p)==107 and row['LD']==107,'AUX_STABLE_PATH_LENGTH')
        need(not set(p)&(D_deleted|{e,rc(e)}),'AUX_STABLE_PATH_RESIDUAL')
    need(data['success_count']==0,'AUX_STABLE_NO_SUCCESS')
    return {'old_8234_proper_RC_subsets_meeting_all_PCRs':8,
            'cyclic':cyclic,'acyclic_but_L107':too_long,
            'fixed_nonstable_M107_plus_one_RC_pair':{'kernel_size':105,'all_rejected_by_nonfixed_paths':True},
            'scope':'Only these frozen-base subproblems. Not new global lower bounds.'}
