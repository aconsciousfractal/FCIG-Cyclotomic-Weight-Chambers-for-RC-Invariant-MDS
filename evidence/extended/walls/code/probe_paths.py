"""Bounded search along the stored paired geodesics in DNA k=8 component 14.

This does NOT enumerate the entire F-component and does NOT establish a lower
bound for all RC-MDS.  Every reached state is checked as a literal RC set and
its exact longest path is certified by integer levels and a matching path.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from core import check, dag, dump, inedges, locate_parent, outedges, pack, parent_data, rc


def probe(parent: Path, results: Path) -> dict:
    atlas, terminals = parent_data(parent, 8)
    by_id = {row['id']: i for i, row in enumerate(atlas['rows'])}
    census = json.loads((results / 'components_k8.json').read_text(encoding='utf-8'))
    representative = census['classes'][14][0]
    paths = [p for p in census['positive_geodesics'] if p['source'] == representative]
    n = 4**7
    expected_size = sum(terminals[by_id[representative]])
    seen: set[bytes] = set()
    states = []
    for record in paths:
        selected = terminals[by_id[representative]][:]
        for step, u in enumerate(record['moves'], start=1):
            v = rc(u, 7)
            iu, ou = set(inedges(u, n)), set(outedges(u, n))
            iv, ov = set(inedges(v, n)), set(outedges(v, n))
            check(not ((iu | ou) & (iv | ov)), 'PROBE_STAR_OVERLAP')
            check(all(selected[e] for e in iu | ov), 'PROBE_ILLEGAL_MOVE')
            for e in (iu - ou) | (ov - iv):
                selected[e] = 0
            for e in (ou - iu) | (iv - ov):
                selected[e] = 1
            key = bytes(selected)
            if key in seen:
                continue
            seen.add(key)
            check(sum(selected) == expected_size, 'PROBE_CARDINALITY')
            check(all(selected[e] == selected[rc(e, 8)] for e in range(4*n)), 'PROBE_RC')
            _, length, path, levels = dag(selected, n)
            check(len(path) == length + 1, 'PROBE_WITNESS_LENGTH')
            check(all(not selected[u*4 + v % 4] and (u*4 + v % 4) % n == v
                      for u, v in zip(path, path[1:])), 'PROBE_BAD_PATH')
            check(all(selected[e] or levels[e % n] > levels[e // 4]
                      for e in range(4*n)), 'PROBE_BAD_LEVELS')
            states.append({'target': record['target'], 'prefix_steps': step,
                           'terminal_sha256': hashlib.sha256(pack(selected)).hexdigest(),
                           'L_lower_edges': length, 'path': path})
    histogram: dict[str, int] = {}
    for state in states:
        key = str(state['L_lower_edges'])
        histogram[key] = histogram.get(key, 0) + 1
    return {'status': 'PASS', 'claim': 'CLN: bounded path sample only', 'k': 8,
            'component_id': 14, 'representative': representative,
            'selection_rule': 'All nonempty prefixes of the stored root-to-row geodesics in class 14; deduplicate literal terminal bit vectors.',
            'distinct_states': len(states), 'histogram': histogram, 'states': states,
            'not_claimed': 'No enumeration of the full component; no global optimality claim; no claim that intermediates have weight realizations.'}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'results')
    args = parser.parse_args()
    parent, temporary = locate_parent(args.parent)
    try:
        result = probe(parent, args.out)
        dump(args.out / 'bounded_path_probe.json', result)
        print('BOUNDED_PROBE_PASS', result['distinct_states'], result['histogram'], flush=True)
    finally:
        if temporary:
            temporary.cleanup()


if __name__ == '__main__':
    main()
