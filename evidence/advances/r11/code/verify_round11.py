#!/usr/bin/env python3
"""Read finite certificates, reconstruct their domains, and check exact conclusions."""
from exact import *
import argparse


def main(out):
    base = read_set(ROOT/'inputs/M107.words.txt', 8230)
    old = read_set(ROOT/'inputs/S8234.words.txt', 8234)
    new = read_set(ROOT/'certificates/S8232.words.txt', 8232)
    need(len(graph(base)[2]) == 107, 'BASE_L107')
    height(load(ROOT/'inputs/S8234.height.json'), old)
    _, finite = height(load(ROOT/'certificates/S8232.height.json'), new)
    need(finite == {'cardinality': 8232, 'L': 105, 'surplus': 2,
                   'PCR_cut_histogram': {'1': 8228, '2': 2}}, 'NEW_COUNTS')
    sp = surplus_profile(new)
    need(sp['paired_excess'] == 1 and sp['stable_excess'] == 0 and sp['selected_RC_fixed_words'] == 128, 'NEW_SURPLUS_PROFILE')
    rp = rank_path(load(ROOT/'certificates/S8232.rank_path.json'), new)
    ir = irredundancy(load(ROOT/'certificates/irredundancy.json'), new)
    patch = load(ROOT/'certificates/patch_and_path.json')
    for key, actual in [('removed_from_8234', old-new), ('added_to_8234', new-old),
                        ('removed_from_M107', base-new), ('added_to_M107', new-base)]:
        need(patch.get(key) == sorted(actual), 'PATCH_'+key)
    need(patch.get('cardinality') == 8232 and patch.get('L') == 105 and patch.get('surplus') == 2,
         'PATCH_VALUES')
    walk(patch['path'], simple=True)
    need(len(patch['path']) == 105 and not set(patch['path']) & new, 'PATCH_PATH')
    changed = [g for g in groups() if (g['words'] & new) != (g['words'] & old)]
    stats = load(ROOT/'certificates/patch_statistics.json')
    oh = load(ROOT/'inputs/S8234.height.json')['heights']
    nh = load(ROOT/'certificates/S8232.height.json')['heights']
    actual = {'removed': len(old-new), 'added': len(new-old),
              'changed_PCRs': sum((set(p)&old) != (set(p)&new) for p in pcrs()),
              'changed_nonstable_groups': sum(not g['stable'] for g in changed),
              'changed_stable_groups': sum(g['stable'] for g in changed),
              'changed_group_representatives': [g['representative'] for g in changed],
              'changed_height_pairs': sum(oh[VI[v]] != nh[VI[v]] for v in VR),
              'max_height_change': max(abs(a-b) for a, b in zip(oh, nh))}
    need(stats == actual, 'PATCH_STATISTICS')
    blocks = load(ROOT/'experiments/pcr_phase1/log.json')
    bm = block_summary(blocks[2], old, new)
    result = {'status': 'PASS', 'base': {'cardinality': 8230, 'L': 107},
              'new_height': finite, 'surplus_profile': sp, 'rank_and_path': rp, 'irredundancy': ir,
              'patch': actual, 'positive_block': bm,
              'claim_boundary': '8232 is not an MDS. No global lower bound above 8230.'}
    save(out/'verification.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True); main(args.out)
