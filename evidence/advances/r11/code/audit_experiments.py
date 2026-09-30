#!/usr/bin/env python3
"""Replay saved endpoints and joins; never treat a solver status as a proof."""
from exact import *
import argparse
import struct
from audit_auxiliary import auxiliary


def fusion_row(rec, input_c, proposal_c, output_c):
    si, _ = height(input_c, full_graph=False)
    sg, _ = height(proposal_c, full_graph=False)
    so, info = height(output_c)
    h, g, z = input_c['heights'], proposal_c['heights'], output_c['heights']
    need(all(z[VI[v]] in (h[VI[v]], g[VI[v]]) for v in VR), 'FUSION_PAIR_CHOICE')
    need(rec.get('before') == len(si) and rec.get('proposal_cost') == len(sg) and
         rec.get('after') == len(so), 'FUSION_COST_JOIN')
    need(len(so) <= len(si), 'FUSION_NONINCREASE')
    need(rec.get('changed_pairs') == sum(h[VI[v]] != z[VI[v]] for v in VR),
         'FUSION_CHANGED_PAIRS')
    need(rec.get('height_variables') == sum(h[VI[v]] != g[VI[v]] for v in VR),
         'FUSION_VARIABLE_COUNT')
    digest = hashlib.sha256(struct.pack('<'+'h'*len(z), *z)).hexdigest()
    need(rec.get('solution_sha256') == digest, 'FUSION_HASH')
    return so, info


def main(out):
    old_c = load(ROOT/'inputs/S8234.height.json')
    old, _ = height(old_c)
    new_c = load(ROOT/'certificates/S8232.height.json')
    new, _ = height(new_c)
    records = load(ROOT/'experiments/fusion/log.json')
    need(len(records) == 24, 'FUSION_RUN_SIZE')
    previous = old_c['heights']
    maps, sets, statuses = set(), set(), Counter()
    fusion_results = []
    for i, row in enumerate(records):
        p = ROOT/'experiments/fusion'
        ic = load(p/f'input_{i:03d}.json')
        pc = load(p/f'proposal_{i:03d}.json')
        oc = load(p/f'output_{i:03d}.json')
        need(ic['heights'] == previous, 'FUSION_INPUT_CHAIN')
        so, info = fusion_row(row, ic, pc, oc)
        need(len(so) == 8234, 'FUSION_RETAINED_CARDINALITY')
        previous = oc['heights']
        maps.add(tuple(previous)); sets.add(word_digest(so)); statuses[row['status']] += 1
        fusion_results.append({'iteration': i, 'verified_cardinality': len(so),
                               'verified_L': info['L'], 'set_sha256': word_digest(so)})
    need(len(maps) == 4 and len(sets) == 1, 'FUSION_DISTINCT_COUNTS')
    need(load(ROOT/'experiments/fusion/latest.json')['heights'] == previous, 'FUSION_LATEST')
    phases = []
    for phase, first_c, count in [('pcr_phase1', old_c, 3), ('pcr_phase2', new_c, 7)]:
        p = ROOT/'experiments'/phase
        rows = load(p/'log.json'); need(len(rows) == count, 'PCR_PHASE_SIZE')
        before, _ = height(first_c)
        results, phase_sets, phase_maps = [], set(), set()
        for i, row in enumerate(rows):
            oc = load(p/f'output_{i:03d}.json'); after, info = height(oc)
            result = block_summary(row, before, after)
            result.update({'iteration': i, 'L': info['L'], 'set_sha256': word_digest(after),
                           'solver_status_metadata_only': row['status']})
            results.append(result); phase_sets.add(word_digest(after));phase_maps.add(tuple(oc['heights']))
            before = after
        need(load(p/'latest.json')['heights'] == oc['heights'], 'PCR_LATEST')
        if phase == 'pcr_phase1':
            need(before == new and [x['after'] for x in results] == [8234,8234,8232], 'POSITIVE_PHASE_JOIN')
        else:
            need(before == new and all(x['after'] == 8232 for x in results), 'NO_FURTHER_WITNESS')
        phases.append({'name': phase, 'calls': count, 'distinct_output_maps': len(phase_maps),
                       'distinct_output_sets': len(phase_sets), 'verified_endpoints': results})
    result = {'status': 'PASS', 'auxiliary': auxiliary(), 'fusion': {'calls': len(records), 'distinct_output_maps': len(maps),
              'distinct_output_sets': len(sets), 'solver_status_counts_metadata_only': dict(statuses),
              'endpoints': fusion_results}, 'pcr_phases': phases,
              'scope': '34 saved post-solve endpoints checked. Internal branch-and-bound states are not replayed; no solver optimality status is a negative certificate.'}
    save(out/'experiment_audit.json', result)
    print('PASS: fusion 24 calls / 4 maps / 1 set; PCR phases 3 + 7 endpoints; positive cost 8232.')

if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);main(a.out)
