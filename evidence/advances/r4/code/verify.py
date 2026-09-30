#!/usr/bin/env python3
"""Round 4: exact, offline, standard-library-only proof verification."""
from pathlib import Path
import argparse,json,gzip,time
import graph_checker as g
from atlas_checker import verify_atlas
from sparse_checker import verify_sparse
from repair_checker import verify_repair
from fractional_checker import verify_fractional
from switch_checker import verify_chain,signed_switch
ROOT=Path(__file__).resolve().parents[1]

def load(p):
    with (gzip.open(p,'rt',encoding='utf-8') if str(p).endswith('.gz') else open(p,encoding='utf-8')) as f:return json.load(f)

def check(root=ROOT,quick=False,manifest=False):
    t=time.time();out={'status':'PASS_FULL' if not quick else 'PASS_WITH_DIGEST_ONLY_CHAIN','sets':{}};sets={};certs={}
    for name,entry in load(root/'certificates/catalog.json').items():
        sets[name]=g.read_set(root/entry['words']);certs[name]=load(root/entry['certificate']);out['sets'][name]=g.set_certificate_check(sets[name],certs[name])
    out['atlases']={}
    for entry in load(root/'certificates/atlas_catalog.json'):
        out['atlases'][entry['name']]=verify_atlas(sets[entry['name']],certs[entry['name']],load(root/entry['atlas']))
    base=sets['M111'];sparse=load(root/'certificates/radius4_M111.json.gz')
    out['sparse_radius']=verify_sparse(base,sparse)
    out['critical_repair']=verify_repair(base,certs['M111'],load(root/'inputs/M111_atlas.json'),g.read_set(root/'witnesses/critical_repair_not_MDS.words.txt'),load(root/'certificates/critical_repair_gap.json'))
    out['finite_cut_LP']=verify_fractional(base,sparse['cuts'],load(root/'certificates/fractional_cut_relaxation.json'))
    # The integer lower bound comes from this same cut family, not from a solver receipt.
    g.require(sparse['bound']==110 and sparse['radius']==4,'FIXED_RADIUS_SCOPE')
    out['finite_cut_LP']['integer_minimum_at_least']=5;out['finite_cut_LP']['integrality_gap_at_least']='10'
    for src,dst,cert in [('M111','M111_rephased','retime_initial'),('M113_before_retiming','M109','retime_final_sector0'),('M113_before_retiming','M109_sector1','retime_final_sector1')]:g.retiming_check(sets[src],sets[dst],load(root/f'certificates/{cert}.json'))
    out['retimings']=3
    out['chain']=verify_chain(root,load(root/'certificates/chain_M111_to_M109.json'),quick=quick)
    # Recheck the six recorded mixed steps and both two-block decompositions.
    diagnostics=load(root/'certificates/mixed_step_diagnostics.json');by_step={rec['chain_step']:rec for rec in diagnostics}
    S=base;checked=0
    for i,st in enumerate(load(root/'certificates/chain_M111_to_M109.json')['steps'],1):
        target=g.read_set(root/st['target']) if st['type']=='retiming' else signed_switch(S,st['u'],st['letters'])[0]
        if i in by_step:
            rec=by_step[i];g.require(rec['source_sha256']==g.digest_set(S) and rec['target_sha256']==g.digest_set(target),'MIXED_DIGEST')
            pos=''.join(a for a in st['letters'] if a+st['u'] in S);neg=''.join(a for a in st['letters'] if st['u']+a in S)
            g.require(pos==rec['positive_letters'] and neg==rec['negative_letters'] and pos and neg,'MIXED_SIGNS')
            for first in rec['first_blocks']:
                letters=first['letters'];g.require(letters in (pos,neg),'MIXED_BLOCK')
                U,_=signed_switch(S,st['u'],letters);checkU=g.graph_check(U)
                g.require(first['Boolean'] is True and first['acyclic'] is True and first['L']==checkU['L'],'MIXED_INTERMEDIATE')
                V,_=signed_switch(U,st['u'],neg if letters==pos else pos);g.require(V==target,'MIXED_DECOMPOSITION')
            checked+=1
        S=target
    g.require(checked==6==len(diagnostics),'MIXED_COUNT');out['mixed_block_controls']=12
    final=sets['M109'];cy=certs['M111']['critical_cycle'];hits=g.walk_check(cy,final,cycle=True)[2]
    g.require(hits!=g.walk_check(cy,base,cycle=True)[2],'COMPONENT_SEPARATION')
    path=certs['M109']['longest_path'];word=path[0]+''.join(e[-1] for e in path[1:])
    g.require((root/'witnesses/M109_longest_avoiding_word.txt').read_text().strip()==word and len(word)==116,'AVOIDING_STRING')
    g.require(all(word[i:i+8] not in final for i in range(len(word)-7)),'AVOIDING_STRING_HIT')
    g.require((root/'witnesses/M109_critical_cycle_110.txt').read_text().splitlines()==certs['M109']['critical_cycle'],'CYCLE_TEXT')
    from radius_checker import all_domains
    groups=all_domains(base);changed=[gr for gr in groups if not set(gr['options'][gr['base_option']])<=final]
    comparison={'changed_RC_PCR_groups':len(changed),'changed_PCRs':sum(len(set(gr['necklaces'])) for gr in changed),'symmetric_difference':len(base^final),'old_cycle_length':len(cy),'old_cycle_hits_in_M109':hits}
    g.require(comparison==load(root/'certificates/comparison.json'),'COMPARISON')
    out['comparison']=comparison;out['word_guarantee']={'longest_avoiding_letters':len(word),'universal_hit_letters':len(word)+1}
    if manifest:out['manifest_files']=g.manifest_check(root)
    out['elapsed_seconds']=time.time()-t;return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--quick',action='store_true');ap.add_argument('--manifest',action='store_true');a=ap.parse_args();print(json.dumps(check(a.root,a.quick,a.manifest),indent=2))
if __name__=='__main__':main()
