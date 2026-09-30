"""Portable certificate checker for Round 4.
No imports from any previous round, search executable, or retiming producer
are used in checking the 72 MDS/cycle/patch certificates. small_controls
separately tests the delivered retiming algorithm against exhaustive spaces.
"""
from __future__ import annotations
from pathlib import Path
from collections import deque
from functools import lru_cache
from copy import deepcopy
import argparse,hashlib,json,zipfile,io
from itertools import combinations

class Invalid(ValueError): pass

def need(ok,code):
    if not ok: raise Invalid(code)

def js(path): return json.loads(path.read_text())
def digest(b): return hashlib.sha256(b).hexdigest()

def rc(x,k,q):
    digits=[]
    for _ in range(k): digits.append(x%q);x//=q
    y=0
    for a in digits: y=q*y+(q-1-a)
    return y

@lru_cache(None)
def carrier(q,k):
    n=q**(k-1);E=q*n;visited=bytearray(E);cycles=[]
    for e in range(E):
        if visited[e]:continue
        orbit=[];x=e
        while not visited[x]:
            visited[x]=1;orbit.append(x)
            # RIGHT rotation; independent of the producer's left rotation.
            x=(x%q)*n+x//q
        need(x==e,'PCR_PARTITION');cycles.append(orbit)
    return n,E,cycles,[rc(e,k,q) for e in range(E)],[rc(v,k-1,q) for v in range(n)]

def unpack(path,size):
    raw=path.read_bytes();need(len(raw)==(size+7)//8,'BIT_LENGTH')
    need(size%8==0 or raw[-1]>>(size%8)==0,'BIT_PADDING')
    return [int(bool(raw[i>>3] & (1<<(i&7)))) for i in range(size)]

def check_indicator(m,q,k):
    n,E,cycles,theta,_=carrier(q,k)
    need(len(m)==E and all(type(x) is int and x in (0,1) for x in m),'INDICATOR_TYPE')
    for C in cycles: need(sum(m[e] for e in C)==1,'PCR_HIT')
    need(all(m[e]==m[theta[e]] for e in range(E)),'RC')
    return len(cycles)

def check_path(m,q,k,path,L):
    n,E,*_=carrier(q,k)
    need(type(L) is int and L>=0 and isinstance(path,list) and len(path)==L,'PATH_LENGTH')
    need(all(type(e) is int and 0<=e<E for e in path),'PATH_RANGE')
    need(all(not m[e] for e in path),'PATH_USES_DELETION')
    need(all(a%n==b//q for a,b in zip(path,path[1:])),'PATH_ADJACENCY')
    vertices=[e//q for e in path]+([path[-1]%n] if path else [])
    need(len(vertices)==len(set(vertices)),'PATH_REPEATS_VERTEX')

def check_rank(m,q,k,rank,L):
    n,E,*_=carrier(q,k)
    need(len(rank)==n and all(type(x) is int and 0<=x<=L for x in rank),'RANK_RANGE')
    need(all(m[e] or rank[e%n]>rank[e//q] for e in range(E)),'RANK_EDGE')

def reverse_dp(m,q,k):
    """Independent sink elimination, opposite to the producer's source DP."""
    n,E,*_=carrier(q,k);out=[0]*n;height=[0]*n
    for e,x in enumerate(m):
        if not x:out[e//q]+=1
    stack=[u for u in range(n) if out[u]==0];seen=0
    while stack:
        v=stack.pop();seen+=1
        for a in range(q):
            e=a*n+v
            if m[e]:continue
            u=e//q;height[u]=max(height[u],height[v]+1);out[u]-=1
            if out[u]==0:stack.append(u)
    need(seen==n,'REVERSE_DP_CYCLE')
    return max(height)

def check_cycle(m,q,k,es,hits=None,L=None):
    n,E,*_=carrier(q,k)
    need(isinstance(es,list) and 1<=len(es)<=n,'CYCLE_LENGTH')
    need(all(type(e) is int and 0<=e<E for e in es),'CYCLE_RANGE')
    need(len({e//q for e in es})==len(es),'CYCLE_NOT_SIMPLE')
    need(all(es[i]%n==es[(i+1)%len(es)]//q for i in range(len(es))),'CYCLE_ADJACENCY')
    actual=sum(m[e] for e in es)
    if hits is not None:need(actual==hits,'CYCLE_HITS')
    if L is not None:need(hits==1 and len(es)==L+1,'CYCLE_LOWER_BOUND')
    return actual

def check_retiming(source,target,q,k,T,epsilon,H,f):
    n,E,_,_,theta=carrier(q,k)
    need(type(T) is int and T>=1 and epsilon in (0,1),'RETIMING_DOMAIN')
    need(len(H)==len(f)==n and all(type(x) is int for x in H+f),'RETIMING_INTEGER')
    need(all(H[v]+H[theta[v]]==T-1-T*epsilon for v in range(n)),'RETIMING_RC_OFFSET')
    need(all(H[e%n]-H[e//q]>=1-T*source[e] for e in range(E)),'RETIMING_INEQUALITY')
    need(all(f[v]==-(H[v]//T) and f[v]+f[theta[v]]==epsilon for v in range(n)),'RETIMING_FLOOR')
    need(all(target[e]==source[e]+f[e//q]-f[e%n] for e in range(E)),'RETIMING_GRADIENT')

def check_patch(source,target,rec,q=4,k=8):
    n,E,*_=carrier(q,k);u=rec['lower_vertex'];B=rec['alphabet_mask']
    need(type(u) is int and 0<=u<n and type(B) is int and 0<B<(1<<q)-1,'PATCH_DOMAIN')
    delta={}
    for a in range(q):
        if not B>>a&1:continue
        for e,sign in [(a*n+u,-1),(q*u+a,1)]:
            delta[e]=delta.get(e,0)+sign
            mate=rc(e,k,q);delta[mate]=delta.get(mate,0)+sign
    expected=source[:]
    for e,d in delta.items():expected[e]+=d
    need(all(x in (0,1) for x in expected),'PATCH_NONBOOLEAN')
    need(target==expected,'PATCH_ENDPOINT')
    removed=[e for e in range(E) if source[e] and not target[e]]
    added=[e for e in range(E) if target[e] and not source[e]]
    need(removed==rec['removed'] and added==rec['added'],'PATCH_SUPPORT')
    need(len(removed)==len(added),'PATCH_CARDINALITY')
    def word(e,length):
        t=[]
        for _ in range(length):t.append('ACGT'[e%4]);e//=4
        return ''.join(reversed(t))
    need([word(e,k) for e in removed]==rec['removed_words'],'PATCH_REMOVED_WORDS')
    need([word(e,k) for e in added]==rec['added_words'],'PATCH_ADDED_WORDS')
    need(word(u,k-1)==rec['lower_word'],'PATCH_LOWER_WORD')

def check_separator(a,b,rec,q=4,k=8):
    C=rec['edge_codes'];need(rec['length']==len(C)<=2*k-1,'SEPARATOR_LENGTH')
    x=check_cycle(a,q,k,C);y=check_cycle(b,q,k,C)
    need(x==rec['hits_source'] and y==rec['hits_target'] and x!=y,'SEPARATOR_HITS')

def check_manifest(root):
    manifest=root/'MANIFEST_SHA256.txt'
    need(manifest.is_file(),'MANIFEST_MISSING')
    listed=set()
    for line in manifest.read_text().splitlines():
        h,p=line.split('  ',1);need(p not in listed,'MANIFEST_DUPLICATE');listed.add(p)
        f=root/p;need(f.is_file() and digest(f.read_bytes())==h,'MANIFEST_HASH:'+p)
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts}
    need(actual==listed,'MANIFEST_FILE_SET')
    return {'status':'PASS','files':len(listed)}

def check_optimum_coverage(rows, certificates):
    for row in rows:
        k,cid=row['k'],row['component']
        expected=([f'k8_component_{cid:02d}_sector_{ep}' for ep in (0,1)]
                  if k==8 else [row['source']])
        need(row['attainers']==expected,'OPTIMUM_ATTAINER_COVERAGE')
        for ep,name in enumerate(expected):
            if k==8:
                r=certificates[name]
                need(r['source']==row['source'] and r['epsilon']==ep,
                     'OPTIMUM_RETIMING_JOIN')

def check_chain_coverage(chain):
    need([r['name'] for r in chain['component_optimality']]==chain['states'],
         'ESCAPE_OPTIMALITY_COVERAGE')
    need(all(r['cycle_length']==len(r['cycle_edges']) for r in chain['component_optimality']),
         'ESCAPE_CYCLE_RECORDED_LENGTH')

def check_separator_coverage(seps,ps,optima,states):
    rows={r['component']:r['source'] for r in optima if r['k']==8}
    need(len(seps)==21 and {r['source_component'] for r in seps}==set(rows),
         'SEPARATOR_COVERAGE')
    need(all(r['source']==rows[r['source_component']] and r['target']=='escape_131' for r in seps),
         'SEPARATOR_SOURCE_JOIN')
    need([(r['source'],r['target']) for r in ps]==list(combinations(states,2)),
         'ESCAPE_PAIR_COVERAGE')

def inherited_custody(root, bases):
    parent=root.parent/'atlas'
    for row in bases:
        expected=parent/'results/certificates'/f"{row['row']}.bits"
        need(expected.read_bytes()==(root/row['input_path']).read_bytes(),
                f"source bytes {row['row']}")
    return {'base_byte_joins':len(bases),'source':'../atlas/results/certificates'}

def run(root,include_small=True):
    certificates={};indicators={};edge_checks=0;row_hashes=[]
    for path in sorted((root/'certificates').glob('*.json')):
        rec=js(path);name=rec['name'];need(name==path.stem and name not in certificates,'CERT_NAME')
        q,k=rec['q'],rec['k'];n,E,*_=carrier(q,k);raw=(root/rec['bits']).read_bytes()
        need(digest(raw)==rec['bits_sha256'],'CERT_BITS_HASH');m=unpack(root/rec['bits'],E)
        N=check_indicator(m,q,k);need(N==rec['size']==sum(m),'MDS_SIZE')
        check_rank(m,q,k,rec['levels'],rec['L']);check_path(m,q,k,rec['path_edges'],rec['L'])
        need(reverse_dp(m,q,k)==rec['L'],'REVERSE_DP_LENGTH')
        certificates[name]=rec;indicators[name]=m;edge_checks+=E;row_hashes.append([name,digest(raw)])
    need(len(certificates)==72,'CERTIFICATE_COVERAGE')
    bases=js(root/'inputs/source_rows.json');need(len(bases)==25,'SOURCE_ROW_COVERAGE')
    for row in bases:
        b=(root/row['input_path']).read_bytes();need(digest(b)==row['input_sha256'],'SOURCE_BYTES')
        need(unpack(root/row['input_path'],4**row['k'])==indicators[row['row']],'SOURCE_JOIN')
    custody=inherited_custody(root,bases)
    retimings=0
    for name,r in certificates.items():
        if 'cut_potential' in r:
            check_retiming(indicators[r['source']],indicators[name],r['q'],r['k'],r['T'],r['epsilon'],r['H'],r['cut_potential']);retimings+=1
    need(retimings==42,'RETIMING_COVERAGE')
    optima=js(root/'results/component_optima.json');need(len(optima)==25,'OPTIMUM_COVERAGE')
    need([(r['k'],r['component']) for r in optima]==[(4,0)]+[(6,i) for i in range(3)]+[(8,i) for i in range(21)],'COMPONENT_INDEX')
    check_optimum_coverage(optima,certificates)
    for row in optima:
        k=row['k'];C=row['cycle_edges'];L=row['optimum'];src=indicators[row['source']]
        need(row['cycle_length']==len(C),'CYCLE_RECORDED_LENGTH');check_cycle(src,4,k,C,row['cycle_hits'],L)
        for name in row['attainers']:
            need(certificates[name]['L']==L,'OPTIMUM_ATTAINMENT');check_cycle(indicators[name],4,k,C,1,L)
    chain=js(root/'results/escape_chain.json');need(chain['states']==['start_141','escape_135','prepare_135','escape_133','escape_131'],'ESCAPE_STATE_COVERAGE')
    for i,r in enumerate(chain['patches']):
        need(r['source']==chain['states'][i] and r['target']==chain['states'][i+1],'PATCH_CHAIN')
        check_patch(indicators[r['source']],indicators[r['target']],r)
        need(r['L_before']==certificates[r['source']]['L'] and r['L_after']==certificates[r['target']]['L'],'PATCH_RECORDED_LENGTH')
    need(len(chain['patches'])==4,'PATCH_COUNT')
    check_chain_coverage(chain)
    for r in chain['component_optimality']:
        need(r['optimum']==certificates[r['name']]['L'],'ESCAPE_OPTIMUM_JOIN')
        check_cycle(indicators[r['name']],4,8,r['cycle_edges'],r['cycle_hits'],r['optimum'])
    a=indicators['start_141'];b=indicators['escape_131']
    need(sum(x!=y for x,y in zip(a,b))==chain['changed_from_start']==20,'ESCAPE_HAMMING')
    need(a==indicators['k8_c094'],'ESCAPE_SOURCE_JOIN')
    def dna(e,k=8):
        digits=[]
        for _ in range(k):digits.append('ACGT'[e%4]);e//=4
        return ''.join(digits[::-1])
    need(chain['final_removed_words']==[dna(e) for e in range(len(a)) if a[e] and not b[e]],'FINAL_REMOVALS_TEXT')
    need(chain['final_added_words']==[dna(e) for e in range(len(a)) if b[e] and not a[e]],'FINAL_ADDITIONS_TEXT')
    need((root/'certificates/escape_131.words.txt').read_text().splitlines()==[dna(e) for e in range(len(b)) if b[e]],'FINAL_WORD_LIST')
    compact=js(root/'results/cycle_132_compact.json');cyclic=compact['cyclic_word']
    need(compact['length']==len(cyclic)==132,'COMPACT_CYCLE_LENGTH')
    need(compact['edge_codes']==chain['component_optimality'][-1]['cycle_edges'],'COMPACT_CYCLE_JOIN')
    need(all(''.join(cyclic[(i+j)%132] for j in range(8))==dna(e) for i,e in enumerate(compact['edge_codes'])),'COMPACT_CYCLE_WINDOWS')
    hits=[i for i,e in enumerate(compact['edge_codes']) if b[e]]
    need(hits==[compact['unique_selected_index']] and dna(compact['edge_codes'][hits[0]])==compact['unique_selected_word'],'COMPACT_CYCLE_HIT')

    seps=js(root/'results/separators_131.json');need({r['source_component'] for r in seps}==set(range(21)) and len(seps)==21,'SEPARATOR_COVERAGE')
    for r in seps:check_separator(indicators[r['source']],indicators[r['target']],r)
    ps=js(root/'results/escape_component_separators.json');need(len(ps)==10,'ESCAPE_SEPARATORS_COVERAGE')
    check_separator_coverage(seps,ps,optima,chain['states'])
    for r in ps:check_separator(indicators[r['source']],indicators[r['target']],r)
    # Mutations target semantics directly, rather than only digest mismatch.
    mutations=[]
    def reject(tag,fn):
        try:fn()
        except (Invalid,ValueError,IndexError,KeyError,TypeError) as exc:mutations.append({'tag':tag,'rejected_by':str(exc)});return
        raise Invalid('MUTATION_ACCEPTED:'+tag)
    r=certificates['escape_131'];m=indicators['escape_131'];n,E,*_=carrier(4,8)
    mm=m[:];mm[0]=0;reject('mandatory_loop_removed',lambda:check_indicator(mm,4,8))
    mm2=m[:];e=next(e for e in range(E) if not m[e]);mm2[e]=1;reject('one_word_added',lambda:check_indicator(mm2,4,8))
    mm3=m[:];e=next(e for e in range(E) if m[e] and e!=rc(e,8,4));mm3[e]=0;reject('rc_mate_not_changed',lambda:check_indicator(mm3,4,8))
    reject('nonbinary_indicator',lambda:check_indicator([2]+m[1:],4,8))
    reject('fractional_indicator',lambda:check_indicator([1.0]+m[1:],4,8))
    reject('flat_rank',lambda:check_rank(m,4,8,[0]*n,r['L']))
    rr=r['levels'][:];rr[0]=-1;reject('negative_rank',lambda:check_rank(m,4,8,rr,r['L']))
    reject('length_lowered_without_new_certificate',lambda:check_rank(m,4,8,r['levels'],130))
    reject('path_truncated',lambda:check_path(m,4,8,r['path_edges'][:-1],r['L']))
    path=r['path_edges'][:];path[0]=0;reject('path_uses_removed_loop',lambda:check_path(m,4,8,path,r['L']))
    cy=chain['component_optimality'][-1]['cycle_edges'];reject('cycle_wrong_hit_count',lambda:check_cycle(m,4,8,cy,2))
    broken=cy[:];broken[0]=0;reject('cycle_broken',lambda:check_cycle(m,4,8,broken,1))
    reject('cycle_wrong_length_bound',lambda:check_cycle(m,4,8,cy,1,130))
    z=certificates['k8_component_07_sector_0'];src=indicators[z['source']];target=indicators[z['name']]
    reject('wrong_rc_sector',lambda:check_retiming(src,target,4,8,z['T'],1,z['H'],z['cut_potential']))
    HH=z['H'][:];HH[0]+=1;reject('one_potential_coordinate',lambda:check_retiming(src,target,4,8,z['T'],0,HH,z['cut_potential']))
    ff=z['cut_potential'][:];ff[0]+=1;reject('cut_floor_wrong',lambda:check_retiming(src,target,4,8,z['T'],0,z['H'],ff))
    patch=deepcopy(chain['patches'][0]);patch['alphabet_mask']=1;reject('patch_wrong_letter',lambda:check_patch(indicators[patch['source']],indicators[patch['target']],patch))
    sep=deepcopy(seps[0]);sep['hits_target']=sep['hits_source'];reject('separator_counts_fused',lambda:check_separator(indicators[sep['source']],m,sep))
    missing=deepcopy(chain);missing['component_optimality'].pop()
    reject('missing_final_component_lower',lambda:check_chain_coverage(missing))
    missing=deepcopy(optima);missing[-1]['attainers'].pop()
    reject('missing_opposite_sector',lambda:check_optimum_coverage(missing,certificates))
    doubled=deepcopy(optima);doubled[-1]['attainers'][1]=doubled[-1]['attainers'][0]
    reject('duplicate_sector_attainer',lambda:check_optimum_coverage(doubled,certificates))
    wrong=deepcopy(seps);wrong[0]['source']=wrong[1]['source']
    reject('separator_wrong_source_join',lambda:check_separator_coverage(wrong,ps,optima,chain['states']))
    doubled=deepcopy(ps);doubled[-1]=doubled[0]
    reject('duplicate_escape_pair',lambda:check_separator_coverage(seps,doubled,optima,chain['states']))
    bad=deepcopy(chain);bad['component_optimality'][-1]['cycle_length']-=1
    reject('cycle_metadata_off_by_one',lambda:check_chain_coverage(bad))
    small=None
    if include_small:
        import small_controls
        small=small_controls.run();need(small==js(root/'results/small_controls.json'),'SMALL_CONTROL_REPLAY')
    return {'status':'PASS','inherited_custody':custody,'mds_certificates':len(certificates),'independent_indicator_edge_checks':edge_checks,'retiming_certificates':retimings,'old_components_optimized':25,'k8_ordinary_components':21,'k8_paired_components':42,'new_chain_L':[certificates[x]['L'] for x in chain['states']],'new_final_size':sum(m),'new_final_component_minimum':131,'literal_change_from_start':20,'long_cycle_lower_certificates':30,'short_separators':31,'semantic_mutations':mutations,'small_controls':small,'indicator_digest':digest(json.dumps(row_hashes,sort_keys=True).encode())}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output',type=Path);p.add_argument('--manifest',action='store_true');a=p.parse_args()
    if a.manifest:r=check_manifest(a.root)
    else:r=run(a.root)
    raw=json.dumps(r,sort_keys=True,indent=2)+'\n'
    if a.output:a.output.write_text(raw)
    print(raw)
if __name__=='__main__':main()
