"""Replay accepted experimental endpoints from sparse patches and saved hashes.
Restart sources are recovered by exact equality, not guessed. Resets are allowed;
only the accepted endpoints (not a connected move chain) were recorded upstream.
"""
from verify_round9 import *
from collections import Counter

def update(source,h):
    fixed={w for w in EDGES if rc(w)==w};m=(source-fixed)|fixed;delta=Counter();u=V[h['u']]
    for a in range(4):
        if h['mask']>>a&1:
            for w,s in [(u+ABC[a],1),(ABC[a]+u,-1)]:delta[w]+=s;delta[rc(w)]+=s
    for w,s in delta.items():
        x=int(w in m)+s
        if x not in (0,1):return None
        if x:m.add(w)
        else:m.discard(w)
    restored={decode(i) for i in h['restored_fixed_edges']}
    if len(restored)!=128 or not restored<=fixed:return None
    m=(m-fixed)|(fixed-restored)
    digest=hashlib.sha256(bytes(int(w in m) for w in EDGES)).hexdigest()
    return m if digest==h['sha256_indicator'] else None

def main():
    start=time.time();base=read_words(ROOT/'inputs/M107.words.txt');m=base;hist=json.loads((ROOT/'logs/stable_walk_history.json').read_text());ids=set();resets=[];profile=Counter();stable_changes=Counter();fixed={w for w in EDGES if rc(w)==w}
    for i,h in enumerate(hist):
        out=update(m,h);source='previous'
        if out is None:out=update(base,h);source='base'
        req(out is not None,'RECOVER_ACCEPTED_ENDPOINT_'+str(i))
        if source=='base':resets.append(i)
        req(len(out)==8230 and {rc(w) for w in out}==out,'ACCEPTED_MEMBERSHIP')
        # Full necklace census and full longest-path traversal, not just hashes.
        necklaces(out);_,_,d,_,path=graph_data(EDGESET-out)
        req(len(path)==h['L'] and len(path)<=114,'ACCEPTED_DAG_LENGTH')
        req(h['sha256_indicator'] not in ids,'REPEATED_ACCEPTED_ENDPOINT');ids.add(h['sha256_indicator'])
        profile[len(path)]+=1;stable_changes[len((out^base)&fixed)//2]+=1;m=out
        if i%200==0:print('checked',i+1,'elapsed',round(time.time()-start,2),flush=True)
    req(m==read_words(ROOT/'certificates/stable_walk_final.words.txt'),'FINAL_WALK_ENDPOINT')
    out={'status':'PASS','accepted_endpoints':len(hist),'distinct_hashes':len(ids),'recovered_reset_indices':resets,'length_profile':dict(sorted(profile.items())),'stable_changes_from_base_profile':dict(sorted(stable_changes.items())),'full_PCR_and_DAG_checks':len(hist),'seconds':time.time()-start,'boundary':'replay of accepted endpoints, not an uninterrupted legal move chain; rejected frontier candidates are not independently replayed'}
    (OUT/'walk_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
