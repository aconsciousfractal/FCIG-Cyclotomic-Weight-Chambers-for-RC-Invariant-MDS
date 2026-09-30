"""Literal mixed-direction companion exchange and exact history replay."""
from collections import Counter
import json
import graph_checker as g


def signed_switch(selected,u,letters):
    g.require(isinstance(u,str) and len(u)==7 and set(u)<=set(g.ALPHABET),'SIGNED_VERTEX')
    g.require(isinstance(letters,str) and 0<len(letters)<=4 and len(set(letters))==len(letters)
              and set(letters)<=set(g.ALPHABET),'SIGNED_LETTERS')
    delta=Counter();orientations=[]
    for a in letters:
        incoming,outgoing=a+u,u+a
        sign=int(incoming in selected)-int(outgoing in selected)
        g.require(sign!=0,'SIGNED_INACTIVE_LETTER')
        orientations.append(sign)
        delta[outgoing]+=sign;delta[incoming]-=sign
        delta[g.rc(outgoing)]+=sign;delta[g.rc(incoming)]-=sign
    result=selected.copy()
    g.require(any(delta.values()),'SIGNED_NULL_SWITCH')
    for w,change in delta.items():
        value=int(w in selected)+change;g.require(value in (0,1),'SIGNED_NONBOOLEAN')
        if value:result.add(w)
        else:result.discard(w)
    return result,{'mixed_directions':len(set(orientations))==2,'letters':len(letters)}


def verify_chain(root,chain,quick=False):
    S=g.read_set(root/chain['initial']);g.require(g.digest_set(S)==chain['initial_sha256'],'CHAIN_INITIAL')
    begin=g.graph_check(S)['L'];rows=[];types=Counter();mixed=0;up=down=equal=0;last=begin
    for k,step in enumerate(chain['steps'],1):
        typ=step['type']
        if typ=='partial_patch':target=g.apply_patch(S,step['u'],step['letters'])
        elif typ=='signed_switch':
            target,info=signed_switch(S,step['u'],step['letters']);mixed+=int(info['mixed_directions'])
        elif typ=='retiming':
            target=g.read_set(root/step['target']);c=json.loads((root/step['certificate']).read_text());g.retiming_check(S,target,c)
        else:raise g.VerificationError('CHAIN_UNKNOWN_TYPE')
        g.require(g.digest_set(target)==step['sha256_after'],'CHAIN_DIGEST')
        # Quick mode skips DAG/longest-path recomputation at intermediate steps only.
        # The output labels this explicitly; it is NOT full chain certification.
        L=step['L_after']
        if not quick:g.require(g.graph_check(target)['L']==L,'CHAIN_LENGTH')
        g.require(type(L) is int and L>=0,'CHAIN_L_TYPE')
        up+=L>last;down+=L<last;equal+=L==last;last=L;types[typ]+=1
        rows.append({'step':k,'type':typ,'L':L});S=target
    g.require(S==g.read_set(root/chain['final']),'CHAIN_FINAL')
    final=g.graph_check(S)['L'];g.require(final==last,'CHAIN_FINAL_L')
    return {'status':'PASS_FULL' if not quick else 'PASS_DIGEST_REPLAY_ONLY',
      'verified_steps':len(rows),'types':dict(types),'mixed_direction_steps':mixed,
      'initial_L':begin,'final_L':final,'maximum_L':max([begin]+[s['L'] for s in rows]),
      'strict_increases':up,'strict_decreases':down,'equal_length_steps':equal,'steps':rows}
