#!/usr/bin/env python3
"""Independent, standard-library-only certificate checker for DNA k=8.

This implementation uses literal strings for rotations, reverse complement,
and graph overlap. It imports neither the discovery engine nor retime.py.
Certificates are untrusted. A false check raises VerificationError.
"""
from __future__ import annotations
from collections import Counter, deque
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json

ALPHABET = 'ACGT'
K = 8
VERTICES = [''.join(t) for t in product(ALPHABET, repeat=K-1)]
EDGES = [''.join(t) for t in product(ALPHABET, repeat=K)]
VI = {v: i for i, v in enumerate(VERTICES)}
TRANSLATION = str.maketrans('ACGT', 'TGCA')
ROOT = Path(__file__).resolve().parents[1]

class VerificationError(ValueError):
    pass

def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)

def rc(word: str) -> str:
    return word.translate(TRANSLATION)[::-1]

def read_set(path: Path) -> set[str]:
    words = path.read_text(encoding='ascii').splitlines()
    require(all(len(w) == K and set(w) <= set(ALPHABET) for w in words), 'INVALID_WORD')
    require(len(set(words)) == len(words), 'DUPLICATE_WORD')
    return set(words)

def digest_set(selected: set[str]) -> str:
    return hashlib.sha256(('\n'.join(sorted(selected))+'\n').encode('ascii')).hexdigest()

def necklaces_check(selected: set[str]) -> int:
    seen: set[str] = set()
    necklaces = 0
    for w in EDGES:
        if w in seen:
            continue
        orbit = {w[i:] + w[:i] for i in range(K)}
        require(len(orbit & selected) == 1, 'NECKLACE_SELECTION')
        necklaces += 1
        seen.update(orbit)
    require(necklaces == 8230, 'NECKLACE_COUNT')
    return necklaces

def graph_check(selected: set[str]) -> dict:
    """Recompute a topological order and longest path from scratch."""
    require(selected <= set(EDGES), 'INVALID_SELECTED_EDGE')
    require(len(selected) == 8230, 'CARDINALITY')
    require({rc(w) for w in selected} == selected, 'RC_CLOSURE')
    necklaces_check(selected)
    incoming = dict.fromkeys(VERTICES, 0)
    height = dict.fromkeys(VERTICES, 0)
    parent: dict[str, str] = {}
    for w in EDGES:
        if w not in selected:
            incoming[w[1:]] += 1
    queue = deque(v for v in VERTICES if incoming[v] == 0)
    processed = 0
    while queue:
        u = queue.popleft()
        processed += 1
        for a in ALPHABET:
            w = u + a
            if w in selected:
                continue
            v = w[1:]
            if height[v] < height[u] + 1:
                height[v] = height[u] + 1
                parent[v] = w
            incoming[v] -= 1
            if incoming[v] == 0:
                queue.append(v)
    require(processed == len(VERTICES), 'RESIDUAL_DIRECTED_CYCLE')
    end = max(VERTICES, key=height.__getitem__)
    path: list[str] = []
    current = end
    while current in parent:
        w = parent[current]
        path.append(w)
        current = w[:-1]
    path.reverse()
    return {'L': height[end], 'rank': [height[v] for v in VERTICES],
            'longest_path': path, 'vertices_processed': processed,
            'residual_edges': len(EDGES)-len(selected), 'necklaces':8230,
            'sha256_words': digest_set(selected)}

def walk_check(words: list[str], selected: set[str], *, cycle: bool = False) -> tuple[str, str, int]:
    require(isinstance(words, list) and len(words) > 0, 'EMPTY_WALK')
    require(all(isinstance(w,str) and len(w) == K and set(w) <= set(ALPHABET) for w in words), 'WALK_WORD')
    require(all(a[1:] == b[:-1] for a,b in zip(words,words[1:])), 'WALK_OVERLAP')
    start, end = words[0][:-1], words[-1][1:]
    if cycle:
        require(start == end, 'CYCLE_NOT_CLOSED')
    return start, end, sum(w in selected for w in words)

def set_certificate_check(selected: set[str], cert: dict, *, recompute: bool = True) -> dict:
    require(digest_set(selected) == cert['sha256_words'], 'SET_DIGEST')
    L, ranks = cert['L'], cert['rank']
    require(type(L) is int and L >= 0, 'L_TYPE')
    require(len(ranks) == len(VERTICES), 'RANK_SIZE')
    require(all(type(x) is int and 0 <= x <= L for x in ranks), 'RANK_RANGE')
    for w in EDGES:
        if w not in selected:
            require(ranks[VI[w[1:]]] >= ranks[VI[w[:-1]]] + 1, 'RANK_EDGE')
    path = cert['longest_path']
    start, end, hits = walk_check(path, selected)
    require(hits == 0 and len(path) == L, 'LONGEST_PATH')
    require(len({w[:-1] for w in path} | {end}) == L+1, 'PATH_NOT_SIMPLE')
    result = graph_check(selected) if recompute else {'L': L}
    require(result['L'] == L, 'RECOMPUTED_LENGTH')
    if 'critical_cycle' in cert:
        cyc = cert['critical_cycle']
        _, _, hits = walk_check(cyc, selected, cycle=True)
        require(hits > 0, 'UNCUT_CYCLE')
        require(len({w[:-1] for w in cyc}) == len(cyc), 'CRITICAL_CYCLE_NOT_SIMPLE')
        require(len(cyc) == cert['critical_cycle_length'] and hits == cert['critical_cycle_hits'], 'CYCLE_COUNTS')
        # Cycle counts stay constant under complete moves. Summing the gaps
        # between hits gives L >= ceil(|Q| / hits) - 1 in the component.
        lower = (len(cyc) + hits - 1)//hits - 1
        require(lower == cert['component_lower_bound'], 'CYCLE_BOUND')
        require(lower <= L, 'INVALID_COMPONENT_BOUND')
        result['component_lower_bound'] = lower
    return {k:v for k,v in result.items() if k not in ('rank','longest_path')}

def apply_patch(selected: set[str], u: str, letters: str) -> set[str]:
    require(len(u) == K-1 and set(u) <= set(ALPHABET), 'PATCH_VERTEX')
    require(0 < len(letters) < 4 and len(set(letters)) == len(letters)
            and set(letters) <= set(ALPHABET), 'PATCH_LETTERS')
    delta: Counter[str] = Counter()
    for a in letters:
        delta[u+a] += 1
        delta[a+u] -= 1
        delta[rc(u+a)] += 1
        delta[rc(a+u)] -= 1
    new = selected.copy()
    for w, change in delta.items():
        value = int(w in selected) + change
        require(value in (0,1), 'PATCH_NONBOOLEAN')
        if value:
            new.add(w)
        else:
            new.discard(w)
    return new

def retiming_check(base: set[str], target: set[str], cert: dict) -> None:
    T, epsilon = cert['T'], cert['epsilon']
    H, f, ranks = cert['H'], cert['cut_potential'], cert['bounded_rank']
    require(type(T) is int and T >= 1 and epsilon in (0,1), 'RETIMING_PARAMETERS')
    require(all(len(a) == len(VERTICES) for a in (H,f,ranks)), 'RETIMING_ARRAY_LENGTH')
    require(all(type(x) is int for a in (H,f,ranks) for x in a), 'RETIMING_INTEGER')
    for i,v in enumerate(VERTICES):
        j = VI[rc(v)]
        require(H[i]+H[j] == T-1-T*epsilon, 'RETIMING_H_RC_SUM')
        require(f[i]+f[j] == epsilon, 'RETIMING_F_RC_SUM')
        require(f[i] == -(H[i]//T), 'RETIMING_QUOTIENT')
        require(ranks[i] == H[i]%T, 'RETIMING_REMAINDER')
    for w in EDGES:
        i,j = VI[w[:-1]], VI[w[1:]]
        require(H[j]-H[i] >= 1-T*int(w in base), 'RETIMING_DIFFERENCE')
        require(int(w in target) == int(w in base)+f[i]-f[j], 'RETIMING_GRADIENT')
        if w not in target:
            require(ranks[j] > ranks[i], 'RETIMING_RANK')

def chain_check(root: Path, chain: dict) -> dict:
    selected = read_set(root/chain['initial'])
    require(digest_set(selected) == chain['initial_sha256'], 'CHAIN_INITIAL')
    graph_check(selected)
    steps = []
    for index,step in enumerate(chain['steps'],1):
        if step['type'] == 'partial_patch':
            target = apply_patch(selected, step['u'], step['letters'])
        elif step['type'] == 'retiming':
            target = read_set(root/step['target'])
            cert = json.loads((root/step['certificate']).read_text())
            retiming_check(selected,target,cert)
        else:
            raise VerificationError('CHAIN_UNKNOWN_STEP')
        require(digest_set(target) == step['sha256_after'], 'CHAIN_DIGEST')
        result = graph_check(target)
        require(result['L'] == step['L_after'], 'CHAIN_LENGTH')
        steps.append({'step':index,'type':step['type'],'L':result['L']})
        selected=target
    final = read_set(root/chain['final'])
    require(selected == final, 'CHAIN_FINAL')
    return {'verified_steps':len(steps),'steps':steps,'final_L':steps[-1]['L']}

def manifest_check(root: Path) -> int:
    path = root/'MANIFEST_SHA256.txt'
    require(path.exists(), 'MANIFEST_MISSING')
    entries=path.read_text().splitlines()
    seen=set()
    for line in entries:
        expected,name=line.split('  ',1)
        p = root/name
        require(p.resolve().is_relative_to(root.resolve()), 'MANIFEST_PATH')
        require(p.is_file() and name not in seen, 'MANIFEST_FILE')
        require(hashlib.sha256(p.read_bytes()).hexdigest() == expected, 'MANIFEST_HASH:'+name)
        seen.add(name)
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()
            and p.name != 'MANIFEST_SHA256.txt' and '__pycache__' not in p.parts}
    require(seen == actual, 'MANIFEST_COVERAGE')
    return len(seen)
