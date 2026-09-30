"""Independent literal-DNA verification primitives (standard library only).
No imports from discovery, upstream, NumPy, SciPy or a solver.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
ABC = 'ACGT'
TRANS = str.maketrans('ACGT', 'TGCA')
VERTICES = tuple(map(''.join, product(ABC, repeat=7)))
WORDS = tuple(v + a for v in VERTICES for a in ABC)
VI = {v: i for i, v in enumerate(VERTICES)}
WI = {w: i for i, w in enumerate(WORDS)}
WORD_SET = frozenset(WORDS)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def rc(word):
    return word.translate(TRANS)[::-1]


VR = tuple(v for v in VERTICES if v < rc(v))
ER = tuple(w for w in WORDS if w <= rc(w))


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def word_digest(selected):
    return hashlib.sha256(('\n'.join(sorted(selected))+'\n').encode('ascii')).hexdigest()


def read_set(path, card=None, symmetric=True):
    words = Path(path).read_text(encoding='ascii').splitlines()
    need(all(w in WORD_SET for w in words), 'WORD_ALPHABET_OR_LENGTH')
    need(len(words) == len(set(words)), 'DUPLICATE_WORD')
    s = frozenset(words)
    if card is not None:
        need(len(s) == card, 'CARDINALITY')
    if symmetric:
        need({rc(w) for w in s} == s, 'RC_SET')
    return s


@lru_cache(None)
def pcrs():
    seen, out = set(), []
    for w in WORDS:
        if w in seen:
            continue
        p = tuple(sorted({w[i:] + w[:i] for i in range(8)}))
        seen.update(p)
        out.append(p)
    need(len(out) == 8230 and seen == WORD_SET, 'PCR_CENSUS')
    return tuple(out)


@lru_cache(None)
def groups():
    seen, out = set(), []
    for p in pcrs():
        if p[0] in seen:
            continue
        q = tuple(sorted(rc(w) for w in p))
        full = frozenset(p) | frozenset(q)
        seen.update(full)
        out.append({'words': full, 'stable': p == q, 'representative': min(full)})
    need(len(out) == 4179 and sum(g['stable'] for g in out) == 128, 'GROUP_CENSUS')
    return tuple(out)


def graph(selected):
    """Rebuild the entire residual lower graph, including isolated vertices."""
    adj = {v: [] for v in VERTICES}
    indeg = dict.fromkeys(VERTICES, 0)
    for w in WORDS:
        if w not in selected:
            adj[w[:-1]].append(w)
            indeg[w[1:]] += 1
    queue = deque(v for v in VERTICES if indeg[v] == 0)
    order, parent = [], {}
    depth = dict.fromkeys(VERTICES, 0)
    while queue:
        u = queue.popleft()
        order.append(u)
        for w in adj[u]:
            v = w[1:]
            if depth[v] < depth[u] + 1:
                depth[v] = depth[u] + 1
                parent[v] = w
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    need(len(order) == len(VERTICES), 'RESIDUAL_DIRECTED_CYCLE')
    suffix = dict.fromkeys(VERTICES, 0)
    for u in reversed(order):
        suffix[u] = max([0] + [suffix[w[1:]]+1 for w in adj[u]])
    end = max(VERTICES, key=depth.get)
    path = []
    while end in parent:
        w = parent[end]
        path.append(w)
        end = w[:-1]
    path.reverse()
    return depth, suffix, path


def walk(words, closed=False, simple=False):
    need(isinstance(words, list) and len(words) > 0, 'EMPTY_WALK')
    need(all(w in WORD_SET for w in words), 'WALK_WORDS')
    need(all(a[1:] == b[:-1] for a, b in zip(words, words[1:])), 'WALK_CONTINUITY')
    vertices = [w[:-1] for w in words]
    if closed:
        need(words[-1][1:] == words[0][:-1], 'WALK_CLOSURE')
    else:
        vertices.append(words[-1][1:])
    if simple:
        need(len(vertices) == len(set(vertices)), 'WALK_NOT_SIMPLE')


def height(c, expected=None, full_graph=True):
    need(c.get('half_height') == 53 and type(c.get('half_height')) is int,
         'HEIGHT_HALF_BOUND')
    need(c.get('target_L') == 106 and type(c.get('target_L')) is int,
         'HEIGHT_TARGET')
    h = c.get('heights')
    need(isinstance(h, list) and len(h) == len(VERTICES), 'HEIGHT_SIZE')
    need(all(type(x) is int for x in h), 'HEIGHT_INTEGRAL')
    need(all(-53 <= x <= 53 for x in h), 'HEIGHT_DOMAIN')
    r = dict(zip(VERTICES, h))
    need(all(r[rc(v)] == -r[v] for v in VERTICES), 'HEIGHT_RC')
    s = frozenset(w for w in WORDS if r[w[1:]] <= r[w[:-1]])
    need(type(c.get('cardinality')) is int and len(s) == c['cardinality'], 'HEIGHT_COST')
    if expected is not None:
        need(s == expected, 'HEIGHT_SET_JOIN')
    need({rc(w) for w in s} == s, 'HEIGHT_SELECTED_RC')
    counts = Counter(sum(w in s for w in p) for p in pcrs())
    need(0 not in counts, 'PCR_UNHIT')
    need(len(s)-8230 == sum((k-1)*v for k, v in counts.items()), 'SURPLUS_IDENTITY')
    result = {'cardinality': len(s), 'surplus': len(s)-8230,
              'PCR_cut_histogram': {str(k): v for k, v in sorted(counts.items())}}
    if full_graph:
        _, _, path = graph(s)
        need(len(path) <= 106, 'HEIGHT_LENGTH')
        result['L'] = len(path)
    return s, result


def rank_path(c, selected):
    need(c.get('cardinality') == len(selected), 'RANK_CARDINALITY')
    need(c.get('selected_sha256') == word_digest(selected), 'RANK_SET_HASH')
    L, d = c.get('L'), c.get('rank')
    need(type(L) is int and L >= 0, 'RANK_LENGTH_TYPE')
    need(isinstance(d, list) and len(d) == len(VERTICES), 'RANK_SIZE')
    need(all(type(x) is int and 0 <= x <= L for x in d), 'RANK_DOMAIN')
    ranks = dict(zip(VERTICES, d))
    for w in WORDS:
        if w not in selected:
            need(ranks[w[1:]] >= ranks[w[:-1]] + 1, 'RANK_EDGE')
    p = c.get('path')
    walk(p, simple=True)
    need(len(p) == L, 'PATH_LENGTH')
    need(not (set(p) & selected), 'PATH_SELECTED_EDGE')
    string = p[0][:-1]+''.join(w[-1] for w in p)
    need(c.get('avoiding_word') == string, 'AVOIDING_WORD_JOIN')
    need(len(string) == L+7 and all(string[i:i+8] not in selected for i in range(L)),
         'AVOIDING_WORD_WINDOWS')
    # A second graph computation checks that the supplied rank is a longest-prefix rank.
    rebuilt, _, found = graph(selected)
    need([rebuilt[v] for v in VERTICES] == d and len(found) == L, 'LONGEST_RANK_REBUILD')
    return {'cardinality': len(selected), 'L': L, 'residual_arcs': 65536-len(selected),
            'avoiding_length': len(string), 'all_words_hit_at_length': len(string)+1}


def irredundancy(c, selected):
    required = []
    for p in pcrs():
        chosen = sorted(set(p) & selected)
        need(chosen, 'IRREDUNDANCY_PCR_UNHIT')
        if len(chosen) > 1:
            required.append({'representative': p[0], 'words': list(p), 'selected': chosen})
    need(c.get('overselected_PCRs') == required, 'IRREDUNDANCY_COMPLETE_PCRS')
    exception = {w for x in required for w in x['selected']}
    cycles = c.get('cycles')
    need(isinstance(cycles, list), 'IRREDUNDANCY_CYCLE_LIST')
    labels = [x.get('selected_edge') for x in cycles]
    need(len(labels) == len(set(labels)) and set(labels) == exception, 'IRREDUNDANCY_COVERAGE')
    for x in cycles:
        p = x['cycle']
        walk(p, closed=True, simple=True)
        need(set(p) & selected == {x['selected_edge']}, 'IRREDUNDANCY_SINGLE_CUT')
    return {'ordinary_inclusion_minimal': True,
            'singleton_PCR_cycles': len(pcrs())-len(required),
            'extra_cycles': len(cycles),
            'extra_cycle_lengths': [len(x['cycle']) for x in cycles]}


def block_summary(c, before, after):
    gs = groups()
    ids = c.get('free_groups')
    need(isinstance(ids, list) and all(type(i) is int and 0 <= i < len(gs) for i in ids),
         'BLOCK_GROUP_IDS')
    need(ids == sorted(set(ids)), 'BLOCK_GROUP_ORDER_OR_DUPLICATE')
    stable = {i for i, g in enumerate(gs) if g['stable']}
    need(stable <= set(ids), 'BLOCK_ALL_STABLE_FREE')
    free = set().union(*(gs[i]['words'] for i in ids))
    need(before ^ after <= free, 'BLOCK_FROZEN_INDICATORS')
    need(c.get('before') == len(before) and c.get('after') == len(after), 'BLOCK_COST')
    if c.get('incumbent_available'):
        need(c.get('selected_cardinality') == len(after), 'BLOCK_POSTPROCESSING_JOIN')
        need(c.get('selected_difference') == len(before ^ after), 'BLOCK_DELTA')
    ns, st = len(set(ids)-stable), len(stable)
    nx = sum(w in free for w in ER)
    rows = sum(w in free or w not in before for w in ER)
    rows += sum(p[0] in free for p in pcrs()) + 1
    need((c.get('selection_variables'), c.get('variables'), c.get('rows')) ==
         (nx, nx+len(VR), rows), 'BLOCK_MODEL_COUNTS')
    need((c.get('stable_groups_free'), c.get('nonstable_groups_free')) == (st, ns),
         'BLOCK_GROUP_COUNTS')
    # Check the rank domains and the simultaneous existence of free real ranks,
    # using doubled exact ranks rather than the floating-point solver values.
    fixed_deleted = before | free
    d0, b0, _ = graph(fixed_deleted)
    d, _, path = graph(after)
    need(c.get('L') == len(path) and len(path) <= 106, 'BLOCK_LENGTH')
    r2 = {v: 106+d[v]-d[rc(v)] for v in VERTICES}
    for v in VERTICES:
        need(2*d0[v] <= r2[v] <= 2*(106-b0[v]), 'BLOCK_REAL_RANK_DOMAIN')
        need(r2[v] + r2[rc(v)] == 212, 'BLOCK_REAL_RANK_RC')
    for w in WORDS:
        if w not in after:
            need(r2[w[1:]]-r2[w[:-1]] >= 2, 'BLOCK_REAL_RANK_EDGE')
    return {'nonstable_groups_free': ns, 'stable_groups_free': st,
            'binary_variables': nx, 'real_ranks': len(VR), 'rows': rows,
            'before': len(before), 'after': len(after), 'changed_words': len(before ^ after)}


def surplus_profile(selected):
    """Exact paired-PCR bookkeeping, with finite stable-domain census."""
    need({rc(w) for w in selected} == set(selected), 'SURPLUS_RC')
    stable_excess = paired_excess = 0
    stable_fixed_total = 0
    for g in groups():
        w = g['representative']
        p = {w[i:]+w[:i] for i in range(8)}
        count = len(p & selected)
        need(count >= 1, 'SURPLUS_PCR_UNHIT')
        if g['stable']:
            need(sum(rc(x) == x for x in p) == 2, 'STABLE_FIXED_CENSUS')
            stable_fixed_total += 2
            stable_excess += count-1
        else:
            q = {rc(x) for x in p}
            need(len(q & selected) == count, 'SURPLUS_PAIR_EQUALITY')
            paired_excess += count-1
    need(len(selected)-8230 == 2*paired_excess+stable_excess, 'PAIRED_SURPLUS_IDENTITY')
    return {'stable_PCRs': 128, 'nonstable_PCR_pairs': 4051,
            'individually_RC_fixed_words': stable_fixed_total,
            'selected_RC_fixed_words': sum(rc(w) == w for w in selected),
            'paired_excess': paired_excess, 'stable_excess': stable_excess,
            'total_excess': 2*paired_excess+stable_excess}
