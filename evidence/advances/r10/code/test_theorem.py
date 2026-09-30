'Exhaustive tiny-graph regression, NOT the proof of the general theorem.'
from itertools import product
from collections import deque
from pathlib import Path
import argparse, json

def length(n, edges):
    deg = [0] * n
    ad = [[] for _ in range(n)]
    d = [0] * n
    for u, v in edges:
        ad[u].append(v)
        deg[v] += 1
    q = deque((i for i in range(n) if deg[i] == 0))
    num = 0
    while q:
        u = q.popleft()
        num += 1
        for v in ad[u]:
            d[v] = max(d[v], d[u] + 1)
            deg[v] -= 1
            if deg[v] == 0:
                q.append(v)
    return (max(d, default=0), d) if num == n else (None, d)

def run():
    results = []
    rounding = 0
    for n in [3, 4]:
        th = lambda i: n - 1 - i
        pairs = [(u, v) for u in range(n) for v in range(n)]
        orbits = []
        seen = set()
        for e in pairs:
            if e in seen:
                continue
            u, v = e
            o = {e, (th(v), th(u))}
            seen |= o
            orbits.append(o)
        size = 1 << len(orbits)
        lens = []
        weights = []
        for mask in range(size):
            es = set().union(*(orbits[i] for i in range(len(orbits)) if mask >> i & 1))
            L, d = length(n, es)
            lens.append(L)
            weights.append(len(es))
            if L is not None:
                hh = []
                for u in range(n):
                    z = d[u] - d[th(u)]
                    hh.append((z + 1) // 2 if z >= 0 else -((-z + 1) // 2))
                if not (all((hh[th(u)] == -hh[u] for u in range(n))) and all((hh[v] > hh[u] for u, v in es))):
                    raise AssertionError('Verification invariant at original line 36')
                if not max(map(abs, hh), default=0) <= (L + 1) // 2:
                    raise AssertionError('Verification invariant at original line 37')
                rounding += 1
        for H in [0, 1, 2]:
            opt = [weights[m] if lens[m] is not None and lens[m] <= 2 * H else -1 for m in range(size)]
            for i in range(len(orbits)):
                for m in range(size):
                    if m >> i & 1:
                        opt[m] = max(opt[m], opt[m ^ 1 << i])
            reps = [u for u in range(n) if u < th(u)]
            heightmasks = []
            for labs in product(range(-H, H + 1), repeat=len(reps)):
                h = [0] * n
                for u, x in zip(reps, labs):
                    h[u] = x
                    h[th(u)] = -x
                mask = sum((1 << i for i, o in enumerate(orbits) if all((h[v] > h[u] for u, v in o))))
                heightmasks.append(mask)
            for m in range(size):
                if not max((weights[m & hm] for hm in heightmasks)) == opt[m]:
                    raise AssertionError('Verification invariant at original line 48')
            results.append({'vertices': n, 'H': H, 'RC_ambient_graphs': size, 'height_assignments': len(heightmasks)})
    return {'status': 'PASS', 'scope': 'all RC-invariant directed ambient graphs on n=3,4 including loops, theta(i)=n-1-i; H=0,1,2', 'cases': results, 'total_ambient_and_bound_pairs': sum((r['RC_ambient_graphs'] for r in results)), 'acyclic_rounding_checks': rounding}
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    r = run()
    (a.out / 'theorem_regression.json').write_text(json.dumps(r, indent=2))
    print(json.dumps(r, indent=2))
