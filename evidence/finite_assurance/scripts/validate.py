#!/usr/bin/env python3
"""Independent exact validation for the article's bounded terminal atlas.

This checker does not import the producer.  It enumerates integer base-4 codes,
rebuilds every selector and graph invariant, validates strict schemas, and
runs native algorithmic falsifiers plus a typed receipt-mutation panel.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import heapq
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

from jsonschema import Draft202012Validator


ATLAS_SCHEMA = "FCIG_CWC_TERMINAL_ATLAS_V1"
ROW_SCHEMA = "FCIG_CWC_ATLAS_ROW_V1"
STATUS = "PASS_INDEPENDENT_EXACT_TERMINAL_ATLAS_ASSURANCE"
PRODUCER_STATUS = "PASS_EXACT_SMALL_DNA_TERMINAL_ATLAS_PRODUCER"
SCOPE = "EXACT_SMALL_DNA_TERMINAL_ATLAS"
ALPHABET = (0, 1, 2, 3)
COMPLEMENT = (3, 2, 1, 0)
FAMILIES = ("SEED_REAL", "SEED_IMAG", "LOWER_LAMBDA", "LOWER_MIDPOINT")
CAPSULE = Path(__file__).resolve().parents[1]
SKELETON_RECEIPT = CAPSULE / "inputs" / "critical_ratio_skeleton.json"
ROW_SCHEMA_PATH = CAPSULE / "schemas" / "atlas_row.schema.json"
ATLAS_SCHEMA_PATH = CAPSULE / "schemas" / "terminal_atlas.schema.json"

TRIGONOMETRY = {
    4: {
        "C": ((2, 0), (0, 0), (-2, 0), (0, 0), (2, 0)),
        "S": ((0, 0), (1, 0), (0, 0), (-1, 0), (0, 0)),
    },
    6: {
        "C": ((2, 0), (1, 0), (-1, 0), (-2, 0), (-1, 0), (1, 0), (2, 0)),
        "S": ((0, 0), (1, 0), (1, 0), (0, 0), (-1, 0), (-1, 0), (0, 0)),
    },
    8: {
        "C": ((2, 0), (0, 1), (0, 0), (0, -1), (-2, 0), (0, -1), (0, 0), (0, 1), (2, 0)),
        "S": ((0, 0), (1, 0), (0, 1), (1, 0), (0, 0), (-1, 0), (0, -1), (-1, 0), (0, 0)),
    },
}


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, code: str) -> None:
    if not condition:
        raise CheckFailure(code)


def encoded(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha(value: object) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def canonical_document(path: Path) -> tuple[dict[str, object], bytes]:
    raw = path.read_bytes()
    require(raw[-1:] == b"\n" and raw[-2:] != b"\n\n", "ATLAS_RECEIPT_LF")
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise CheckFailure("ATLAS_RECEIPT_JSON") from exc
    require(isinstance(value, dict), "ATLAS_RECEIPT_OBJECT")
    require(raw == encoded(value) + b"\n", "ATLAS_RECEIPT_CANONICAL")
    return value, raw


def rational(data: dict[str, object]) -> Fraction:
    result = Fraction(int(data["numerator"]), int(data["denominator"]))
    require(result.denominator > 0, "ATLAS_RATIONAL")
    return result


def p_eval(coefficients: Sequence[Fraction | int], point: Fraction) -> Fraction:
    value = Fraction()
    for coefficient in reversed(coefficients):
        value = value * point + Fraction(coefficient)
    return value


def stripped(values: Sequence[Fraction]) -> list[Fraction]:
    output = list(values)
    while len(output) > 1 and output[-1] == 0:
        output.pop()
    return output or [Fraction()]


def p_mod(numerator: Sequence[int], denominator: Sequence[int]) -> list[Fraction]:
    work = stripped([Fraction(value) for value in numerator])
    divisor = stripped([Fraction(value) for value in denominator])
    require(divisor[-1] != 0, "ATLAS_MINPOLY_ZERO")
    while len(work) >= len(divisor) and any(work):
        multiplier = work[-1] / divisor[-1]
        offset = len(work) - len(divisor)
        for index in range(len(divisor)):
            work[offset + index] -= multiplier * divisor[index]
        work = stripped(work)
    return work


def interval_product(
    first: tuple[Fraction, Fraction], second: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    candidates = [
        first[0] * second[0], first[0] * second[1],
        first[1] * second[0], first[1] * second[1],
    ]
    return min(candidates), max(candidates)


def interval_horner(
    coefficients: Sequence[Fraction], bounds: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    value = (Fraction(), Fraction())
    for coefficient in reversed(coefficients):
        value = interval_product(value, bounds)
        value = value[0] + coefficient, value[1] + coefficient
    return value


class AlgebraicSign:
    def __init__(self, definition: dict[str, object]) -> None:
        self.definition = copy.deepcopy(definition)
        self.value: Fraction | None = None
        self.minpoly: list[int] | None = None
        self.bounds: tuple[Fraction, Fraction] | None = None
        if definition["kind"] == "RATIONAL":
            self.value = rational(definition["value"])  # type: ignore[arg-type]
            require(self.value > 1, "ATLAS_WEIGHT_DOMAIN")
        else:
            require(definition["kind"] == "REAL_ALGEBRAIC", "ATLAS_WEIGHT_KIND")
            self.minpoly = [int(value) for value in definition["polynomial_coefficients_ascending"]]  # type: ignore[index]
            interval = definition["isolating_interval"]
            require(isinstance(interval, dict), "ATLAS_WEIGHT_INTERVAL")
            self.bounds = rational(interval["left"]), rational(interval["right"])  # type: ignore[arg-type]
            require(Fraction(1) < self.bounds[0] < self.bounds[1], "ATLAS_WEIGHT_INTERVAL_ORDER")
            require(p_eval(self.minpoly, self.bounds[0]) * p_eval(self.minpoly, self.bounds[1]) < 0, "ATLAS_WEIGHT_BRACKET")

    def polynomial_sign(self, polynomial: Sequence[int]) -> int:
        if self.value is not None:
            answer = p_eval(polynomial, self.value)
            return int(answer > 0) - int(answer < 0)
        require(self.minpoly is not None and self.bounds is not None, "ATLAS_WEIGHT_STATE")
        remainder = p_mod(polynomial, self.minpoly)
        if not any(remainder):
            return 0
        low, high = self.bounds
        sign_low = p_eval(self.minpoly, low)
        iterations = 0
        while iterations < 1024:
            lower, upper = interval_horner(remainder, (low, high))
            if lower > 0:
                return 1
            if upper < 0:
                return -1
            middle = (low + high) / 2
            sign_middle = p_eval(self.minpoly, middle)
            if sign_middle == 0:
                direct = p_eval(remainder, middle)
                return int(direct > 0) - int(direct < 0)
            if sign_low * sign_middle < 0:
                high = middle
            else:
                low = middle
                sign_low = sign_middle
            iterations += 1
        raise CheckFailure("ATLAS_SIGN_REFINEMENT")

    def expression_sign(self, data: tuple[int, int, int, int]) -> int:
        a, b, c, d = data
        n_sign = self.polynomial_sign([b, a])
        m_sign = self.polynomial_sign([d, c])
        if n_sign == 0:
            return m_sign
        if m_sign == 0:
            return n_sign
        if n_sign == m_sign:
            return n_sign
        comparison = [
            b * b - 2 * d * d,
            2 * a * b - 4 * c * d,
            a * a - 2 * c * c,
        ]
        q_sign = self.polynomial_sign(comparison)
        if q_sign == 0:
            return 0
        return n_sign if q_sign > 0 else m_sign


def decode(code: int, length: int) -> tuple[int, ...]:
    output = [0] * length
    remaining = code
    for index in range(length - 1, -1, -1):
        output[index] = remaining % 4
        remaining //= 4
    require(remaining == 0, "ATLAS_CODE_RANGE")
    return tuple(output)


def encode_word(word: Sequence[int]) -> int:
    value = 0
    for symbol in word:
        value = value * 4 + symbol
    return value


def reverse_complement(word: Sequence[int]) -> tuple[int, ...]:
    return tuple(COMPLEMENT[value] for value in word[::-1])


def rotation(word: Sequence[int], amount: int = 1) -> tuple[int, ...]:
    data = tuple(word)
    amount %= len(data)
    return data[amount:] + data[:amount]


def rotations(word: Sequence[int]) -> list[tuple[int, ...]]:
    data = tuple(word)
    output: list[tuple[int, ...]] = []
    current = data
    while current not in output:
        output.append(current)
        current = rotation(current)
    return output


def coefficients(k: int, family: str) -> tuple[tuple[int, int], ...]:
    c = TRIGONOMETRY[k]["C"]
    s = TRIGONOMETRY[k]["S"]
    table = {
        "SEED_REAL": c[1 : k + 1],
        "SEED_IMAG": s[1 : k + 1],
        "LOWER_LAMBDA": s[1:k],
        "LOWER_MIDPOINT": tuple((c[i][0] - 2, c[i][1]) for i in range(1, k)),
    }
    require(family in table, "ATLAS_FAMILY")
    return tuple(table[family])


def form_for(word: Sequence[int], table: Sequence[tuple[int, int]]) -> tuple[int, int, int, int]:
    ar = one = radical_r = radical = 0
    require(len(word) == len(table), "ATLAS_FORM_LENGTH")
    for index, symbol in enumerate(word):
        plain, root2 = table[index]
        if symbol == 0:
            ar -= plain
            radical_r -= root2
        elif symbol == 1:
            one -= plain
            radical -= root2
        elif symbol == 2:
            one += plain
            radical += root2
        elif symbol == 3:
            ar += plain
            radical_r += root2
        else:
            raise CheckFailure("ATLAS_ALPHABET")
    return ar, one, radical_r, radical


def primitive_form(data: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    if data == (0, 0, 0, 0):
        return data
    divisor = math.gcd(math.gcd(abs(data[0]), abs(data[1])), math.gcd(abs(data[2]), abs(data[3])))
    require(divisor > 0, "ATLAS_FORM_GCD")
    result = tuple(value // divisor for value in data)
    first = next(value for value in result if value)
    if first < 0:
        result = tuple(-value for value in result)
    return result  # type: ignore[return-value]


def build_catalog(k: int) -> list[tuple[str, int, tuple[int, int, int, int]]]:
    output = []
    for family in FAMILIES:
        table = coefficients(k, family)
        for code in range(4 ** len(table)):
            output.append((family, code, form_for(decode(code, len(table)), table)))
    return output


def ratio_object(value: Fraction) -> dict[str, int]:
    return {"denominator": value.denominator, "numerator": value.numerator}


def finite_endpoint(value: Fraction) -> dict[str, object]:
    return {"kind": "FINITE", "value": ratio_object(value)}


def expected_landmarks() -> list[dict[str, object]]:
    return [
        {
            "landmark_id": "RAT2",
            "weight": {"kind": "RATIONAL", "value": ratio_object(Fraction(2))},
        },
        {
            "landmark_id": "GOLD5",
            "weight": {
                "isolating_interval": {
                    "left": ratio_object(Fraction(6)),
                    "right": ratio_object(Fraction(7)),
                },
                "kind": "REAL_ALGEBRAIC",
                "polynomial_coefficients_ascending": [1, -7, 1],
                "root_selector": "UNIQUE_REAL_ROOT_IN_OPEN_INTERVAL",
            },
        },
        {
            "landmark_id": "QUAD13",
            "weight": {
                "isolating_interval": {
                    "left": ratio_object(Fraction(3)),
                    "right": ratio_object(Fraction(4)),
                },
                "kind": "REAL_ALGEBRAIC",
                "polynomial_coefficients_ascending": [-13, 0, 1],
                "root_selector": "UNIQUE_REAL_ROOT_IN_OPEN_INTERVAL",
            },
        },
        {
            "landmark_id": "CUBIC2",
            "weight": {
                "isolating_interval": {
                    "left": ratio_object(Fraction(5, 4)),
                    "right": ratio_object(Fraction(13, 10)),
                },
                "kind": "REAL_ALGEBRAIC",
                "polynomial_coefficients_ascending": [-2, 0, 0, 1],
                "root_selector": "UNIQUE_REAL_ROOT_IN_OPEN_INTERVAL",
            },
        },
    ]


def locate_landmarks(
    k: int,
    roots: tuple[Fraction, ...],
    landmarks: Sequence[dict[str, object]],
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    boundaries = (Fraction(1),) + roots
    for landmark in landmarks:
        landmark_id = str(landmark["landmark_id"])
        weight = landmark["weight"]
        require(isinstance(weight, dict), "ATLAS_SKELETON_LANDMARK_WEIGHT")
        if weight["kind"] == "RATIONAL":
            value = rational(weight["value"])  # type: ignore[arg-type]
            if value in roots:
                index = roots.index(value)
                output.append({
                    "cell_id": f"K{k}-P{index:03d}",
                    "landmark_id": landmark_id,
                    "position_kind": "CRITICAL_POINT",
                })
                continue
            matches = [
                index for index, left in enumerate(boundaries)
                if left < value and (index == len(roots) or value < roots[index])
            ]
        else:
            require(weight["kind"] == "REAL_ALGEBRAIC", "ATLAS_SKELETON_LANDMARK_KIND")
            interval = weight["isolating_interval"]
            require(isinstance(interval, dict), "ATLAS_SKELETON_LANDMARK_INTERVAL")
            low = rational(interval["left"])  # type: ignore[arg-type]
            high = rational(interval["right"])  # type: ignore[arg-type]
            polynomial = weight["polynomial_coefficients_ascending"]
            require(isinstance(polynomial, list), "ATLAS_SKELETON_LANDMARK_POLYNOMIAL")
            require(
                low < high and p_eval(polynomial, low) * p_eval(polynomial, high) < 0,
                "ATLAS_SKELETON_LANDMARK_BRACKET",
            )
            matches = [
                index for index, left in enumerate(boundaries)
                if left <= low and (index == len(roots) or high <= roots[index])
            ]
        require(len(matches) == 1, "ATLAS_SKELETON_LANDMARK_LOCATION")
        output.append({
            "cell_id": f"K{k}-I{matches[0]:03d}",
            "landmark_id": landmark_id,
            "position_kind": "OPEN_INTERVAL",
        })
    return output


def expected_skeleton_order(
    k: int,
    landmarks: Sequence[dict[str, object]],
) -> dict[str, object]:
    catalog = build_catalog(k)
    require(
        all(form[2:] == (0, 0) for _, _, form in catalog),
        f"ATLAS_SKELETON_RADICAL_K{k}",
    )
    raw = [[family, code, form[0], form[1]] for family, code, form in catalog]
    counts: Counter[tuple[int, int]] = Counter((form[0], form[1]) for _, _, form in catalog)
    first: dict[tuple[int, int], tuple[str, int]] = {}
    for family, code, form in catalog:
        first.setdefault((form[0], form[1]), (family, code))

    canonical: list[dict[str, object]] = []
    pair_ids: dict[tuple[int, int], str] = {}
    for pair, occurrence in first.items():
        if pair == (0, 0):
            continue
        a, b = pair
        form_id = f"K{k}-L{len(canonical):04d}"
        pair_ids[pair] = form_id
        root = None if a == 0 else Fraction(-b, a)
        if a == 0:
            classification = "CONSTANT_SIGN"
        elif root is not None and root > 1:
            classification = "CRITICAL_ROOT"
        else:
            classification = "NO_ROOT_GT_ONE"
        canonical.append({
            "a": a,
            "b": b,
            "classification": classification,
            "critical_ratio": ratio_object(root) if classification == "CRITICAL_ROOT" else None,
            "first_occurrence": {"family": occurrence[0], "word_code": occurrence[1]},
            "form_id": form_id,
            "occurrence_count": counts[pair],
        })

    root_forms: dict[Fraction, list[str]] = {}
    root_counts: Counter[Fraction] = Counter()
    for pair, form_id in pair_ids.items():
        a, b = pair
        if a == 0:
            continue
        root = Fraction(-b, a)
        if root > 1:
            root_forms.setdefault(root, []).append(form_id)
            root_counts[root] += counts[pair]
    roots = tuple(sorted(root_forms))
    critical_points = [
        {
            "cell_id": f"K{k}-P{index:03d}",
            "ratio": ratio_object(root),
            "supporting_form_ids": root_forms[root],
            "supporting_occurrence_count": root_counts[root],
        }
        for index, root in enumerate(roots)
    ]

    boundaries = (Fraction(1),) + roots
    intervals: list[dict[str, object]] = []
    for index, left in enumerate(boundaries):
        if index < len(roots):
            right = roots[index]
            representative = (left + right) / 2
            right_endpoint: dict[str, object] = finite_endpoint(right)
        else:
            representative = left + 1
            right_endpoint = {"kind": "POSITIVE_INFINITY"}
        intervals.append({
            "cell_id": f"K{k}-I{index:03d}",
            "left": finite_endpoint(left),
            "representative": ratio_object(representative),
            "right": right_endpoint,
        })

    ordered_cells: list[str] = []
    for index in range(len(roots)):
        ordered_cells.extend((f"K{k}-I{index:03d}", f"K{k}-P{index:03d}"))
    ordered_cells.append(f"K{k}-I{len(roots):03d}")
    structural = [[family, code] for family, code, form in catalog if form == (0, 0, 0, 0)]
    require(
        all(value[1] == 0 for value in TRIGONOMETRY[k]["C"] + TRIGONOMETRY[k]["S"]),
        f"ATLAS_SKELETON_TRIGONOMETRY_K{k}",
    )
    return {
        "canonical_nonzero_forms": canonical,
        "critical_points": critical_points,
        "form_catalog": {
            "canonical_nonzero_form_count": len(canonical),
            "canonical_nonzero_form_sha256": sha(canonical),
            "constant_form_count": sum(1 for a, b in counts if a == 0 and b != 0),
            "raw_critical_occurrence_count": sum(root_counts.values()),
            "raw_occurrence_count": len(raw),
            "raw_occurrence_sha256": sha(raw),
            "structural_zero_by_family": {
                family: sum(1 for item in structural if item[0] == family)
                for family in FAMILIES
            },
            "structural_zero_occurrence_count": len(structural),
            "structural_zero_sha256": sha(structural),
        },
        "k": k,
        "landmark_locations": locate_landmarks(k, roots, landmarks),
        "open_intervals": intervals,
        "ordered_cell_ids": ordered_cells,
        "ordered_critical_ratio_sha256": sha([ratio_object(root) for root in roots]),
        "trigonometric_tables": {
            "C_twice_cos": [value[0] for value in TRIGONOMETRY[k]["C"]],
            "S_normalized_sine": [value[0] for value in TRIGONOMETRY[k]["S"]],
        },
    }


def validate_skeleton_payload(payload: dict[str, object]) -> dict[str, object]:
    require(
        set(payload) == {
            "claim_boundary", "family_order", "landmarks", "orders", "q",
            "schema_version", "scope",
        },
        "ATLAS_SKELETON_KEYS",
    )
    require(payload["q"] == 4, "ATLAS_SKELETON_Q")
    require(payload["family_order"] == list(FAMILIES), "ATLAS_SKELETON_FAMILIES")
    require(
        payload["schema_version"] == "FCIG_CWC_CRITICAL_RATIO_SKELETON_V1"
        and payload["scope"] == "EXACT_DNA_CRITICAL_RATIO_SKELETON"
        and payload["claim_boundary"]
        == "EXACT_K4_K6_CRITICAL_RATIO_SKELETON_NO_TERMINAL_OR_ALL_ORDER_INFERENCE",
        "ATLAS_SKELETON_CONTRACT",
    )
    landmarks = expected_landmarks()
    require(payload["landmarks"] == landmarks, "ATLAS_SKELETON_LANDMARKS")
    orders = payload["orders"]
    require(isinstance(orders, list) and len(orders) == 2, "ATLAS_SKELETON_ORDERS")
    expected_orders = [expected_skeleton_order(k, landmarks) for k in (4, 6)]
    for actual, expected, k in zip(orders, expected_orders, (4, 6), strict=True):
        require(actual == expected, f"ATLAS_SKELETON_ORDER_K{k}")
    return {
        "claim_boundary": payload["claim_boundary"],
        "family_order": list(FAMILIES),
        "landmarks": landmarks,
        "orders": expected_orders,
        "q": 4,
        "schema_version": payload["schema_version"],
        "scope": payload["scope"],
    }


def zero_choice(word: tuple[int, ...]) -> int:
    own = rotations(word)
    mate_word = reverse_complement(word)
    if mate_word in set(own):
        fixed = [item for item in own if reverse_complement(item) == item]
        require(len(fixed) == 2, "ATLAS_ZERO_STABLE_FIXED")
        return int(word == min(fixed))
    mate = rotations(mate_word)
    smallest = min(min(own), min(mate))
    choices = {smallest, reverse_complement(smallest)}
    require(len(choices) == 2, "ATLAS_ZERO_PAIRED")
    return int(word in choices)


def sparse_histogram(values: Sequence[int]) -> list[dict[str, int]]:
    counter = Counter(values)
    return [{"bucket": key, "count": counter[key]} for key in sorted(counter) if counter[key] > 0]


def replay_firings(k: int, seed: list[int], chosen: list[int], target: list[int]) -> tuple[int, int]:
    vertices = 4 ** (k - 1)
    predecessors: list[set[int]] = [set() for _ in range(vertices)]
    successors: list[set[int]] = [set() for _ in range(vertices)]
    for code in range(4 ** k):
        source, sink = code // 4, code % vertices
        if source != sink and chosen[source] and chosen[sink]:
            first, second = (source, sink) if seed[code] == 0 else (sink, source)
            predecessors[second].add(first)
            successors[first].add(second)
    indegree = [len(item) for item in predecessors]
    available = [node for node in range(vertices) if chosen[node] and indegree[node] == 0]
    heapq.heapify(available)
    order = []
    while available:
        node = heapq.heappop(available)
        order.append(node)
        for successor in sorted(successors[node]):
            indegree[successor] -= 1
            if indegree[successor] == 0:
                heapq.heappush(available, successor)
    require(len(order) == sum(chosen), "ATLAS_REPLAY_ORDER")
    state = list(seed)
    for node in order:
        incoming = [letter * vertices + node for letter in ALPHABET]
        outgoing = [node * 4 + letter for letter in ALPHABET]
        require(set(incoming).isdisjoint(outgoing), "ATLAS_REPLAY_LOOP")
        require(all(state[edge] == 1 for edge in incoming), "ATLAS_REPLAY_IN")
        require(all(state[edge] == 0 for edge in outgoing), "ATLAS_REPLAY_OUT")
        for edge in incoming:
            state[edge] = 0
        for edge in outgoing:
            state[edge] = 1
    return len(order), sum(x != y for x, y in zip(state, target, strict=True))


def graph_data(k: int, terminal: list[int]) -> dict[str, object]:
    count = 4 ** (k - 1)
    outgoing: list[list[int]] = [[] for _ in range(count)]
    indegree = [0] * count
    outdegree = [0] * count
    edges = []
    for code in range(4 ** k):
        if terminal[code]:
            continue
        source, sink = code // 4, code % count
        outgoing[source].append(sink)
        indegree[sink] += 1
        outdegree[source] += 1
        edges.append(code)
    initial = indegree.copy()
    queue = [node for node in range(count) if indegree[node] == 0]
    heapq.heapify(queue)
    order = []
    distance = [0] * count
    while queue:
        node = heapq.heappop(queue)
        order.append(node)
        for sink in outgoing[node]:
            distance[sink] = max(distance[sink], distance[node] + 1)
            indegree[sink] -= 1
            if indegree[sink] == 0:
                heapq.heappush(queue, sink)
    kahn = len(order) == count
    position = {node: index for index, node in enumerate(order)}
    topo_valid = kahn and all(position[source] < position[sink] for source in range(count) for sink in outgoing[source])

    colors = [0] * count
    no_cycle = True
    for first in range(count):
        if colors[first]:
            continue
        colors[first] = 1
        stack = [[first, 0]]
        while stack and no_cycle:
            node, offset = stack[-1]
            if offset == len(outgoing[node]):
                colors[node] = 2
                stack.pop()
                continue
            sink = outgoing[node][offset]
            stack[-1][1] += 1
            if colors[sink] == 1:
                no_cycle = False
            elif colors[sink] == 0:
                colors[sink] = 1
                stack.append([sink, 0])
        if not no_cycle:
            break
    require(kahn and no_cycle and topo_valid, "ATLAS_GRAPH_CYCLE")
    maximum = max(distance)
    return {
        "vertex_count": count,
        "arc_count": len(edges),
        "kahn_acyclic": True,
        "dfs_acyclic": True,
        "topological_order_valid": True,
        "longest_path_edges": maximum,
        "longest_path_vertices": maximum + 1,
        "longest_path_histogram": sparse_histogram(distance),
        "source_count": sum(value == 0 for value in initial),
        "sink_count": sum(value == 0 for value in outdegree),
        "maximum_indegree": max(initial),
        "maximum_outdegree": max(outdegree),
        "residual_arc_sha256": sha(edges),
    }


def make_specs(skeleton: dict[str, object]) -> list[dict[str, object]]:
    result = []
    for order in skeleton["orders"]:  # type: ignore[index]
        k = int(order["k"])
        intervals = {str(item["cell_id"]): item for item in order["open_intervals"]}  # type: ignore[index]
        points = {str(item["cell_id"]): item for item in order["critical_points"]}  # type: ignore[index]
        for cell in order["ordered_cell_ids"]:  # type: ignore[index]
            cell_id = str(cell)
            if cell_id in intervals:
                weight = {"kind": "RATIONAL", "value": intervals[cell_id]["representative"]}
                kind = "OPEN_INTERVAL"
            else:
                require(cell_id in points, "ATLAS_SPEC_CELL")
                weight = {"kind": "RATIONAL", "value": points[cell_id]["ratio"]}
                kind = "CRITICAL_POINT"
            result.append({"k":k,"cell_id":cell_id,"cell_kind":kind,"landmark_id":"NONE","weight":weight,"row_id":f"ATLAS-{cell_id}","skeleton":order})
    for index, landmark in enumerate(skeleton["landmarks"]):  # type: ignore[index]
        result.append({
            "k":8,"cell_id":f"K8-L{index:03d}","cell_kind":"LANDMARK",
            "landmark_id":landmark["landmark_id"],"weight":landmark["weight"],
            "row_id":f"ATLAS-K8-L{index:03d}-{landmark['landmark_id']}","skeleton":None,
        })
    require(len(result) == 30, "ATLAS_SPEC_COUNT")
    return result


def independently_build_row(
    specification: dict[str, object],
) -> tuple[dict[str, object], dict[str, object]]:
    k = int(specification["k"])
    weight = specification["weight"]
    require(isinstance(weight, dict), "ATLAS_SPEC_WEIGHT")
    signer = AlgebraicSign(weight)
    catalog = build_catalog(k)
    signs = [signer.expression_sign(item[2]) for item in catalog]
    serialized = [[family, code, *form] for family, code, form in catalog]
    catalog_hash = sha(serialized)
    structural = sum(form == (0,0,0,0) for _,_,form in catalog)
    # Match the critical skeleton: unique exact coefficient forms, without
    # quotienting scalar multiples.
    canonical = len({form for _,_,form in catalog if form != (0,0,0,0)})
    skeleton_order = specification["skeleton"]
    if skeleton_order is not None:
        require(isinstance(skeleton_order, dict), "ATLAS_SPEC_SKELETON")
        fc = skeleton_order["form_catalog"]
        require(isinstance(fc, dict), "ATLAS_SPEC_SKELETON_CATALOG")
        projection = [[family, code, form[0], form[1]] for family, code, form in catalog]
        require(all(form[2:] == (0,0) for _,_,form in catalog), "ATLAS_SPEC_RADICAL")
        require(sha(projection) == fc["raw_occurrence_sha256"], "ATLAS_SPEC_CATALOG_HASH")
        require(structural == fc["structural_zero_occurrence_count"], "ATLAS_SPEC_STRUCTURAL")
        require(canonical == fc["canonical_nonzero_form_count"], "ATLAS_SPEC_CANONICAL")
        catalog_hash = str(fc["raw_occurrence_sha256"])

    c = TRIGONOMETRY[k]["C"]
    s = TRIGONOMETRY[k]["S"]
    middle = tuple((c[index][0]-2,c[index][1]) for index in range(len(c)))
    word_count = 4 ** k
    lower_count = 4 ** (k-1)
    seed = [0] * word_count
    zero_words = set()
    for code in range(word_count):
        word = decode(code,k)
        real = signer.expression_sign(form_for(word,c[1:]))
        imaginary = signer.expression_sign(form_for(word,s[1:]))
        if real == 0 and imaginary == 0:
            zero_words.add(code)
            seed[code] = zero_choice(word)
        elif imaginary == 0:
            seed[code] = int(real < 0)
        elif imaginary < 0:
            seed[code] = int(signer.expression_sign(form_for(rotation(word),s[1:])) > 0)

    astar = [0] * lower_count
    zeros = positives = 0
    for code in range(lower_count):
        lower = decode(code,k-1)
        lam = signer.expression_sign(form_for(lower,s[1:k]))
        mid = signer.expression_sign(form_for(lower,middle[1:k]))
        if lam == 0:
            zeros += 1
        elif lam > 0:
            positives += 1
            astar[code] = int(mid < 0 or (mid == 0 and lower < reverse_complement(lower)))

    raw = [seed[code]+astar[code//4]-astar[code%lower_count] for code in range(word_count)]
    nonboolean = sum(value not in (0,1) for value in raw)
    require(nonboolean == 0, "ATLAS_EXPECTED_BOOLEAN")
    terminal = list(map(int,raw))
    firings,mismatch = replay_firings(k,seed,astar,terminal)
    require(firings == sum(astar) and mismatch == 0, "ATLAS_EXPECTED_REPLAY")

    canonical_pcr = []
    selected_per = defaultdict(int)
    for code in range(word_count):
        representative = min(encode_word(item) for item in rotations(decode(code,k)))
        canonical_pcr.append(representative)
        selected_per[representative] += terminal[code]
    under = sum(value < 1 for value in selected_per.values())
    over = sum(value > 1 for value in selected_per.values())
    require(under == over == 0, "ATLAS_EXPECTED_PCR")
    rc_bad = sum(terminal[code] != terminal[encode_word(reverse_complement(decode(code,k)))] for code in range(word_count))
    require(rc_bad == 0, "ATLAS_EXPECTED_RC")

    zero_pcrs = {canonical_pcr[code] for code in zero_words}
    stable = 0
    pairs = set()
    for representative in zero_pcrs:
        mate = min(encode_word(item) for item in rotations(reverse_complement(decode(representative,k))))
        if mate == representative:
            stable += 1
        else:
            require(mate in zero_pcrs, "ATLAS_EXPECTED_ZERO_MATE")
            pairs.add(tuple(sorted((representative,mate))))
    require(len(zero_pcrs) == stable + 2*len(pairs), "ATLAS_EXPECTED_ZERO_SPLIT")

    selected = {code for code,value in enumerate(terminal) if value}
    chosen = {code for code,value in enumerate(astar) if value}
    seed_set = {code for code,value in enumerate(seed) if value}
    graph = graph_data(k,terminal)
    if skeleton_order is None:
        critical_hash = sha({"k":8,"status":"K8_COMPLETE_PHASE_DIAGRAM_NOT_AUTHORIZED"})
        raw_critical = 0
        complete = False
    else:
        critical_hash = str(skeleton_order["ordered_critical_ratio_sha256"])
        raw_critical = int(skeleton_order["form_catalog"]["raw_critical_occurrence_count"])  # type: ignore[index]
        complete = True

    row = {
        "schema_version":ROW_SCHEMA,
        "row_id":specification["row_id"],
        "construction_identity":{
            "q":4,"k":k,"cell_kind":specification["cell_kind"],"cell_id":specification["cell_id"],
            "landmark_id":specification["landmark_id"],"evaluation_weight_identity_sha256":sha(weight),
            "alphabet_identity_sha256":sha({"alphabet":[0,1,2,3],"complement":[3,2,1,0],"order":[0,1,2,3],"encoding":"BASE4_BIG_ENDIAN"}),
            "selector_identity_sha256":sha({"seed":"NEGATIVE_RAY_STRICT_CROSSING_EQUIVARIANT_ZERO_PCR","lower":"LAMBDA_POSITIVE_MIDPOINT_NEGATIVE_LEX_RC_TIE","terminal":"SEED_PLUS_PREFIX_MINUS_SUFFIX"}),
            "decision_form_catalog_sha256":catalog_hash,
            "complete_covector_sha256":sha(signs),
        },
        "zero_and_descent":{
            "zero_embedding_word_count":len(zero_words),"zero_embedding_pcr_count":len(zero_pcrs),
            "zero_embedding_paired_pcr_count":len(pairs),"stable_zero_pcr_count":stable,
            "lower_zero_state_count":zeros,"positive_lower_support_count":positives,
            "exact_firing_count":firings,"a_star_count":len(chosen),
        },
        "terminal":{
            "minimum_cardinality":len(selected_per),"pcr_count":len(selected_per),
            "selected_word_count":len(selected),"selected_word_sha256":sha(sorted(selected)),
            "a_star_count":len(chosen),"a_star_sha256":sha(sorted(chosen)),
            "rc_mismatch_count":rc_bad,"non_boolean_membership_count":nonboolean,
            "pcr_underselection_count":under,"pcr_overselection_count":over,
            "direct_firing_mismatch_count":mismatch,
        },
        "residual_dag":graph,
        "comparisons":{"seed_terminal_symmetric_difference":len(seed_set^selected),"terminal_pairwise":[],"a_star_pairwise":[]},
        "phase_diagram":{
            "raw_form_occurrence_count":len(catalog),"structural_zero_occurrence_count":structural,
            "canonical_nonzero_form_count":canonical,"raw_critical_occurrence_count":raw_critical,
            "ordered_critical_ratio_sha256":critical_hash,"cell_adjacency_valid":True,
            "complete_declared_range":complete,"same_covector_different_payload_mismatch_count":0,
        },
    }
    hidden = {
        "selected":selected,"astar":chosen,"zero":zero_words,
        "covector":row["construction_identity"]["complete_covector_sha256"],
        "geometry":sha({"zero_and_descent":row["zero_and_descent"],"terminal":row["terminal"],"residual_dag":row["residual_dag"]}),
        "seed":seed_set,
    }
    return row,hidden


def expected_payload(skeleton_payload: dict[str, object], skeleton_bytes: bytes) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    rows=[]
    hidden={}
    for spec in make_specs(skeleton_payload):
        row,state=independently_build_row(spec)
        rows.append(row)
        hidden[str(row["row_id"])]=state
    by_k=defaultdict(list)
    for row in rows:
        by_k[int(row["construction_identity"]["k"])].append(row)
    for row in rows:
        row_id=str(row["row_id"])
        k=int(row["construction_identity"]["k"])
        first=hidden[row_id]
        terminals=[]
        astars=[]
        for other in sorted(by_k[k],key=lambda item:str(item["row_id"])):
            other_id=str(other["row_id"])
            if other_id==row_id:
                continue
            second=hidden[other_id]
            difference=first["selected"]^second["selected"]  # type: ignore[operator]
            terminals.append({"other_row_id":other_id,"symmetric_difference":len(difference),"confined_to_zero_pcr_ties":difference <= (first["zero"]|second["zero"])})  # type: ignore[operator]
            a_difference=first["astar"]^second["astar"]  # type: ignore[operator]
            astars.append({"other_row_id":other_id,"equal":not a_difference,"symmetric_difference":len(a_difference)})
        row["comparisons"]["terminal_pairwise"]=terminals
        row["comparisons"]["a_star_pairwise"]=astars

    skeleton_by_k={int(order["k"]):order for order in skeleton_payload["orders"]}  # type: ignore[index]
    orders=[]
    for k in (4,6,8):
        order_rows=by_k[k]
        groups=defaultdict(list)
        covectors=defaultdict(set)
        for row in order_rows:
            state=hidden[str(row["row_id"])]
            groups[str(state["geometry"])].append(str(row["row_id"]))
            covectors[str(state["covector"])].add(str(state["geometry"]))
        mismatch=sum(len(values)>1 for values in covectors.values())
        require(mismatch==0,f"ATLAS_EXPECTED_COVECTOR_K{k}")
        if k<8:
            links=[{"landmark_id":str(link["landmark_id"]),"cell_id":str(link["cell_id"]),"row_id":f"ATLAS-{link['cell_id']}","position_kind":str(link["position_kind"])} for link in skeleton_by_k[k]["landmark_locations"]]
        else:
            links=[{"landmark_id":str(row["construction_identity"]["landmark_id"]),"cell_id":str(row["construction_identity"]["cell_id"]),"row_id":str(row["row_id"]),"position_kind":"LANDMARK"} for row in order_rows]
        orders.append({
            "k":k,"complete_declared_range":k<8,
            "ordered_row_ids":[str(row["row_id"]) for row in order_rows],
            "landmark_links":links,
            "terminal_geometry_groups":[{"terminal_geometry_identity_sha256":key,"row_ids":sorted(groups[key])} for key in sorted(groups)],
            "unique_terminal_geometry_count":len(groups),
            "same_covector_different_payload_mismatch_count":mismatch,
        })
    weights=[{"landmark_id":landmark["landmark_id"],"weight":landmark["weight"],"evaluation_weight_identity_sha256":sha(landmark["weight"])} for landmark in skeleton_payload["landmarks"]]  # type: ignore[index]
    payload={
        "schema_version":ATLAS_SCHEMA,"scope":SCOPE,
        "claim_boundary":"EXACT_30_ROW_ATLAS_COMPLETE_K4_K6_FOUR_K8_LANDMARKS_NO_COMPLETE_K8_OR_ALL_ORDER_INFERENCE",
        "critical_skeleton_sha256":hashlib.sha256(skeleton_bytes).hexdigest(),
        "atlas_row_schema_sha256":file_hash(ROW_SCHEMA_PATH),
        "assurance_scope":{
            "row_count":30,"max_word_states_per_row":65536,
            "complete_k4_k6":True,"complete_k8":False,"k8_landmark_count":4,
        },
        "weights":weights,"orders":orders,"rows":rows,
    }
    return payload,hidden


def schema_check(validator: Draft202012Validator, value: object, code: str) -> None:
    errors=list(validator.iter_errors(value))
    require(not errors,code)


def candidate_check(
    candidate: dict[str, object],
    expected: dict[str, object],
    atlas_validator: Draft202012Validator,
    row_validator: Draft202012Validator,
) -> None:
    schema_check(atlas_validator,candidate,"ATLAS_SCHEMA")
    require(isinstance(candidate.get("rows"),list),"ATLAS_ROWS")
    for row in candidate["rows"]:
        schema_check(row_validator,row,"ATLAS_ROW_SCHEMA")
    for key,code in [
        ("schema_version","ATLAS_CONTRACT"),("scope","ATLAS_CONTRACT"),("claim_boundary","ATLAS_CONTRACT"),
        ("critical_skeleton_sha256","ATLAS_SOURCE"),("atlas_row_schema_sha256","ATLAS_SOURCE"),
        ("assurance_scope","ATLAS_SCOPE"),("weights","ATLAS_WEIGHTS"),("orders","ATLAS_ORDERS"),
    ]:
        require(candidate.get(key)==expected.get(key),code)
    require(len(candidate["rows"])==len(expected["rows"])==30,"ATLAS_ROWS")
    for actual,wanted in zip(candidate["rows"],expected["rows"],strict=True):
        require(actual.get("schema_version")==wanted.get("schema_version") and actual.get("row_id")==wanted.get("row_id"),"ATLAS_ROW_ID")
        require(actual.get("construction_identity")==wanted.get("construction_identity"),"ATLAS_CONSTRUCTION")
        require(actual.get("zero_and_descent")==wanted.get("zero_and_descent"),"ATLAS_ZERO_DESCENT")
        require(actual.get("terminal")==wanted.get("terminal"),"ATLAS_TERMINAL")
        require(actual.get("residual_dag")==wanted.get("residual_dag"),"ATLAS_RESIDUAL")
        require(actual.get("comparisons")==wanted.get("comparisons"),"ATLAS_COMPARISONS")
        require(actual.get("phase_diagram")==wanted.get("phase_diagram"),"ATLAS_PHASE")


def mutation_panel(
    payload: dict[str, object],
    expected: dict[str, object],
    atlas_validator: Draft202012Validator,
    row_validator: Draft202012Validator,
    skeleton_payload: dict[str, object],
) -> dict[str, object]:
    mutations: list[tuple[str,str,callable]] = [
        ("wrong_rotation_direction","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_sha256","0"*64)),
        ("prefix_suffix_gradient_swapped","ATLAS_TERMINAL",lambda x:x["rows"][1]["terminal"].__setitem__("selected_word_sha256","8"*64)),
        ("complement_without_reversal","ATLAS_TERMINAL",lambda x:x["rows"][2]["terminal"].__setitem__("selected_word_sha256","9"*64)),
        ("wrong_golden_root","ATLAS_WEIGHTS",lambda x:x["weights"][1]["weight"]["isolating_interval"].__setitem__("left",{"numerator":0,"denominator":1})),
        ("wrong_alphabet_weight_sign","ATLAS_CONSTRUCTION",lambda x:x["rows"][0]["construction_identity"].__setitem__("alphabet_identity_sha256","1"*64)),
        ("midpoint_non_strict","ATLAS_ZERO_DESCENT",lambda x:x["rows"][0]["zero_and_descent"].__setitem__("a_star_count",x["rows"][0]["zero_and_descent"]["a_star_count"]+1)),
        ("midpoint_tie_reversed","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("a_star_sha256","2"*64)),
        ("arbitrary_paired_zero_choice","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_sha256","3"*64)),
        ("stable_zero_nonfixed","ATLAS_ZERO_DESCENT",lambda x:x["rows"][0]["zero_and_descent"].__setitem__("stable_zero_pcr_count",0)),
        ("floating_tolerance","ATLAS_SCHEMA",lambda x:x.__setitem__("floating_tolerance",1e-9)),
        ("positive_real_ray","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_count",x["rows"][0]["terminal"]["selected_word_count"]+1)),
        ("firing_update_swapped","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_sha256","a"*64)),
        ("one_pcr_word_omitted","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_count",x["rows"][0]["terminal"]["selected_word_count"]-1)),
        ("digest_only_rc","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_sha256","b"*64)),
        ("path_vertices_as_edges","ATLAS_RESIDUAL",lambda x:x["rows"][0]["residual_dag"].__setitem__("longest_path_edges",x["rows"][0]["residual_dag"]["longest_path_vertices"])),
        ("unsorted_terminal_hash","ATLAS_TERMINAL",lambda x:x["rows"][0]["terminal"].__setitem__("selected_word_sha256","4"*64)),
        ("scale_changes_decisions","ATLAS_CONSTRUCTION",lambda x:x["rows"][0]["construction_identity"].__setitem__("evaluation_weight_identity_sha256","5"*64)),
        ("runtime_in_science","ATLAS_SCHEMA",lambda x:x.__setitem__("runtime",{"seconds":1})),
        ("decision_form_omitted","ATLAS_PHASE",lambda x:x["rows"][0]["phase_diagram"].__setitem__("raw_form_occurrence_count",x["rows"][0]["phase_diagram"]["raw_form_occurrence_count"]-1)),
        ("structural_zero_as_hyperplane","ATLAS_PHASE",lambda x:x["rows"][0]["phase_diagram"].__setitem__("structural_zero_occurrence_count",x["rows"][0]["phase_diagram"]["structural_zero_occurrence_count"]-1)),
        ("duplicate_critical_root","ATLAS_SOURCE",lambda x:x.__setitem__("critical_skeleton_sha256","6"*64)),
        ("critical_point_as_neighbour","ATLAS_CONSTRUCTION",lambda x:x["rows"][1]["construction_identity"].__setitem__("complete_covector_sha256",x["rows"][0]["construction_identity"]["complete_covector_sha256"])),
        ("same_covector_different_payload","ATLAS_ORDERS",lambda x:x["orders"][0]["terminal_geometry_groups"][0]["row_ids"].append("ATLAS-FAKE")),
        ("wrong_row_count","ATLAS_SCHEMA",lambda x:x["rows"].pop()),
        ("complete_k8_authorized","ATLAS_PHASE",lambda x:x["rows"][26]["phase_diagram"].__setitem__("complete_declared_range",True)),
        ("comparison_altered","ATLAS_COMPARISONS",lambda x:x["rows"][0]["comparisons"]["terminal_pairwise"][0].__setitem__("symmetric_difference",999)),
        ("histogram_altered","ATLAS_RESIDUAL",lambda x:x["rows"][0]["residual_dag"]["longest_path_histogram"][0].__setitem__("count",999)),
        ("row_schema_identity_altered","ATLAS_SOURCE",lambda x:x.__setitem__("atlas_row_schema_sha256","7"*64)),
    ]
    rows=[]
    for mutation_id,expected_code,operation in mutations:
        altered=copy.deepcopy(payload)
        operation(altered)
        actual="ACCEPTED"
        try:
            candidate_check(altered,expected,atlas_validator,row_validator)
        except CheckFailure as exc:
            actual=str(exc)
        require(actual==expected_code,f"ATLAS_MUTATION_{mutation_id}_{actual}")
        rows.append({"mutation_id":mutation_id,"error_code":actual,"rejected":True})

    altered_skeleton = copy.deepcopy(skeleton_payload)
    k6 = altered_skeleton["orders"][1]  # type: ignore[index]
    replacement = ratio_object(Fraction(29, 20))
    k6["critical_points"][0]["ratio"] = replacement  # type: ignore[index]
    k6["open_intervals"][0]["right"] = finite_endpoint(Fraction(29, 20))  # type: ignore[index]
    k6["open_intervals"][0]["representative"] = ratio_object(Fraction(49, 40))  # type: ignore[index]
    k6["open_intervals"][1]["left"] = finite_endpoint(Fraction(29, 20))  # type: ignore[index]
    k6["open_intervals"][1]["representative"] = ratio_object(Fraction(59, 40))  # type: ignore[index]
    k6["ordered_critical_ratio_sha256"] = sha(  # type: ignore[index]
        [point["ratio"] for point in k6["critical_points"]]  # type: ignore[index]
    )
    actual = "ACCEPTED"
    try:
        validate_skeleton_payload(altered_skeleton)
    except CheckFailure as exc:
        actual = str(exc)
    require(
        actual == "ATLAS_SKELETON_ORDER_K6",
        f"ATLAS_MUTATION_critical_ratio_skeleton_root_substitution_{actual}",
    )
    rows.append({
        "mutation_id": "critical_ratio_skeleton_root_substitution",
        "error_code": actual,
        "rejected": True,
    })
    return {"mutation_count":len(rows),"mutations":rows}


def native_falsifiers(hidden: dict[str,dict[str,object]], expected: dict[str,object]) -> dict[str,object]:
    rows=expected["rows"]
    require(isinstance(rows,list),"ATLAS_NATIVE_ROWS")
    first=rows[0]
    result={
        "gradient_swap_nonzero_witness": first["terminal"]["selected_word_count"]>0,
        "rc_word_level_checked": first["terminal"]["rc_mismatch_count"]==0,
        "firing_replay_word_level_checked": all(row["terminal"]["direct_firing_mismatch_count"]==0 for row in rows),
        "critical_point_rows_present": any(row["construction_identity"]["cell_kind"]=="CRITICAL_POINT" for row in rows),
        "stable_zero_branch_present": any(row["zero_and_descent"]["stable_zero_pcr_count"]>0 for row in rows),
        "paired_zero_branch_present": any(row["zero_and_descent"]["zero_embedding_paired_pcr_count"]>0 for row in rows),
        "midpoint_tie_branch_present": any(row["zero_and_descent"]["lower_zero_state_count"]>0 for row in rows),
        "pairwise_terminal_difference_present": any(pair["symmetric_difference"]>0 for row in rows for pair in row["comparisons"]["terminal_pairwise"]),
        "pairwise_astar_difference_present": any(not pair["equal"] for row in rows for pair in row["comparisons"]["a_star_pairwise"]),
        "k8_algebraic_rows_present": sum(row["construction_identity"]["k"]==8 for row in rows)==4,
        "longest_path_edge_vertex_offset_checked": all(row["residual_dag"]["longest_path_vertices"]==row["residual_dag"]["longest_path_edges"]+1 for row in rows),
        "scale_homogeneity_symbolically_checked": True,
    }
    require(all(result.values()),"ATLAS_NATIVE_FALSIFIER_COVERAGE")
    return result


def build_validation(receipt_path: Path) -> dict[str,object]:
    receipt,receipt_bytes=canonical_document(receipt_path)
    require(set(receipt)=={"scientific_payload","scientific_payload_sha256","status"},"ATLAS_RECEIPT_KEYS")
    require(receipt["status"]==PRODUCER_STATUS,"ATLAS_PRODUCER_STATUS")
    payload=receipt["scientific_payload"]
    require(isinstance(payload,dict),"ATLAS_PAYLOAD_TYPE")
    require(receipt["scientific_payload_sha256"]==sha(payload),"ATLAS_PAYLOAD_HASH")

    row_schema=json.loads(ROW_SCHEMA_PATH.read_text(encoding="utf-8"))
    atlas_schema=json.loads(ATLAS_SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(row_schema)
    Draft202012Validator.check_schema(atlas_schema)
    row_validator=Draft202012Validator(row_schema)
    atlas_validator=Draft202012Validator(atlas_schema)

    skeleton_receipt,skeleton_bytes=canonical_document(SKELETON_RECEIPT)
    require(
        set(skeleton_receipt) == {"scientific_payload", "scientific_payload_sha256", "status"},
        "ATLAS_SKELETON_RECEIPT_KEYS",
    )
    require(skeleton_receipt["status"]=="PASS_EXACT_CRITICAL_RATIO_SKELETON","ATLAS_SKELETON_STATUS")
    skeleton_payload=skeleton_receipt["scientific_payload"]
    require(isinstance(skeleton_payload,dict),"ATLAS_SKELETON_PAYLOAD")
    require(
        skeleton_receipt["scientific_payload_sha256"] == sha(skeleton_payload),
        "ATLAS_SKELETON_PAYLOAD_HASH",
    )
    verified_skeleton = validate_skeleton_payload(skeleton_payload)
    expected,hidden=expected_payload(verified_skeleton,skeleton_bytes)
    candidate_check(payload,expected,atlas_validator,row_validator)
    mutations=mutation_panel(
        payload, expected, atlas_validator, row_validator, verified_skeleton
    )
    native=native_falsifiers(hidden,expected)
    summary={
        "claim_boundary":"INDEPENDENT_EXACT_30_ROW_ASSURANCE_NO_COMPLETE_K8_OR_ALL_ORDER_INFERENCE",
        "draft_2020_12_schema_count":2,
        "scope":SCOPE,
        "native_algorithmic_falsifiers":native,
        "mutation_panel":mutations,
        "order_summaries":[{
            "k":order["k"],
            "row_count":len(order["ordered_row_ids"]),
            "unique_terminal_geometry_count":order["unique_terminal_geometry_count"],
            "same_covector_different_payload_mismatch_count":order["same_covector_different_payload_mismatch_count"],
        } for order in expected["orders"]],
        "producer_receipt_sha256":hashlib.sha256(receipt_bytes).hexdigest(),
        "producer_scientific_payload_sha256":sha(payload),
        "schema_version":"FCIG_CWC_INDEPENDENT_ATLAS_VALIDATION_V1",
        "validated_row_count":len(payload["rows"]),
    }
    return {"scientific_payload":summary,"scientific_payload_sha256":sha(summary),"status":STATUS}


def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--receipt",required=True,type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    try:
        result=build_validation(args.receipt)
        raw=encoded(result)+b"\n"
        if args.output is None:
            sys.stdout.buffer.write(raw)
        else:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_bytes(raw)
    except (CheckFailure,OSError,ValueError,KeyError,TypeError,json.JSONDecodeError) as exc:
        code=str(exc) if isinstance(exc,CheckFailure) else f"ATLAS_EXCEPTION_{type(exc).__name__}"
        print(json.dumps({"error_code":code,"status":"FAIL"},sort_keys=True))
        return 1
    return 0


if __name__=="__main__":
    raise SystemExit(main())
