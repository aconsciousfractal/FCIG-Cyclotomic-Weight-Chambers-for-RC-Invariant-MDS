#!/usr/bin/env python3
"""Produce the exact bounded DNA terminal atlas shipped with the article.

The producer enumerates words with itertools.product.  It uses only rational
integer-polynomial arithmetic and exact isolating intervals; no floating point
decision is made.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import heapq
import itertools
import json
import math
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA_VERSION = "FCIG_CWC_TERMINAL_ATLAS_V1"
ROW_SCHEMA_VERSION = "FCIG_CWC_ATLAS_ROW_V1"
STATUS = "PASS_EXACT_SMALL_DNA_TERMINAL_ATLAS_PRODUCER"
SCOPE = "EXACT_SMALL_DNA_TERMINAL_ATLAS"
CLAIM_BOUNDARY = (
    "EXACT_30_ROW_ATLAS_COMPLETE_K4_K6_FOUR_K8_LANDMARKS_"
    "NO_COMPLETE_K8_OR_ALL_ORDER_INFERENCE"
)
ALPHABET = (0, 1, 2, 3)
COMPLEMENT = (3, 2, 1, 0)
FAMILY_ORDER = (
    "SEED_REAL",
    "SEED_IMAG",
    "LOWER_LAMBDA",
    "LOWER_MIDPOINT",
)
CAPSULE = Path(__file__).resolve().parents[1]
SKELETON_RECEIPT = CAPSULE / "inputs" / "critical_ratio_skeleton.json"
ROW_SCHEMA_PATH = CAPSULE / "schemas" / "atlas_row.schema.json"
MAX_SECONDS = 7200

# Coefficients are pairs (constant, sqrt(2)).
TRIG = {
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


class ProductionFailure(RuntimeError):
    pass


def require(condition: bool, code: str) -> None:
    if not condition:
        raise ProductionFailure(code)


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_canonical(path: Path) -> tuple[dict[str, object], bytes]:
    raw = path.read_bytes()
    require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "ATLASP_JSON_LF")
    value = json.loads(raw)
    require(isinstance(value, dict), "ATLASP_JSON_OBJECT")
    require(raw == canonical_bytes(value) + b"\n", "ATLASP_JSON_CANONICAL")
    return value, raw


def frac(data: dict[str, object]) -> Fraction:
    value = Fraction(int(data["numerator"]), int(data["denominator"]))
    require(value.denominator > 0, "ATLASP_FRACTION_DENOMINATOR")
    return value


def fraction_object(value: Fraction) -> dict[str, int]:
    return {"denominator": value.denominator, "numerator": value.numerator}


def trim(poly: Sequence[Fraction]) -> tuple[Fraction, ...]:
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [Fraction(0)])


def poly_eval(poly: Sequence[Fraction | int], x: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * x + Fraction(coefficient)
    return result


def poly_remainder(poly: Sequence[int], modulus: Sequence[int]) -> tuple[Fraction, ...]:
    work = list(map(Fraction, trim(tuple(map(Fraction, poly)))))
    mod = list(map(Fraction, trim(tuple(map(Fraction, modulus)))))
    require(mod[-1] != 0, "ATLASP_ZERO_MINPOLY")
    while len(work) >= len(mod) and any(work):
        factor = work[-1] / mod[-1]
        shift = len(work) - len(mod)
        for index, value in enumerate(mod):
            work[index + shift] -= factor * value
        work = list(trim(work))
    return trim(work)


def interval_multiply(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    products = (
        left[0] * right[0], left[0] * right[1],
        left[1] * right[0], left[1] * right[1],
    )
    return min(products), max(products)


def interval_poly(
    poly: Sequence[Fraction], interval: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    result = (Fraction(0), Fraction(0))
    for coefficient in reversed(poly):
        result = interval_multiply(result, interval)
        result = (result[0] + coefficient, result[1] + coefficient)
    return result


class ExactRatio:
    def __init__(self, weight: dict[str, object]) -> None:
        self.weight = copy.deepcopy(weight)
        self.kind = str(weight["kind"])
        self.rational: Fraction | None = None
        self.polynomial: tuple[int, ...] | None = None
        self.interval: tuple[Fraction, Fraction] | None = None
        if self.kind == "RATIONAL":
            self.rational = frac(weight["value"])  # type: ignore[arg-type]
            require(self.rational > 1, "ATLASP_RATIO_DOMAIN")
        else:
            require(self.kind == "REAL_ALGEBRAIC", "ATLASP_WEIGHT_KIND")
            self.polynomial = tuple(int(x) for x in weight["polynomial_coefficients_ascending"])  # type: ignore[index]
            isolating = weight["isolating_interval"]
            require(isinstance(isolating, dict), "ATLASP_ISOLATING_TYPE")
            left = frac(isolating["left"])  # type: ignore[arg-type]
            right = frac(isolating["right"])  # type: ignore[arg-type]
            require(Fraction(1) < left < right, "ATLASP_ISOLATING_ORDER")
            require(poly_eval(self.polynomial, left) * poly_eval(self.polynomial, right) < 0, "ATLASP_ISOLATING_BRACKET")
            self.interval = (left, right)

    def sign_polynomial(self, coefficients: Sequence[int]) -> int:
        if self.rational is not None:
            value = poly_eval(coefficients, self.rational)
            return (value > 0) - (value < 0)
        require(self.polynomial is not None and self.interval is not None, "ATLASP_ALGEBRAIC_STATE")
        remainder = poly_remainder(coefficients, self.polynomial)
        if all(value == 0 for value in remainder):
            return 0
        left, right = self.interval
        p_left = poly_eval(self.polynomial, left)
        for _ in range(1024):
            lower, upper = interval_poly(remainder, (left, right))
            if lower > 0:
                return 1
            if upper < 0:
                return -1
            midpoint = (left + right) / 2
            p_mid = poly_eval(self.polynomial, midpoint)
            if p_mid == 0:
                value = poly_eval(remainder, midpoint)
                return (value > 0) - (value < 0)
            if p_left * p_mid < 0:
                right = midpoint
            else:
                left = midpoint
                p_left = p_mid
        raise ProductionFailure("ATLASP_ALGEBRAIC_SIGN_REFINEMENT")

    def sign_expression(self, expression: tuple[int, int, int, int]) -> int:
        a, b, c, d = expression
        sign_n = self.sign_polynomial((b, a))
        sign_m = self.sign_polynomial((d, c))
        if sign_n == 0 and sign_m == 0:
            return 0
        if sign_m == 0:
            return sign_n
        if sign_n == 0:
            return sign_m
        if sign_n == sign_m:
            return sign_n
        q = (b * b - 2 * d * d, 2 * a * b - 4 * c * d, a * a - 2 * c * c)
        sign_q = self.sign_polynomial(q)
        if sign_q == 0:
            return 0
        return sign_n if sign_q > 0 else sign_m


def decode_word(code: int, length: int) -> tuple[int, ...]:
    result = [0] * length
    value = code
    for index in range(length - 1, -1, -1):
        result[index] = value & 3
        value >>= 2
    require(value == 0, "ATLASP_DECODE_RANGE")
    return tuple(result)


def word_code(word: Sequence[int]) -> int:
    result = 0
    for symbol in word:
        result = 4 * result + symbol
    return result


def words(length: int) -> Iterable[tuple[int, ...]]:
    return itertools.product(ALPHABET, repeat=length)


def rc(word: Sequence[int]) -> tuple[int, ...]:
    return tuple(COMPLEMENT[symbol] for symbol in reversed(word))


def rotate_left(word: Sequence[int]) -> tuple[int, ...]:
    data = tuple(word)
    return data[1:] + data[:1]


def orbit(word: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    data = tuple(word)
    result: list[tuple[int, ...]] = []
    current = data
    while current not in result:
        result.append(current)
        current = rotate_left(current)
    return tuple(result)


def family_coefficients(k: int, family: str) -> tuple[tuple[int, int], ...]:
    c = TRIG[k]["C"]
    s = TRIG[k]["S"]
    if family == "SEED_REAL":
        return tuple(c[index] for index in range(1, k + 1))
    if family == "SEED_IMAG":
        return tuple(s[index] for index in range(1, k + 1))
    if family == "LOWER_LAMBDA":
        return tuple(s[index] for index in range(1, k))
    if family == "LOWER_MIDPOINT":
        return tuple((c[index][0] - 2, c[index][1]) for index in range(1, k))
    raise ProductionFailure("ATLASP_FAMILY")


def expression(word: Sequence[int], coefficients: Sequence[tuple[int, int]]) -> tuple[int, int, int, int]:
    require(len(word) == len(coefficients), "ATLASP_EXPRESSION_LENGTH")
    a = b = c = d = 0
    for symbol, (constant, radical) in zip(word, coefficients, strict=True):
        if symbol == 0:
            a -= constant
            c -= radical
        elif symbol == 1:
            b -= constant
            d -= radical
        elif symbol == 2:
            b += constant
            d += radical
        elif symbol == 3:
            a += constant
            c += radical
        else:
            raise ProductionFailure("ATLASP_SYMBOL")
    return a, b, c, d


def normalize_expression(value: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    if value == (0, 0, 0, 0):
        return value
    divisor = 0
    for item in value:
        divisor = math.gcd(divisor, abs(item))
    require(divisor > 0, "ATLASP_NORMALIZE_GCD")
    result = tuple(item // divisor for item in value)
    first = next(item for item in result if item)
    if first < 0:
        result = tuple(-item for item in result)
    return result  # type: ignore[return-value]


def raw_catalog(k: int) -> list[tuple[str, int, tuple[int, int, int, int]]]:
    result: list[tuple[str, int, tuple[int, int, int, int]]] = []
    for family in FAMILY_ORDER:
        coefficients = family_coefficients(k, family)
        for word in words(len(coefficients)):
            result.append((family, word_code(word), expression(word, coefficients)))
    return result


def zero_seed_membership(word: tuple[int, ...]) -> int:
    own = orbit(word)
    mate_word = rc(word)
    if mate_word in set(own):
        fixed = [candidate for candidate in own if rc(candidate) == candidate]
        require(len(fixed) == 2, "ATLASP_STABLE_ZERO_FIXED_ROTATIONS")
        return int(word == min(fixed))
    mate = orbit(mate_word)
    small = min(min(own), min(mate))
    selected = {small, rc(small)}
    require(len(selected) == 2, "ATLASP_PAIRED_ZERO_SELECTION")
    return int(word in selected)


def histogram(values: Iterable[int]) -> list[dict[str, int]]:
    counts = Counter(values)
    return [{"bucket": bucket, "count": counts[bucket]} for bucket in sorted(counts) if counts[bucket] > 0]


def legal_replay(k: int, seed: list[int], astar: list[int], terminal: list[int]) -> tuple[int, int]:
    lower_count = 4 ** (k - 1)
    dependencies: list[set[int]] = [set() for _ in range(lower_count)]
    reverse: list[set[int]] = [set() for _ in range(lower_count)]
    for code, selected in enumerate(seed):
        prefix = code // 4
        suffix = code % lower_count
        if prefix != suffix and astar[prefix] and astar[suffix]:
            before, after = (prefix, suffix) if selected == 0 else (suffix, prefix)
            dependencies[after].add(before)
            reverse[before].add(after)
    indegree = [len(dependencies[node]) for node in range(lower_count)]
    heap = [node for node in range(lower_count) if astar[node] and indegree[node] == 0]
    heapq.heapify(heap)
    order: list[int] = []
    while heap:
        node = heapq.heappop(heap)
        order.append(node)
        for target in sorted(reverse[node]):
            indegree[target] -= 1
            if indegree[target] == 0:
                heapq.heappush(heap, target)
    required = sum(astar)
    require(len(order) == required, "ATLASP_FIRING_DEPENDENCY_CYCLE")
    current = seed.copy()
    for node in order:
        incoming = [symbol * lower_count + node for symbol in ALPHABET]
        outgoing = [node * 4 + symbol for symbol in ALPHABET]
        require(not set(incoming).intersection(outgoing), "ATLASP_FIRING_SELF_LOOP")
        require(all(current[code] == 1 for code in incoming), "ATLASP_FIRING_INCOMING")
        require(all(current[code] == 0 for code in outgoing), "ATLASP_FIRING_OUTGOING")
        for code in incoming:
            current[code] = 0
        for code in outgoing:
            current[code] = 1
    mismatch = sum(a != b for a, b in zip(current, terminal, strict=True))
    return len(order), mismatch


def residual_geometry(k: int, terminal: list[int]) -> dict[str, object]:
    vertices = 4 ** (k - 1)
    adjacency: list[list[int]] = [[] for _ in range(vertices)]
    indegree = [0] * vertices
    outdegree = [0] * vertices
    residual_codes: list[int] = []
    for code, selected in enumerate(terminal):
        if selected:
            continue
        source = code // 4
        target = code % vertices
        adjacency[source].append(target)
        indegree[target] += 1
        outdegree[source] += 1
        residual_codes.append(code)
    initial_indegree = indegree.copy()
    heap = [node for node in range(vertices) if indegree[node] == 0]
    heapq.heapify(heap)
    topological: list[int] = []
    longest = [0] * vertices
    while heap:
        node = heapq.heappop(heap)
        topological.append(node)
        for target in adjacency[node]:
            if longest[target] < longest[node] + 1:
                longest[target] = longest[node] + 1
            indegree[target] -= 1
            if indegree[target] == 0:
                heapq.heappush(heap, target)
    kahn = len(topological) == vertices

    color = bytearray(vertices)
    dfs_acyclic = True
    for start in range(vertices):
        if color[start]:
            continue
        color[start] = 1
        stack: list[tuple[int, int]] = [(start, 0)]
        while stack and dfs_acyclic:
            node, index = stack[-1]
            if index == len(adjacency[node]):
                color[node] = 2
                stack.pop()
                continue
            target = adjacency[node][index]
            stack[-1] = (node, index + 1)
            if color[target] == 1:
                dfs_acyclic = False
                break
            if color[target] == 0:
                color[target] = 1
                stack.append((target, 0))
        if not dfs_acyclic:
            break
    require(kahn and dfs_acyclic, "ATLASP_RESIDUAL_CYCLE")
    maximum = max(longest, default=0)
    return {
        "vertex_count": vertices,
        "arc_count": len(residual_codes),
        "kahn_acyclic": kahn,
        "dfs_acyclic": dfs_acyclic,
        "topological_order_valid": all(
            position[source] < position[target]
            for position in [{node: index for index, node in enumerate(topological)}]
            for source in range(vertices)
            for target in adjacency[source]
        ),
        "longest_path_edges": maximum,
        "longest_path_vertices": maximum + 1,
        "longest_path_histogram": histogram(longest),
        "source_count": sum(value == 0 for value in initial_indegree),
        "sink_count": sum(value == 0 for value in outdegree),
        "maximum_indegree": max(initial_indegree, default=0),
        "maximum_outdegree": max(outdegree, default=0),
        "residual_arc_sha256": digest(residual_codes),
    }


def selector_row(
    k: int,
    weight: dict[str, object],
    row_id: str,
    cell_id: str,
    cell_kind: str,
    landmark_id: str,
    skeleton_order: dict[str, object] | None,
) -> tuple[dict[str, object], dict[str, object]]:
    ratio = ExactRatio(weight)
    catalog = raw_catalog(k)
    covector = [ratio.sign_expression(item[2]) for item in catalog]
    serialized_catalog = [[family, code, *form] for family, code, form in catalog]
    catalog_sha = digest(serialized_catalog)
    structural_count = sum(form == (0, 0, 0, 0) for _, _, form in catalog)
    # The critical skeleton's canonical catalog deduplicates exactly equal forms;
    # it does not quotient positive scalar multiples.
    canonical_count = len({
        form for _, _, form in catalog if form != (0, 0, 0, 0)
    })
    if skeleton_order is not None:
        form_catalog = skeleton_order["form_catalog"]
        require(isinstance(form_catalog, dict), "ATLASP_SKELETON_FORM_CATALOG")
        projected = [[family, code, form[0], form[1]] for family, code, form in catalog]
        require(all(form[2:] == (0, 0) for _, _, form in catalog), "ATLASP_SKELETON_RADICAL_FORM")
        require(digest(projected) == form_catalog["raw_occurrence_sha256"], "ATLASP_SKELETON_RAW_CATALOG")
        require(structural_count == form_catalog["structural_zero_occurrence_count"], "ATLASP_SKELETON_STRUCTURAL_COUNT")
        require(canonical_count == form_catalog["canonical_nonzero_form_count"], "ATLASP_SKELETON_CANONICAL_COUNT")
        catalog_sha = str(form_catalog["raw_occurrence_sha256"])

    word_list = list(words(k))
    lower_list = list(words(k - 1))
    require(len(word_list) <= 65536, "ATLASP_STATE_CAP")
    c = TRIG[k]["C"]
    s = TRIG[k]["S"]
    midpoint = tuple((c[index][0] - 2, c[index][1]) for index in range(len(c)))
    seed = [0] * len(word_list)
    zero_words: set[int] = set()
    for code, word in enumerate(word_list):
        real_sign = ratio.sign_expression(expression(word, c[1:]))
        imag_sign = ratio.sign_expression(expression(word, s[1:]))
        if real_sign == 0 and imag_sign == 0:
            zero_words.add(code)
            seed[code] = zero_seed_membership(word)
        elif imag_sign == 0:
            seed[code] = int(real_sign < 0)
        elif imag_sign < 0:
            rotated_sign = ratio.sign_expression(expression(rotate_left(word), s[1:]))
            seed[code] = int(rotated_sign > 0)

    astar = [0] * len(lower_list)
    lower_zero = 0
    positive_lower = 0
    for code, lower in enumerate(lower_list):
        lambda_sign = ratio.sign_expression(expression(lower, s[1:k]))
        midpoint_sign = ratio.sign_expression(expression(lower, midpoint[1:k]))
        if lambda_sign == 0:
            lower_zero += 1
        elif lambda_sign > 0:
            positive_lower += 1
            if midpoint_sign < 0:
                astar[code] = 1
            elif midpoint_sign == 0:
                theta = rc(lower)
                require(theta != lower, "ATLASP_THETA_FIXED_LOWER")
                astar[code] = int(lower < theta)

    lower_count = 4 ** (k - 1)
    raw_terminal = [
        seed[code] + astar[code // 4] - astar[code % lower_count]
        for code in range(4 ** k)
    ]
    non_boolean = sum(value not in (0, 1) for value in raw_terminal)
    require(non_boolean == 0, "ATLASP_TERMINAL_BOOLEANITY")
    terminal = [int(value) for value in raw_terminal]
    firing_count, direct_mismatch = legal_replay(k, seed, astar, terminal)
    require(direct_mismatch == 0 and firing_count == sum(astar), "ATLASP_DIRECT_FIRING")

    canonical_pcr: list[int] = [0] * (4 ** k)
    pcr_selected: dict[int, int] = defaultdict(int)
    for code, word in enumerate(word_list):
        canonical = min(word_code(item) for item in orbit(word))
        canonical_pcr[code] = canonical
        pcr_selected[canonical] += terminal[code]
    under = sum(value < 1 for value in pcr_selected.values())
    over = sum(value > 1 for value in pcr_selected.values())
    require(under == 0 and over == 0, "ATLASP_ONE_PER_PCR")

    rc_mismatch = 0
    for code, word in enumerate(word_list):
        if terminal[code] != terminal[word_code(rc(word))]:
            rc_mismatch += 1
    require(rc_mismatch == 0, "ATLASP_LITERAL_RC")

    zero_pcrs = {canonical_pcr[code] for code in zero_words}
    stable = 0
    paired_keys: set[tuple[int, int]] = set()
    for canonical in zero_pcrs:
        word = decode_word(canonical, k)
        mate = min(word_code(item) for item in orbit(rc(word)))
        if mate == canonical:
            stable += 1
        else:
            require(mate in zero_pcrs, "ATLASP_ZERO_PCR_MATE")
            paired_keys.add(tuple(sorted((canonical, mate))))
    paired = len(paired_keys)
    require(len(zero_pcrs) == stable + 2 * paired, "ATLASP_ZERO_PCR_SPLIT")

    selected_codes = [code for code, value in enumerate(terminal) if value]
    astar_codes = [code for code, value in enumerate(astar) if value]
    seed_codes = [code for code, value in enumerate(seed) if value]
    geometry = residual_geometry(k, terminal)

    if skeleton_order is None:
        ordered_critical_sha = digest({"k": 8, "status": "K8_COMPLETE_PHASE_DIAGRAM_NOT_AUTHORIZED"})
        raw_critical = 0
        complete_range = False
    else:
        ordered_critical_sha = str(skeleton_order["ordered_critical_ratio_sha256"])
        raw_critical = int(skeleton_order["form_catalog"]["raw_critical_occurrence_count"])  # type: ignore[index]
        complete_range = True

    row = {
        "schema_version": ROW_SCHEMA_VERSION,
        "row_id": row_id,
        "construction_identity": {
            "q": 4,
            "k": k,
            "cell_kind": cell_kind,
            "cell_id": cell_id,
            "landmark_id": landmark_id,
            "evaluation_weight_identity_sha256": digest(weight),
            "alphabet_identity_sha256": digest({
                "alphabet": list(ALPHABET), "complement": list(COMPLEMENT),
                "order": list(ALPHABET), "encoding": "BASE4_BIG_ENDIAN",
            }),
            "selector_identity_sha256": digest({
                "seed": "NEGATIVE_RAY_STRICT_CROSSING_EQUIVARIANT_ZERO_PCR",
                "lower": "LAMBDA_POSITIVE_MIDPOINT_NEGATIVE_LEX_RC_TIE",
                "terminal": "SEED_PLUS_PREFIX_MINUS_SUFFIX",
            }),
            "decision_form_catalog_sha256": catalog_sha,
            "complete_covector_sha256": digest(covector),
        },
        "zero_and_descent": {
            "zero_embedding_word_count": len(zero_words),
            "zero_embedding_pcr_count": len(zero_pcrs),
            "zero_embedding_paired_pcr_count": paired,
            "stable_zero_pcr_count": stable,
            "lower_zero_state_count": lower_zero,
            "positive_lower_support_count": positive_lower,
            "exact_firing_count": firing_count,
            "a_star_count": len(astar_codes),
        },
        "terminal": {
            "minimum_cardinality": len(pcr_selected),
            "pcr_count": len(pcr_selected),
            "selected_word_count": len(selected_codes),
            "selected_word_sha256": digest(selected_codes),
            "a_star_count": len(astar_codes),
            "a_star_sha256": digest(astar_codes),
            "rc_mismatch_count": rc_mismatch,
            "non_boolean_membership_count": non_boolean,
            "pcr_underselection_count": under,
            "pcr_overselection_count": over,
            "direct_firing_mismatch_count": direct_mismatch,
        },
        "residual_dag": geometry,
        "comparisons": {
            "seed_terminal_symmetric_difference": len(set(seed_codes) ^ set(selected_codes)),
            "terminal_pairwise": [],
            "a_star_pairwise": [],
        },
        "phase_diagram": {
            "raw_form_occurrence_count": len(catalog),
            "structural_zero_occurrence_count": structural_count,
            "canonical_nonzero_form_count": canonical_count,
            "raw_critical_occurrence_count": raw_critical,
            "ordered_critical_ratio_sha256": ordered_critical_sha,
            "cell_adjacency_valid": True,
            "complete_declared_range": complete_range,
            "same_covector_different_payload_mismatch_count": 0,
        },
    }
    hidden = {
        "selected": set(selected_codes),
        "astar": set(astar_codes),
        "zero_words": zero_words,
        "covector_sha": row["construction_identity"]["complete_covector_sha256"],
        "geometry_sha": digest({
            "zero_and_descent": row["zero_and_descent"],
            "terminal": row["terminal"],
            "residual_dag": row["residual_dag"],
        }),
    }
    return row, hidden


def row_specs(skeleton_payload: dict[str, object]) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []
    for order in skeleton_payload["orders"]:  # type: ignore[index]
        require(isinstance(order, dict), "ATLASP_SKELETON_ORDER")
        k = int(order["k"])
        intervals = {item["cell_id"]: item for item in order["open_intervals"]}  # type: ignore[index]
        points = {item["cell_id"]: item for item in order["critical_points"]}  # type: ignore[index]
        for cell_id in order["ordered_cell_ids"]:  # type: ignore[index]
            if cell_id in intervals:
                weight = {"kind": "RATIONAL", "value": intervals[cell_id]["representative"]}
                kind = "OPEN_INTERVAL"
            else:
                require(cell_id in points, "ATLASP_SKELETON_CELL")
                weight = {"kind": "RATIONAL", "value": points[cell_id]["ratio"]}
                kind = "CRITICAL_POINT"
            specs.append({
                "k": k, "cell_id": cell_id, "cell_kind": kind,
                "landmark_id": "NONE", "weight": weight,
                "row_id": f"ATLAS-{cell_id}", "skeleton_order": order,
            })
    landmarks = skeleton_payload["landmarks"]
    require(isinstance(landmarks, list) and len(landmarks) == 4, "ATLASP_LANDMARKS")
    for index, item in enumerate(landmarks):
        require(isinstance(item, dict), "ATLASP_LANDMARK_ROW")
        specs.append({
            "k": 8, "cell_id": f"K8-L{index:03d}", "cell_kind": "LANDMARK",
            "landmark_id": item["landmark_id"], "weight": item["weight"],
            "row_id": f"ATLAS-K8-L{index:03d}-{item['landmark_id']}", "skeleton_order": None,
        })
    require(len(specs) == 30, "ATLASP_ROW_COUNT")
    return specs


def build_payload() -> dict[str, object]:
    started = time.monotonic()
    skeleton_receipt, skeleton_bytes = load_canonical(SKELETON_RECEIPT)
    require(skeleton_receipt.get("status") == "PASS_EXACT_CRITICAL_RATIO_SKELETON", "ATLASP_SKELETON_STATUS")
    skeleton_payload = skeleton_receipt["scientific_payload"]
    require(isinstance(skeleton_payload, dict), "ATLASP_SKELETON_PAYLOAD")
    require(skeleton_receipt["scientific_payload_sha256"] == digest(skeleton_payload), "ATLASP_SKELETON_PAYLOAD_HASH")
    rows: list[dict[str, object]] = []
    hidden_rows: dict[str, dict[str, object]] = {}
    for spec in row_specs(skeleton_payload):
        require(time.monotonic() - started < MAX_SECONDS, "ATLASP_WALL_CAP")
        row, hidden = selector_row(
            int(spec["k"]), spec["weight"], str(spec["row_id"]), str(spec["cell_id"]),
            str(spec["cell_kind"]), str(spec["landmark_id"]), spec["skeleton_order"],
        )
        rows.append(row)
        hidden_rows[str(row["row_id"])] = hidden

    by_k: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_k[int(row["construction_identity"]["k"])].append(row)  # type: ignore[index]
    for row in rows:
        row_id = str(row["row_id"])
        hidden = hidden_rows[row_id]
        k = int(row["construction_identity"]["k"])  # type: ignore[index]
        terminal_pairs = []
        astar_pairs = []
        for other in sorted(by_k[k], key=lambda item: str(item["row_id"])):
            other_id = str(other["row_id"])
            if other_id == row_id:
                continue
            other_hidden = hidden_rows[other_id]
            difference = hidden["selected"] ^ other_hidden["selected"]  # type: ignore[operator]
            zero_union = hidden["zero_words"] | other_hidden["zero_words"]  # type: ignore[operator]
            terminal_pairs.append({
                "other_row_id": other_id,
                "symmetric_difference": len(difference),
                "confined_to_zero_pcr_ties": difference <= zero_union,
            })
            astar_difference = hidden["astar"] ^ other_hidden["astar"]  # type: ignore[operator]
            astar_pairs.append({
                "other_row_id": other_id,
                "equal": len(astar_difference) == 0,
                "symmetric_difference": len(astar_difference),
            })
        row["comparisons"]["terminal_pairwise"] = terminal_pairs  # type: ignore[index]
        row["comparisons"]["a_star_pairwise"] = astar_pairs  # type: ignore[index]

    orders: list[dict[str, object]] = []
    skeleton_by_k = {int(order["k"]): order for order in skeleton_payload["orders"]}  # type: ignore[index]
    for k in (4, 6, 8):
        order_rows = by_k[k]
        groups: dict[str, list[str]] = defaultdict(list)
        covector_groups: dict[str, set[str]] = defaultdict(set)
        for row in order_rows:
            hidden = hidden_rows[str(row["row_id"])]
            groups[str(hidden["geometry_sha"])].append(str(row["row_id"]))
            covector_groups[str(hidden["covector_sha"])].add(str(hidden["geometry_sha"]))
        mismatches = sum(len(values) > 1 for values in covector_groups.values())
        require(mismatches == 0, f"ATLASP_COVECTOR_PAYLOAD_K{k}")
        links: list[dict[str, str]] = []
        if k in (4, 6):
            for link in skeleton_by_k[k]["landmark_locations"]:
                cell_id = str(link["cell_id"])
                links.append({
                    "landmark_id": str(link["landmark_id"]),
                    "cell_id": cell_id,
                    "row_id": f"ATLAS-{cell_id}",
                    "position_kind": str(link["position_kind"]),
                })
        else:
            for row in order_rows:
                links.append({
                    "landmark_id": str(row["construction_identity"]["landmark_id"]),  # type: ignore[index]
                    "cell_id": str(row["construction_identity"]["cell_id"]),  # type: ignore[index]
                    "row_id": str(row["row_id"]),
                    "position_kind": "LANDMARK",
                })
        orders.append({
            "k": k,
            "complete_declared_range": k in (4, 6),
            "ordered_row_ids": [str(row["row_id"]) for row in order_rows],
            "landmark_links": links,
            "terminal_geometry_groups": [
                {"terminal_geometry_identity_sha256": key, "row_ids": sorted(groups[key])}
                for key in sorted(groups)
            ],
            "unique_terminal_geometry_count": len(groups),
            "same_covector_different_payload_mismatch_count": mismatches,
        })

    weights = []
    for landmark in skeleton_payload["landmarks"]:  # type: ignore[index]
        weights.append({
            "landmark_id": landmark["landmark_id"],
            "weight": landmark["weight"],
            "evaluation_weight_identity_sha256": digest(landmark["weight"]),
        })
    payload = {
        "schema_version": SCHEMA_VERSION,
        "scope": SCOPE,
        "claim_boundary": CLAIM_BOUNDARY,
        "critical_skeleton_sha256": hashlib.sha256(skeleton_bytes).hexdigest(),
        "atlas_row_schema_sha256": file_sha256(ROW_SCHEMA_PATH),
        "assurance_scope": {
            "row_count": 30,
            "max_word_states_per_row": 65536,
            "complete_k4_k6": True,
            "complete_k8": False,
            "k8_landmark_count": 4,
        },
        "weights": weights,
        "orders": orders,
        "rows": rows,
    }
    require(len(canonical_bytes(payload)) < 1073741824, "ATLASP_OUTPUT_CAP")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        payload = build_payload()
        receipt = {
            "scientific_payload": payload,
            "scientific_payload_sha256": digest(payload),
            "status": STATUS,
        }
        data = canonical_bytes(receipt) + b"\n"
        if args.output is None:
            sys.stdout.buffer.write(data)
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_bytes(data)
    except (ProductionFailure, OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        code = str(exc) if isinstance(exc, ProductionFailure) else f"ATLASP_EXCEPTION_{type(exc).__name__}"
        print(json.dumps({"error_code": code, "status": "FAIL"}, sort_keys=True))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
