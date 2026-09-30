#!/usr/bin/env python3
"""Independent finite checks for the single-odd-prime embedding formula.

The manuscript supplies the general proof.  This standard-library checker
tests the cyclotomic composition identity, exact binary word counts, one
independent orbit count, and a two-odd-prime negative control.  It does not
claim to prove the all-order theorem by enumeration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path


SCHEMA = "FCIG_CWC_SINGLE_ODD_PRIME_EMBEDDING_V1"
STATUS = "PASS_SINGLE_ODD_PRIME_EMBEDDING_FINITE_ASSURANCE"


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, code: str) -> None:
    if not condition:
        raise CheckFailure(code)


def canonical(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def trim(values: list[int]) -> list[int]:
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def divisors(n: int) -> list[int]:
    return [value for value in range(1, n + 1) if n % value == 0]


def totient(n: int) -> int:
    return sum(math.gcd(value, n) == 1 for value in range(1, n + 1))


def polynomial_division_exact(numerator: list[int], denominator: list[int]) -> list[int]:
    work = trim(list(numerator))
    divisor = trim(list(denominator))
    require(divisor[-1] == 1, "CYCLOTOMIC_DIVISOR_MONIC")
    if len(work) < len(divisor):
        raise CheckFailure("CYCLOTOMIC_DIVISION_DEGREE")
    quotient = [0] * (len(work) - len(divisor) + 1)
    while len(work) >= len(divisor):
        coefficient = work[-1]
        offset = len(work) - len(divisor)
        quotient[offset] = coefficient
        for index, value in enumerate(divisor):
            work[offset + index] -= coefficient * value
        trim(work)
    require(work == [0], "CYCLOTOMIC_DIVISION_REMAINDER")
    return trim(quotient)


@lru_cache(maxsize=None)
def cyclotomic(n: int) -> tuple[int, ...]:
    polynomial = [-1] + [0] * (n - 1) + [1]
    for divisor in divisors(n)[:-1]:
        polynomial = polynomial_division_exact(polynomial, list(cyclotomic(divisor)))
    return tuple(polynomial)


def compose_power(polynomial: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    output = [0] * ((len(polynomial) - 1) * exponent + 1)
    for index, value in enumerate(polynomial):
        output[index * exponent] = value
    return tuple(trim(output))


def polynomial_remainder(numerator: list[int], denominator: tuple[int, ...]) -> tuple[int, ...]:
    work = trim(list(numerator))
    divisor = list(denominator)
    degree = len(divisor) - 1
    while len(work) > degree:
        coefficient = work[-1]
        offset = len(work) - len(divisor)
        for index, value in enumerate(divisor):
            work[offset + index] -= coefficient * value
        trim(work)
    work.extend([0] * (degree - len(work)))
    return tuple(work)


def monomial_remainders(k: int) -> list[tuple[int, ...]]:
    phi = cyclotomic(k)
    return [polynomial_remainder([0] * index + [1], phi) for index in range(k)]


def half_sums(vectors: list[tuple[int, ...]]) -> Counter[tuple[int, ...]]:
    degree = len(vectors[0])
    base = [-sum(vector[column] for vector in vectors) for column in range(degree)]
    output: Counter[tuple[int, ...]] = Counter()
    for mask in range(1 << len(vectors)):
        value = list(base)
        for index, vector in enumerate(vectors):
            if mask & (1 << index):
                for column in range(degree):
                    value[column] += 2 * vector[column]
        output[tuple(value)] += 1
    return output


def binary_zero_count(k: int) -> int:
    vectors = monomial_remainders(k)
    split = k // 2
    left = half_sums(vectors[:split])
    right = half_sums(vectors[split:])
    return sum(
        multiplicity * left.get(tuple(-entry for entry in value), 0)
        for value, multiplicity in right.items()
    )


def binary_zero_words(k: int) -> list[tuple[int, ...]]:
    vectors = monomial_remainders(k)
    degree = len(vectors[0])
    output: list[tuple[int, ...]] = []
    for mask in range(1 << k):
        word = tuple(1 if mask & (1 << index) else -1 for index in range(k))
        total = [0] * degree
        for sign, vector in zip(word, vectors):
            for column in range(degree):
                total[column] += sign * vector[column]
        if not any(total):
            output.append(word)
    return output


def canonical_rotation(word: tuple[int, ...]) -> tuple[int, ...]:
    return min(word[index:] + word[:index] for index in range(len(word)))


def formula_word_count(q: int, p: int, a: int, b: int) -> int:
    block = q**p + q + (2 ** (p - 1)) * q * (q - 2)
    exponent = (2 ** (b - 1)) * (p ** (a - 1))
    return block**exponent


def formula_orbit_count(q: int, k: int, zero_words: int) -> int:
    fixed = sum(totient(k // d) * q**d for d in divisors(k) if d < k)
    require((zero_words + fixed) % k == 0, "ORBIT_FORMULA_INTEGRALITY")
    return (zero_words + fixed) // k


def build_payload() -> dict[str, object]:
    identities: list[dict[str, int]] = []
    for p in (3, 5, 7):
        for a in (1, 2):
            for b in (1, 2, 3):
                k = (2**b) * (p**a)
                exponent = (2 ** (b - 1)) * (p ** (a - 1))
                require(
                    cyclotomic(k) == compose_power(cyclotomic(2 * p), exponent),
                    f"CYCLOTOMIC_COMPOSITION_{p}_{a}_{b}",
                )
                identities.append({"a": a, "b": b, "k": k, "p": p, "r": exponent})

    positive: list[dict[str, int]] = []
    for p, a, b in ((3, 1, 1), (5, 1, 1), (3, 1, 2), (3, 2, 1), (5, 1, 2), (3, 1, 3)):
        k = (2**b) * (p**a)
        observed = binary_zero_count(k)
        expected = formula_word_count(2, p, a, b)
        require(observed == expected, f"BINARY_WORD_COUNT_{k}")
        positive.append({
            "a": a, "b": b, "expected": expected, "k": k,
            "observed": observed, "p": p,
        })

    words_12 = binary_zero_words(12)
    observed_orbits_12 = len({canonical_rotation(word) for word in words_12})
    expected_orbits_12 = formula_orbit_count(2, 12, len(words_12))
    require(observed_orbits_12 == expected_orbits_12 == 19, "K12_ORBIT_COUNT")

    observed_30 = binary_zero_count(30)
    require(observed_30 == 146854, "K30_EXACT_NEGATIVE_CONTROL")
    mutated_p3 = formula_word_count(2, 3, 1, 1) ** 5
    mutated_p5 = formula_word_count(2, 5, 1, 1) ** 3
    require(observed_30 not in {mutated_p3, mutated_p5}, "K30_DOMAIN_GUARD")
    pcr_30 = formula_orbit_count(2, 30, observed_30)
    require(pcr_30 == 6069, "K30_ORBIT_CONTROL")

    return {
        "claim_boundary": (
            "FINITE_EXACT_ASSURANCE_FOR_THE_SINGLE_ODD_PRIME_FORMULA;"
            "GENERAL_AUTHORITY_REMAINS_THE_MANUSCRIPT_PROOF"
        ),
        "cyclotomic_composition_cases": identities,
        "negative_controls": [
            {
                "control": "missing_two_power_block_factor_at_k12",
                "mutated": formula_word_count(2, 3, 1, 1),
                "observed": len(words_12),
                "rejected": formula_word_count(2, 3, 1, 1) != len(words_12),
            },
            {
                "control": "two_odd_primes_miscast_as_p3_at_k30",
                "mutated": mutated_p3, "observed": observed_30, "rejected": True,
            },
            {
                "control": "two_odd_primes_miscast_as_p5_at_k30",
                "mutated": mutated_p5, "observed": observed_30, "rejected": True,
            },
        ],
        "orbit_controls": [
            {"k": 12, "observed_pcrs": observed_orbits_12, "predicted_pcrs": expected_orbits_12},
            {"k": 30, "observed_zero_words": observed_30, "predicted_pcrs": pcr_30},
        ],
        "positive_binary_word_cases": positive,
    }


def document() -> dict[str, object]:
    payload = build_payload()
    return {
        "schema_version": SCHEMA,
        "scientific_payload": payload,
        "scientific_payload_sha256": hashlib.sha256(canonical(payload)).hexdigest(),
        "status": STATUS,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    value = document()
    encoded = canonical(value) + b"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(encoded)
    print(json.dumps({
        "composition_cases": len(value["scientific_payload"]["cyclotomic_composition_cases"]),
        "payload_sha256": value["scientific_payload_sha256"],
        "status": value["status"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
