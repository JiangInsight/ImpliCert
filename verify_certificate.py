import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import re
import sys


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def binary64(bits):

    demand(isinstance(bits, str) and re.fullmatch(r"[0-9a-fA-F]{16}", bits),
           "expected 16 hex digits per weight")
    word = int(bits, 16)
    exponent = (word >> 52) & 2047
    significand = word & ((1 << 52) - 1)
    demand(word >> 63 == 0 and exponent != 2047 and word != 0,
           "weights must be strictly positive and finite")
    if exponent:
        significand += 1 << 52
    return significand * Fraction(2) ** (max(1, exponent) - 1023 - 52)


def verify(checkpoint, epsilon):
    demand(isinstance(checkpoint, dict) and
           checkpoint.get("schema") == "implicit-kclique-certificate-v1", "unknown checkpoint format")
    n, k = checkpoint["n"], checkpoint["k"]
    demand(type(n) is int and type(k) is int and 3 <= k <= n <= 12,
           "this exhaustive example verifier requires 3 <= k <= n <= 12")
    demand(isinstance(checkpoint["epsilon"], str) and Fraction(checkpoint["epsilon"]) == epsilon
           and epsilon > 0, "checkpoint tolerance differs from requested tolerance")
    vertices = range(n)

    edges = set()
    for e in checkpoint["edges"]:
        demand(isinstance(e, list) and len(e) == 2 and all(type(v) is int for v in e)
               and 0 <= e[0] < e[1] < n and tuple(e) not in edges, "invalid or duplicate edge")
        edges.add(tuple(e))
    selected = checkpoint["witness"]
    demand(isinstance(selected, list) and selected and
           all(type(v) is int and v in vertices for v in selected) and len(set(selected)) == len(selected),
           "invalid witness")
    witness = set(selected)
    demand(isinstance(checkpoint["weight_bits"], list) and len(checkpoint["weight_bits"]) == n,
           "expected one weight per vertex")
    weights = list(map(binary64, checkpoint["weight_bits"]))
    demand(isinstance(checkpoint["lower"], str) and isinstance(checkpoint["upper"], str),
           "bounds must be rational strings")
    lower, upper = Fraction(checkpoint["lower"]), Fraction(checkpoint["upper"])

    cliques = [c for c in combinations(vertices, k) if all(e in edges for e in combinations(c, 2))]
    count = sum(witness.issuperset(c) for c in cliques)
    demand(type(checkpoint["witness_count"]) is int and checkpoint["witness_count"] == count
           and lower == Fraction(count, len(witness)), "incorrect witness count or lower bound")
    loads = [Fraction(0) for _ in vertices]
    for c in cliques:
        denominator = sum(weights[v] for v in c)
        for v in c:
            loads[v] += weights[v] / denominator
    max_load = max(loads)
    demand(max_load <= upper, "claimed upper bound falls below the maximum exact load")
    demand(upper <= (1 + epsilon) * lower, "claimed bounds exceed requested approximation ratio")



    masks = [sum(1 << v for v in c) for c in cliques]
    optimum = max(Fraction(sum((m & s) == m for m in masks), s.bit_count()) for s in range(1, 1 << n))
    demand(lower <= optimum <= max_load, "bounds contradict exhaustive subset optimum")
    return len(cliques), lower, max_load, upper, optimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoints", nargs="+", type=Path)
    parser.add_argument("--epsilon", type=Fraction, default=Fraction(1, 100))
    args = parser.parse_args()
    try:
        for path in args.checkpoints:
            total, lower, max_load, upper, optimum = verify(json.loads(path.read_text()), args.epsilon)
            print(f"PASS {path.name}: {total} cliques; LB={lower}; max load={max_load}; "
                  f"claimed UB={upper}; exhaustive optimum={optimum}")
        return 0
    except (ValueError, KeyError, TypeError, ZeroDivisionError, OSError) as error:
        print(f"REJECT: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
