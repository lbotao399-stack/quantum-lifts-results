"""Verify all proper-product channels of the seven-letter SU(3) singlet."""
import argparse
from functools import cache
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ALG = HERE.parents[1] / "dual_singlets_20261118"
sys.path.insert(0, str(ALG))
from mixed33_exact import exact
from sp3_omega_products import all_products
from su3_bmn_parity_products import labelled
from trace_descent import dimensions


@cache
def dimension(length, flavour, plus, minus, odd):
    group = "SU-" if odd else "SU+"
    return sum(dimensions((plus, minus, flavour, 0, 0), length, 3, group).values())


@cache
def graviton_weight(length, flavour, plus, minus):
    if minus < 0:
        return 0
    recipes = all_products(3, (plus, minus, flavour, 0, 0), length, include_odd_powers=True)
    columns = [labelled(recipe, True) for family in recipes.values() for recipe in family]
    return len(exact(columns)[0])


def factor_pairs():
    # Tracelessness rules out length one. A spin singlet couples equal spins.
    for length in (2, 3):
        for flavour in range(4):
            for dots in range(9):
                for spin2 in range(dots % 2, min(dots, 8 - dots) + 1, 2):
                    left = (length, flavour, (dots + spin2) // 2, (dots - spin2) // 2)
                    right = (7-length, 3-flavour, (8-dots + spin2) // 2, (8-dots - spin2) // 2)
                    if all(f + p + m >= n and dimension(n, f, p, m, n % 2)
                           for n, f, p, m in (left, right)):
                        yield left, right


def verify(rebuild):
    archive = json.loads(gzip.decompress((HERE / "outgoing_factor_minors.json.gz").read_bytes()))
    images = {tuple(k): v for k, v in archive["matrices"]}
    pairs = list(factor_pairs())
    assert pairs == [(tuple(a), tuple(b)) for a, b in archive["pairs"]]
    factors = sorted({key for pair in pairs for key in pair})
    expected = {tuple(k): tuple(v) for k, v in archive["bounds"]}
    assert set(expected) == set(factors)
    required = {(n, f, p, m) for y, f, p, m in factors for n in (y-1, y)
                if all(any(dimension(l, f, p, m, parity) for parity in (0, 1))
                       for l in (n, n+1))}
    assert set(images) == required
    with TemporaryDirectory() as tmp:
        binary = Path(tmp) / "unitary_parity"
        if rebuild:
            src = ALG / "taylor_hooks"
            subprocess.run(["c++", "-O3", "-std=c++17", "-pthread", "-I"+str(src/"native"),
                            str(src/"unitary_parity.cpp"), "-o", str(binary)], check=True)
        for (length, flavour, plus, minus), data in images.items():
            sector = [flavour, 0, 0, plus, minus]
            assert data["mode"] == "complete_colour_image" and data["gauge"] == "SU3"
            assert data["sector"] == sector and data["length"] == length
            assert data["prime"] == 2147483647
            size = [sum(dimension(l, flavour, plus, minus, c) for c in (0, 1))
                    for l in (length, length+1)]
            assert size == [data["cochains"], data["target_evaluations"]]
            for sign in ("even", "odd"):
                minor = data[sign]
                assert 0 < minor["det"] < data["prime"]
                for key, bound in zip(("columns", "rows"), size):
                    indices = minor[key]
                    assert len(indices) == len(set(indices)) == minor["rank"]
                    assert all(0 <= i < bound for i in indices)
            if rebuild:
                args = [*sector, length, data["seed"], data["prime"], "complete-image"]
                raw = subprocess.check_output([str(binary), "SU3", *map(str, args)],
                                              env=dict(os.environ, PAR_THREADS="2"))
                assert json.loads(raw) == data

    for key in factors:
        length, flavour, plus, minus = key
        sign = "odd" if length % 2 else "even"
        adjacent = [(l, flavour, plus, minus) for l in (length-1, length)]
        upper = dimension(*key, length % 2) - sum(images[k][sign]["rank"] for k in adjacent if k in images)
        lower = graviton_weight(*key) - graviton_weight(length, flavour, plus+1, minus-1)
        assert (lower, upper) == expected[key] and lower == upper
    live = sum(all(expected[key][0] for key in pair) for pair in pairs)
    assert (len(images), len(factors), len(pairs), live) == (114, 112, 56, 26)
    return dict(classical_matrices=len(images), factor_sectors=len(factors),
                factor_pairs=len(pairs), pairs_with_both_factors_nonzero=live)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="regenerate every native differential minor")
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(verify(args.rebuild), indent=2))
