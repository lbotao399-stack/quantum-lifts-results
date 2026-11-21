"""Extend the rational GL(2) wall primitive to homogeneous three-colour trace words."""
import argparse
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path

from su3_level32_centralizer import evaluate, restrict_words
from su3_level32_layers import DATA, content, formal, joint
from su3_levi_equation import allowed_targets, first_residual
from su3_level32_scalar_descent import restored, serialized
from verify import add_scaled


def x_power(exponent):
    if exponent < 0 or exponent == 1:
        raise ValueError("scalar power cannot be represented by quadratic and cubic traces")
    cubics = exponent % 2
    pairs = (exponent - 3*cubics)//2
    return (bytes(2),)*pairs + (bytes(3),)*cubics, Fraction((-1)**cubics, 6**(pairs+cubics))


def expanded_factors(centre, traces):
    # Always keep central Grassmann coefficients to the left of matrix traces.
    states = [(centre, (), Fraction(1))]
    for trace in traces:
        if len(trace) == 2:
            options = [(b"", (trace,), Fraction(1)), (trace, (), Fraction(-6))]
        elif len(trace) == 3:
            i, j, k = trace
            koszul = -1 if formal.ODD[j] and formal.ODD[k] else 1
            options = [(b"", (trace,), Fraction(1,2)),
                       (b"", (bytes((i,k,j)),), Fraction(-koszul,2))]
        else:
            raise ValueError("expected a pair trace or a triple trace")
        updated = []
        for old_centre, old_traces, old_scale in states:
            for new_centre, new_traces, scale in options:
                exponent = sum(formal.par(t) for t in old_traces)*formal.par(new_centre)
                updated.append((old_centre+new_centre, old_traces+new_traces,
                                old_scale * scale * (-1)**exponent))
        states = updated
    return states


def add_trace_product(vector, traces, coefficient):
    canonical = []
    for trace in traces:
        phase, term = formal.canon(trace)
        if not phase:
            return
        coefficient *= phase
        canonical.append(term)
    phase, key = formal.multi(tuple(canonical))
    if phase:
        add_scaled(vector, {key: phase}, coefficient)


def extension(degree, centre, traces):
    scalar, norm = x_power(degree)
    answer = {}
    for central_word, other, scale in expanded_factors(centre, traces):
        if not central_word:
            add_trace_product(answer, scalar + other, scale*norm)
            continue
        odd = [k for k in central_word if formal.ODD[k]]
        if len(set(odd)) != len(odd):
            continue
        if degree == 2:
            scale /= 3 * (-2)**len(central_word)
            add_trace_product(answer, (bytes(2)+central_word,) + other, scale)
            add_trace_product(answer, scalar + (central_word,) + other, -scale*norm)
        else:
            lower, lower_norm = x_power(degree-1)
            scale /= 3 * (-2)**len(central_word)
            add_trace_product(answer, scalar + (central_word,) + other, scale*norm)
            add_trace_product(answer, lower + (bytes((0,))+central_word,) + other, -scale*lower_norm)
    return answer


def assemble(level, certificates):
    vector = {}
    counts = {}
    for path in certificates:
        saved = joint.s.read(path)
        assert saved["level"] == level
        degree = saved["scalar_degree"]
        if degree in counts:
            raise ValueError("duplicate scalar degree")
        counts[degree] = 0
        for centre, block in saved["central_coefficients"].items():
            assert block["kind"] == "exact_solution"
            for traces, value in block["terms"]:
                lifted = extension(degree, bytes.fromhex(centre), joint.s.unpack_hex(traces))
                add_scaled(vector, lifted, Fraction(value)/saved["integer_rhs_denominator"])
                counts[degree] += 1
    assert set(counts) == {t.count(0) for t in allowed_targets(level)}
    return vector, counts


def check(level, data, vector):
    checks = {}
    for t, row in first_residual(level, data, vector, lambda _: True).items():
        denominator = lcm(*(c.denominator for c in row.values()))
        restricted = restrict_words(row, denominator)
        polynomial = evaluate(restricted)
        if polynomial:
            raise ValueError(f"nonzero residual for {t.hex()}: {min(polynomial.items())}")
        checks[t.hex()] = dict(formal_residual_terms=len(row), restricted_trace_terms=len(restricted),
                               nonzero_polynomial_coefficients=0)
    assert len(checks) == (293 if level == 32 else 500)
    return dict(level=level, correction_terms=len(vector), terms_by_scalar_degree=dict(sorted(Counter(
                content(w).count(0) for w in vector).items())), tested_contents=len(checks), contents=checks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", type=int, choices=(32,34), default=32)
    p.add_argument("--data", type=Path)
    p.add_argument("--blocks", nargs="+", type=Path)
    p.add_argument("--correction", required=True, type=Path)
    p.add_argument("--check", required=True, type=Path)
    a = p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    if a.blocks:
        vector, counts = assemble(a.level, a.blocks)
        joint.s.write(a.correction, serialized(vector))
        print("GL(2) primitive terms by degree", counts, flush=True)
    vector = restored(joint.s.read(a.correction))
    data = a.data or (DATA if a.level == 32 else Path(__file__).resolve().parents[1]/"su3_j34_20261118")
    report = check(a.level, data, vector)
    joint.s.write(a.check, report)
    print(report["correction_terms"], report["tested_contents"], flush=True)
