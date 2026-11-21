"""Evaluate all 293 first equations of the level-32 SU(3) source for a given correction."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from su3_level32_centralizer import general_colour_test
from su3_level32_layers import DATA, content, joint, physical, residual
from su3_level32_scalar_descent import restored


def check(data, change, engine, work):
    work.mkdir(exist_ok=True, parents=True)
    first = restored(joint.s.read(data/"scalar_layers/levi_correction.json.gz"))
    second = restored(joint.s.read(change))
    assert len(first) == 13400 and all(content(w).count(0) == 5 for w in first)
    assert len(second) == 238 and all(content(w).count(0) == 4 for w in second)
    first.update(second)
    rows = residual(data, first, lambda _: True)
    assert len(rows) == 293
    checks = {}
    by_power = defaultdict(Counter)
    for target, row in rows.items():
        power = target.count(0)
        result = {"formal_terms": len(row)}
        if not row:
            result["kind"] = "formal_zero"
        else:
            for seed in (894611, 620003):
                value, _ = general_colour_test(target, row, seed)
                if Fraction(value):
                    result.update(kind="nonzero", seed=seed, coefficient=value)
                    break
            else:
                if power >= 4:
                    test = physical.polynomials([row], engine, work)[0]
                    result.update(kind="nonzero" if test["nonzero_monomials"] else "polynomial_zero",
                                  certificate=test)
                else:
                    result["kind"] = "unresolved"
        checks[target.hex()] = result
        by_power[power][result["kind"]] += 1
    assert set(by_power) == {3,4,5,6,7}
    return dict(correction_terms=len(first), by_scalar_degree={str(k):dict(v) for k,v in sorted(by_power.items())},
                contents=checks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--delta", type=Path, required=True)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    result = check(a.data, a.delta, a.engine.resolve(), a.work)
    joint.s.write(a.work/"all_contents.json", result)
    print(result["by_scalar_degree"], flush=True)
