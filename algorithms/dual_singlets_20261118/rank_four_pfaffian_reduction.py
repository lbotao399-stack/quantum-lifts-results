"""Reduce terminal USp(8) traces to four factors by polarizing the rank-eight Pfaffian."""
import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

from flint import fmpz_mat

from rank_four_terminal import bounds, quotient, relation
from terminal_commutators import audit_rules, words


def layout(word):
    species, ends, fixed = [], [], []
    for j, tr in enumerate(word):
        first = len(species)
        species.extend(i for i, n in enumerate(tr) for _ in range(n))
        last = len(species) - 1
        assert first < last
        for k in range(first, last):
            fixed.extend((2*k+1, 2*k+2))
        if j < 5:
            ends.extend((2*first, 2*last+1))
        else:
            fixed.extend((2*last+1, 2*first))
    assert len(ends) == 10
    return species, ends+fixed


def run(output):
    checks = audit_rules()
    sectors, witness = [], []
    for t in range(6):
        case = bounds(t)
        case["classical"] = case["cartan"]
        sectors.append(case)
        if not case["formal"]:
            continue
        weight = case["weight"]
        basis = words(weight)
        indices = {w: i for i, w in enumerate(basis)}
        excess = [w for w in basis if len(w) > 4]
        digest = hashlib.sha256()
        stream = Counter()
        projected = []
        if t == 3:
            singlet_basis, cov = quotient(weight, case["formal"])
            assert list(basis) == list(singlet_basis)
        for w in excess:
            species, slots = layout(w)
            row = relation(species, slots)
            leading = row[w]
            assert abs(leading) == 1
            assert all(len(v) < len(w) for v in row if v != w)
            terms = sorted((indices[v], c) for v, c in row.items())
            stream[len(w)] += 1
            digest.update((json.dumps(terms, separators=(",", ":"))+"\n").encode())
            if t == 3 and len(projected) < case["formal"]-case["cartan"]:
                values = [sum(c*int(cov[i, j]) for i, c in terms) for j in range(cov.ncols())]
                if fmpz_mat(projected+[values]).rank() > len(projected):
                    projected.append(values)
                    witness.append(dict(word=indices[w], species=species, slots=slots,
                                        raw=terms, singlet=values))
        case.update(weight_words=len(basis), triangular_rows=len(excess),
                    leading_rows_by_trace_count=dict(sorted(stream.items())),
                    triangular_row_sha256=digest.hexdigest())
        if t == 3:
            assert len(projected) == 7
            for cols in combinations(range(case["formal"]), len(projected)):
                det = int(fmpz_mat([[v[j] for j in cols] for v in projected]).det())
                if det:
                    case["singlet_minor"] = dict(columns=cols, determinant=det)
                    break
            assert "singlet_minor" in case
        print(t, len(basis), len(excess), case["classical"], flush=True)
    output.write_text(json.dumps(dict(group="USp(8)", level=30, field="Q",
                                      pair_rules=checks, sectors=sectors, mixed_relations=witness),
                                 separators=(",", ":"))+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    assert __debug__
    run(p.parse_args().out)
