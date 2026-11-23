"""One exact classical row for both crossings and for additional sources in this diagram."""
import argparse
from fractions import Fraction
from functools import reduce
import gzip
from itertools import chain
import json
from math import gcd, lcm
from pathlib import Path
from zipfile import ZipFile

from flint import fmpz_mat, nmod_mat
from other_crossing import COLOURS, NORMALIZATION, POLES, ROOT
from solve_section import ANCHORS, NESTED

HERE = Path(__file__).resolve().parent
TARGET_SCALE = 28 * 252 * 84**2


def load(path):
    content = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    return json.loads(content)


def rational_classical_row(folder, saved, regenerate):
    selected = load(folder / "section-result.json")["selected_points"]
    rows = list(chain.from_iterable(load(f) for f in sorted(folder.glob("points-*.json.gz"))))
    assert [i for i, values in rows] == list(range(720))
    basis = load(folder / "selected-1009.json")["rows"]
    trace = load(folder / "trace_C4_ordered.json")["rows"]
    values = [sum(Fraction(a) * b for a, b in zip(row, ANCHORS, strict=True)) for row in trace]
    values[basis.index(["nested", 0])] = NESTED
    scales = {"reference": 1, "four": 28 * 252, "nested": 28 * 252 * 84}
    scores = [v * scales[tag] for v, (tag, _) in zip(values, basis, strict=True)]
    mat = [rows[i][1][:-1] for i in selected]
    assert len(mat) == len(mat[0]) == len(scores) == 703
    assert int(nmod_mat(mat, 1009).det()) == 134
    if regenerate:
        scale = lcm(*(v.denominator for v in scores))
        solution = fmpz_mat(mat).transpose().solve(fmpz_mat([[int(scale*v)] for v in scores]))
        ratios = [Fraction(int(v.numerator), int(v.denominator)) / scale for v in solution.entries()]
        den = lcm(*(v.denominator for v in ratios))
        coef = [int(v * den) for v in ratios]
        data = json.dumps({"denominator": hex(den), "numerators": list(map(hex, coef))},
                          separators=(",", ":")).encode()
        saved.write_bytes(gzip.compress(data, mtime=0))
    record = load(saved)
    den = int(record["denominator"], 16)
    coef = [int(i, 16) for i in record["numerators"]]
    assert len(coef) == len(mat) and den > 0 and reduce(gcd, coef, den) == 1
    for col, value in zip(zip(*mat, strict=True), scores, strict=True):
        assert sum(a*b for a,b in zip(coef,col, strict=True)) == den * value
    return selected, coef, den, rows


def replay(args):
    selected, coef, den, rows = rational_classical_row(args.folder, args.dual, args.regenerate)
    new = load(args.columns)
    assert len(new) == len(rows)
    assert all(i == j and old[-1] == b for (i,old),(j,a,b) in zip(rows,new,strict=True))
    section = [Fraction(sum(coef[j] * new[i][k] for j,i in enumerate(selected)), den * TARGET_SCALE)
               for k in (1,2)]
    assert section[1] == Fraction(-299015516814300,49)
    result = dict(classical_pivots_checked=len(coef), common_points=len(rows),
                  nonzero_A_evaluations=sum(bool(v[1]) for v in new),
                  rational_classical_A=str(section[0]), rational_classical_B=str(section[1]))
    if args.packet:
        with ZipFile(args.packet) as archive:
            for p in (1009,1013):
                name = ROOT + "/" + COLOURS + "/coordinates-%d.json" % p
                official = json.loads(archive.read(name))["moments"]
                mod = [v.numerator * pow(v.denominator, -1, p) % p for v in section]
                assert mod == [row[0] for row in official]
    if args.poles:
        inspected = load(args.poles)
        weighted = {"left": {}, "right": {}}
        for row in inspected:
            term = row["weight"] * sum(w * Fraction(s) for w,s in zip(POLES,row["ordered_poles"],strict=True))
            for side,i in row["entries"]:
                assert i not in weighted[side]
                weighted[side][i] = term
        assert sorted(weighted["left"]) == list(range(18))
        assert sorted(weighted["right"]) == list(range(30))
        values = {side: sum(entries.values()) / NORMALIZATION for side,entries in weighted.items()}
        assert values["left"] == values["right"]
        result.update(A_graph_entries={s:len(v) for s,v in weighted.items()},
                      A_ordered_graph_scores={s:str(v) for s,v in values.items()},
                      A_corrected_pairing=str(values["left"] - section[0]))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder", type=Path, default=HERE)
    p.add_argument("--dual", type=Path, default=HERE/"classical-score-row.json.gz")
    p.add_argument("--columns", type=Path, default=HERE/"reverse-columns.json")
    p.add_argument("--poles", type=Path)
    p.add_argument("--packet", type=Path)
    p.add_argument("--regenerate", action="store_true")
    print(json.dumps(replay(p.parse_args()), indent=2))
