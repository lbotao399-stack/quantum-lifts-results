"""Evaluate SO(7) source recipes, every Q0 component, and the saved SF covector."""
import argparse
import gzip
import hashlib
import json
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from math import ceil, comb
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
P = 1000003
MU = (3, 3, 2, 2, 2)


def signature(content):
    from sp_full_profiles import profile
    return profile(content)


def field_profile(record):
    return tuple((tuple(pair), count) for pair, count in record["profile"])


def multi_profile(word):
    return signature(tuple(letter for trace in word for letter in trace))


def setup(runtime, native, plan, functional):
    global DATA
    sys.path[:0] = [native, runtime]
    import sp_environment
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents, sector
    import native_full_projection_sparse_alphabet_v1 as projector
    import native_so7_terminal_trace_v1 as evaluator
    import so7_q1_native_v1 as quantum
    plan = json.loads(Path(plan).read_text())
    src, tgt = plan["sources"], plan["targets"]
    profiles = {field_profile(r): r["index"] for r in tgt}
    def make_frame(row, length, seed):
        pattern = field_profile(row)
        choices = [w for w in contents(MU, length) if signature(w) == pattern]
        alphabet = sorted({l for w in choices for l in w})
        d = row["fermions"]
        k = d + 3
        return Frame(3, ceil((row["dimension"] + 24) / comb(k, d)), k, [d], alphabet,
                     p=P, seed=seed+row["index"], group="SO_odd")
    closure = {r["index"]: make_frame(r, 11, 81118300) for r in tgt}
    choices = {r["index"]: [w for w in contents(MU, 10) if signature(w) == field_profile(r)] for r in src}
    for r in src:
        if r["fermions"] == 4:
            assert all(len(l[2]) <= 1 for w in choices[r["index"]] for l in w)
    sf_frames = {}
    weight = (4, 4, 1, 1, 1)
    from fortuity.full.evalg import nferm
    for degree, k, points in ((7, 11, 5), (9, 12, 1)):
        alphabet = sorted({l for w in contents(weight, 10) if nferm((w,)) == degree for l in w})
        sf_frames[degree] = Frame(3, points, k, [degree], alphabet, p=P,
                                  seed=2026111643+900+degree, group="SO_odd")
    ell = np.load(functional, allow_pickle=False).astype(np.int64)
    assert hashlib.sha256(ell.astype("<u4").tobytes()).hexdigest() == "92cd8ac99187d25f2378dd724a3d62c5a6ff3936a0aee5226eb64f3348129318"
    assert [sf_frames[d].ncols for d in (7, 9)] == [1650, 220]
    assert ell.shape == (1870,) and np.all(ell[1650:] == 0)
    DATA = dict(src=src, profiles=profiles, closure=closure, choices=choices, source_frames={},
                make_frame=make_frame, sf=sf_frames, ell=ell, projector=projector, evaluator=evaluator,
                quantum=quantum, canonical=sp_environment.canonical, fw=sp_environment.fw,
                sector=sector(MU, 10), nferm=nferm, inputs={})


def candidate(index, rng):
    content = list(rng.choice(DATA["choices"][index]))
    rng.shuffle(content)
    traces = []
    while content:
        lengths = [i for i in range(2, len(content)+1) if len(content)-i != 1]
        size = rng.choice(lengths)
        traces.append(tuple(content[:size]))
        del content[:size]
    raw = DATA["canonical"]({tuple(traces): 1})
    return next(iter(raw), None)


def project(word):
    return DATA["projector"].project({word: 1}, DATA["sector"], P)


def images(expression, index):
    record = DATA["src"][index]
    classical = {}
    for w, c in DATA["canonical"](DATA["fw"].Q_expr(expression)).items():
        c %= P
        if c:
            target = DATA["profiles"].get(multi_profile(w))
            # The three target profiles exhaust this particular weight at length 11.
            assert target in record["targets"], (index, target, w)
            classical.setdefault(target, {})[w] = c
    quantum = {}
    if record["fermions"] >= 6:
        if index not in DATA["inputs"]:
            DATA["inputs"][index] = sorted({l for w in DATA["choices"][index] for l in w})
        quantum = DATA["quantum"].apply(expression, P, DATA["inputs"][index])
    return classical, quantum


def paired(expressions):
    # The final 220 covector components vanish identically in the published frame.
    part = [{w: c for w, c in e.items() if DATA["nferm"](w) == 7} for e in expressions]
    value = evaluated(part, DATA["sf"][7])
    return value @ DATA["ell"][:1650] % P


def evaluated(expressions, frame):
    # A projected trace may contain thousands of terms. Bound the temporary
    # word-by-coordinate matrix while preserving the full target frame.
    budget = max(32, (48 << 20) // (8 * frame.ncols))
    if sum(map(len, expressions)) <= budget:
        return DATA["evaluator"].evaluate(expressions, frame)
    rows = []
    for expression in expressions:
        terms = list(expression.items())
        result = np.zeros(frame.ncols, dtype=np.int64)
        for start in range(0, len(terms), budget):
            piece = dict(terms[start:start+budget])
            result = (result + DATA["evaluator"].evaluate([piece], frame)[0]) % P
        rows.append(result)
    return np.asarray(rows)


def evaluate(index, recipes, include_source):
    record = DATA["src"][index]
    expr = [project(w) for w in recipes]
    components = [images(e, index) for e in expr]
    values = {str(t): evaluated([a.get(t, {}) for a, b in components], DATA["closure"][t])
              for t in record["targets"]}
    values["ell_Q1"] = paired([b for a, b in components])
    if include_source:
        if index not in DATA["source_frames"]:
            DATA["source_frames"][index] = DATA["make_frame"](record, 10, 81118400)
        values["source"] = evaluated(expr, DATA["source_frames"][index])
    return {key: v.astype(np.uint32) for key, v in values.items()}


def sample(task):
    index, serial, count, output, include_source = task
    output = Path(output)
    prefix = output / ("p%02d_b%07d" % (index, serial))
    receipt = prefix.with_suffix(".json")
    if receipt.exists():
        existing = json.loads(receipt.read_text())
        assert existing["count"] == count and existing["include_source"] == include_source
        for suffix, key in ((".npz", "matrix_sha256"), (".json.gz", "recipes_sha256")):
            assert hashlib.sha256(prefix.with_suffix(suffix).read_bytes()).hexdigest() == existing[key]
        return existing
    rng = random.Random(81118500 + serial*104729 + index*8191)
    words, seen = [], set()
    while len(words) < count:
        word = candidate(index, rng)
        if word and word not in seen and project(word):
            seen.add(word)
            words.append(word)
    values = evaluate(index, words, include_source)
    path = prefix.with_suffix(".npz")
    tmp = prefix.with_suffix(".tmp.npz")
    np.savez(tmp, **values)
    tmp.replace(path)
    raw = gzip.compress(json.dumps(words, separators=(",", ":")).encode(), mtime=0)
    trace_path = prefix.with_suffix(".json.gz")
    trace_path.write_bytes(raw)
    record = dict(index=index, serial=serial, count=count, include_source=include_source,
                  matrix_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  recipes_sha256=hashlib.sha256(raw).hexdigest(),
                  nonzero_pairings=int(np.count_nonzero(values["ell_Q1"])),
                  frame_widths={key: list(v.shape) for key, v in values.items()})
    tmp = receipt.with_suffix(".tmp")
    tmp.write_text(json.dumps(record) + "\n")
    tmp.replace(receipt)
    return record


def control():
    from sp_full_projected_basis import evaluate as independent
    from q1d_sp_exact import Q1d_multi
    import full_extremal
    def residues(expression):
        out = {}
        for word, coefficient in expression.items():
            c = Fraction(coefficient)
            value = c.numerator * pow(c.denominator, -1, P) % P
            if value:
                out[word] = value
        return out
    reports = []
    for index in range(len(DATA["src"])):
        rng = random.Random(index + 71118)
        while True:
            word = candidate(index, rng)
            if word and (e := project(word)) and len(e) <= 36:
                break
        exact = full_extremal.project({word: 1}, DATA["sector"])
        assert residues(exact) == e
        a, b = images(e, index)
        total = {}
        for w, coefficient in exact.items():
            for term, value in Q1d_multi(w, 7, "SO_odd").items():
                total[term] = total.get(term, 0) + coefficient * value
        expected = residues(DATA["canonical"](total))
        if DATA["src"][index]["fermions"] >= 6:
            assert expected == b
        else:
            assert all(DATA["nferm"](w) != 7 for w in expected)
        assert not residues(DATA["canonical"](DATA["fw"].Q_expr(DATA["canonical"](DATA["fw"].Q_expr(exact)))))
        # Compare every closure block at one independent Python-evaluated point.
        import copy
        checked = 0
        for t in DATA["src"][index]["targets"]:
            frame = copy.copy(DATA["closure"][t])
            frame.P, frame.ncols = 1, frame.cpp
            frame.L = {l: (d, v[:1]) for l, (d, v) in frame.L.items()}
            part = [a.get(t, {})]
            assert np.array_equal(DATA["evaluator"].evaluate(part, frame), independent(part, frame))
            checked += frame.ncols
        pair = int(paired([b])[0])
        reports.append(dict(index=index, projected_terms=len(e), closure_coordinates=checked,
                            Q1_terms=len(expected), pairing=pair))
    archived = json.loads(gzip.decompress((HERE / "representatives.json.gz").read_bytes()))
    letters = [(a, b, tuple(theta)) for a, b, theta in archived["letters"]]
    sf = {tuple(tuple(letters[i] for i in trace) for trace in w): Fraction(n, d)
          for n, d, w in archived["SF"]}
    sf = residues(sf)
    values = np.concatenate([DATA["evaluator"].evaluate(
        [{w: c for w, c in sf.items() if DATA["nferm"](w) == degree}], DATA["sf"][degree])[0]
        for degree in (7, 9)])
    assert hashlib.sha256(values.astype("<u4").tobytes()).hexdigest() == "02935df21cba2bf673a8650a2659eaa37f0e3457eb4710c0a2757ccf470e0b39"
    assert int(values @ DATA["ell"] % P) == 1
    return dict(source_controls=reports, reconstructed_SF_coordinates=len(values), SF_pairing=1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runtime", required=True)
    p.add_argument("--native", required=True)
    p.add_argument("--plan", default=str(HERE / "source_profiles.json"))
    p.add_argument("--functional", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--control", action="store_true")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--batches", type=int, default=1)
    p.add_argument("--first", type=int, default=0)
    p.add_argument("--rows", type=int, default=4)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--include-source", action="store_true")
    p.add_argument("--dimension-cap", action="store_true")
    args = p.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    init = (str(Path(args.runtime).resolve()), str(Path(args.native).resolve()), args.plan, args.functional)
    if args.control:
        setup(*init)
        result = control()
        print(json.dumps(result))
        (args.output / "controls.json").write_text(json.dumps(result, indent=2) + "\n")
        return
    assert 0 <= args.shard < args.shards and args.workers >= 1
    tasks = []
    plan = json.loads(Path(args.plan).read_text())
    for serial in range(args.first, args.first+args.batches):
        for index in range(10):
            if args.dimension_cap and serial*args.rows >= (5*plan["sources"][index]["dimension"]+3)//4 + 32:
                continue
            if (11*serial+index) % args.shards == args.shard:
                tasks.append((index, serial, args.rows, str(args.output), args.include_source))
    with ProcessPoolExecutor(max_workers=args.workers, initializer=setup, initargs=init) as pool:
        for result in pool.map(sample, tasks):
            print(json.dumps({k: result[k] for k in ("index", "serial", "count", "nonzero_pairings")}), flush=True)


if __name__ == "__main__":
    main()
