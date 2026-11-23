"""Run unchanged repository calls with the opt-in exact enumerators.

Usage: python benchmark.py --baseline /path/to/pinned/checkout --out benchmark.json
The degree-nine baseline runs in its own bounded subprocess; degree ten is
only enumerated by the new method, without a claimed baseline timing.
"""
import argparse
from contextlib import contextmanager
import gc
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import platform
import resource
import subprocess
import sys
from time import perf_counter

import numpy as np

import dihedral_traces
import native_catalogue

if not __debug__:
    raise RuntimeError("Verification requires assertions; do not use -O or PYTHONOPTIMIZE")


PIN = "0b2e1a6724361fb2666089bfbe0c1a36a23d48a1"
RELATIVE = Path("algorithms/dual_singlets_20261118/f4/j36_source")


def load_baseline(root):
    source = root / RELATIVE
    sys.path.insert(0, str(source))
    sys.path.insert(0, str(source / "first_page"))
    import mixed_colour_relations
    import single_trace_flag
    return mixed_colour_relations, single_trace_flag


def digest_words(words):
    digest = sha256()
    for word in words:
        digest.update(repr(word).encode() + b"\n")
    return digest.hexdigest()


def digest_rows(rows):
    digest = sha256()
    coordinates = 0
    for indices, coefficients in rows:
        digest.update(len(indices).to_bytes(8, "little"))
        digest.update(indices.astype("<i4").tobytes())
        digest.update(coefficients.astype("<i4").tobytes())
        coordinates += len(indices)
    return digest.hexdigest(), coordinates


@contextmanager
def optimized_calls(trace_module, source_module):
    """Opt in for actual imported callers, restoring their APIs on exit."""
    old = trace_module.trace_words, source_module.catalogue, source_module.canonical
    trace_module.trace_words = dihedral_traces.trace_words
    source_module.catalogue = native_catalogue.catalogue
    source_module.canonical = native_catalogue.canonical14
    try:
        yield
    finally:
        trace_module.trace_words, source_module.catalogue, source_module.canonical = old


def worker(root, method, degree):
    trace, _ = load_baseline(root)
    trace.colour_trace.cache_clear()
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    call = trace.trace_words if method == "baseline" else dihedral_traces.trace_words
    started = perf_counter()
    words = call(degree)
    seconds = perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert len(words) == dihedral_traces.trace_count(degree)
    return dict(method=method, degree=degree, coordinates=len(words), seconds=seconds,
                pre_enumeration_peak_rss_kib=before, peak_rss_kib=peak,
                output_sha256=digest_words(words),
                original_permutation_visits=factorial(degree) * len(trace.partitions(degree, 2)))


def run_worker(root, method, degree):
    command = [sys.executable, str(Path(__file__).resolve()), "--baseline", str(root),
               "--worker", method, "--degree", str(degree)]
    process = subprocess.run(command, check=True, text=True, capture_output=True, timeout=45)
    return json.loads(process.stdout)


def run(root):
    sha = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if sha != PIN:
        raise ValueError(f"baseline commit must be {PIN}; got {sha}")
    trace, source = load_baseline(root)
    result = dict(baseline_commit=sha, python=sys.version, platform=platform.platform(),
                  timings_are_observations=True, scope="exact trace enumeration and source recipe replay",
                  compiler_command="c++ -O3 -std=c++17 -shared -fPIC dihedral_catalogue.cpp -o dihedral.so")
    labelled = []
    for degree in (6, 7, 8, 9):
        before = run_worker(root, "baseline", degree)
        after = run_worker(root, "optimized", degree)
        assert before["output_sha256"] == after["output_sha256"]
        labelled.append(dict(degree=degree, baseline=before, optimized=after,
                             speed_ratio=before["seconds"] / after["seconds"]))
        print("labelled", degree, before["seconds"], after["seconds"], flush=True)
    result["labelled_benchmarks"] = labelled
    result["degree_ten_new_only"] = run_worker(root, "optimized", 10)

    # Independently compare every signed adjacent-transposition coordinate.
    symmetry_checks = []
    for degree in (6, 7, 8):
        old = trace.trace_words(degree)
        new = dihedral_traces.trace_words(degree)
        assert old == new
        old_action = list(trace.signed_transpositions(old))
        new_action = list(trace.signed_transpositions(new))
        for (a, sa), (b, sb) in zip(old_action, new_action):
            assert a == b and np.array_equal(sa, sb)
        symmetry_checks.append(dict(degree=degree, checked_coordinates=len(old) * (degree - 1)))
    result["signed_action_checks"] = symmetry_checks

    # Actual native F4 26-dimensional colour matrices, with both published primes.
    representation = source.native26()
    evaluation_checks = []
    for degree in (6, 7, 8):
        before = trace.trace_words(degree)
        after = dihedral_traces.trace_words(degree)
        for prime in (1000003, 1000033):
            colours = trace.sample_colours(representation, degree, 0, prime)
            a = trace.evaluate_words(before, colours, prime)
            b = trace.evaluate_words(after, colours, prime)
            assert a == b
            evaluation_checks.append(dict(degree=degree, prime=prime, coordinates=len(a),
                                           sha256=sha256(np.array(a, dtype="<i8").tobytes()).hexdigest()))
    result["native_F4_evaluation_checks"] = evaluation_checks

    # This is a one-off startup cost, not counted in warm native timings.
    started = perf_counter()
    native_catalogue.library()
    result["native_library_setup_seconds"] = perf_counter() - started
    catalogues = []
    for repeat in range(3):
        source.catalogue.cache_clear()
        native_catalogue.catalogue.cache_clear()
        started = perf_counter()
        old = source.catalogue()
        before = perf_counter() - started
        started = perf_counter()
        new = native_catalogue.catalogue()
        after = perf_counter() - started
        assert old.dtype == new.dtype and np.array_equal(old, new)
        catalogues.append(dict(repetition=repeat, baseline_seconds=before,
                               optimized_seconds=after, speed_ratio=before/after))
    result["F4_degree_fourteen_catalogue"] = dict(input_content=[4,4,4,2],
        original_word_visits=3153150, exact_dihedral_coordinates=len(old),
        sha256=sha256(old.astype("<i4").tobytes()).hexdigest(), timings=catalogues)

    # Run all 474 archived source recipes through the original single(seed)
    # caller. Every sparse coordinate is retained and compared, with no slice.
    archive = json.loads((source.HERE / "single_trace_formal.json").read_text())
    seeds = archive["recipes"]
    assert len(seeds) == 474
    source.flags()
    rounds = []
    for repeat in range(3):
        gc.collect()
        started = perf_counter()
        old_rows = [source.single(seed) for seed in seeds]
        before = perf_counter() - started
        with optimized_calls(trace, source):
            started = perf_counter()
            new_rows = [source.single(seed) for seed in seeds]
            after = perf_counter() - started
        assert all(np.array_equal(a, b) and np.array_equal(c, d)
                   for (a,c), (b,d) in zip(old_rows, new_rows))
        old_digest, entries = digest_rows(old_rows)
        new_digest, _ = digest_rows(new_rows)
        assert old_digest == new_digest
        rounds.append(dict(repetition=repeat, baseline_seconds=before,
                           optimized_seconds=after, speed_ratio=before/after))
        del old_rows, new_rows
    result["all_archived_F4_source_rows"] = dict(recipes=len(seeds),
        exact_sparse_coordinates=entries, sha256=old_digest, timings=rounds,
        includes="source.single(seed) with catalogue and canonicalization; warm catalogue caches",
        excludes="later elimination, finite-colour rank jobs, and differential construction")

    # Repeated-content call sites used by the degree-eight F4 colour study.
    repeated = []
    words = trace.trace_words(8)
    for counts in ((5,3), (4,4), (2,2,2,2)):
        labels = tuple(i for i,n in enumerate(counts) for _ in range(n))
        started = perf_counter()
        old = sorted({word for w in words for word, sign in [trace.identify(w, labels)] if sign})
        before = perf_counter() - started
        started = perf_counter()
        new = dihedral_traces.trace_words_fixed_content(counts)
        after = perf_counter() - started
        assert old == new
        repeated.append(dict(content=counts, coordinates=len(old), baseline_specialization_seconds=before,
                              direct_seconds=after, speed_ratio=before/after,
                              excludes_baseline_labelled_enumeration=True))
    result["repeated_content_front_end"] = repeated
    files = [Path(__file__), Path(__file__).with_name("dihedral_traces.py"),
             Path(__file__).with_name("native_catalogue.py"),
             Path(__file__).with_name("dihedral_catalogue.cpp"),
             root / RELATIVE / "mixed_colour_relations.py",
             root / RELATIVE / "first_page/single_trace_flag.py",
             root / RELATIVE / "first_page/single_trace_formal.json"]
    result["source_sha256"] = {str(p.name): sha256(p.read_bytes()).hexdigest() for p in files}
    result["claims"] = dict(new_physical_lifts=0, new_fortuitous_generators=0,
        omitted_coordinates=0, changed_trace_orientations=0,
        degree_ten_baseline_not_run=True,
        no_claim_of_whole_solver_speedup=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--worker", choices=("baseline", "optimized"))
    parser.add_argument("--degree", type=int)
    arguments = parser.parse_args()
    if arguments.worker:
        print(json.dumps(worker(arguments.baseline.resolve(), arguments.worker, arguments.degree)))
    else:
        data = run(arguments.baseline.resolve())
        if arguments.out:
            arguments.out.write_text(json.dumps(data, indent=2) + "\n")
        else:
            print(json.dumps(data, indent=2))
