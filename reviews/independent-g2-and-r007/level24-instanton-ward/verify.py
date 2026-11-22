"""Check the lowered G2 source and the return on its proposed instanton target."""
import argparse
import contextlib
import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from return_product import run as check_return


HERE = Path(__file__).resolve().parent
PRIMES = (2147483647, 1000033)


def lowered_weights():
    # Sym^3(3) tensor exterior^2(3) has components [3,1] and [2,0].
    for t in range(2, 5):
        for weight in ((7-t, 4-t, 3-t), (6-t, 4-t, 4-t)):
            if min(weight) >= 0:
                yield (*weight, t, t-2), 10-t


def stored_chain(key, prime):
    name = f"chain-{key}" + ("" if prime == PRIMES[0] else f"-{prime}")
    return json.loads((HERE / (name + ".json")).read_bytes())


def run_json(program, arguments):
    return json.loads(subprocess.check_output([str(program), *map(str, arguments)]))


def verify(repo, endpoint, charged, workers, regenerate, native):
    result = []
    sys.path.insert(0, str(charged))
    if regenerate:
        from build_g2_full_chain import build
        from g2_product_words import run as product_words
        from rebuild import certify
        with contextlib.redirect_stdout(sys.stderr):
            build(repo, native)
        executable = native / "span"
        subprocess.run(["c++", "-std=c++17", "-O3", "-pthread", "-I"+str(native),
                        str(charged / "g2_product_span.cpp"), "-o", str(executable)], check=True)

    for weight, length in lowered_weights():
        key = "".join(map(str, weight))
        recipes = HERE / f"prod-{key}-all.json.gz"
        if recipes.exists():
            packed = json.loads(gzip.decompress(recipes.read_bytes()))
            assert packed["weight"] == list(weight[3:] + weight[:3]) and packed["length"] == length
            table = native / (key + ".dat")
            if regenerate:
                with contextlib.redirect_stdout(sys.stderr):
                    product_words(argparse.Namespace(repo=repo, recipes=recipes, out=table))

        baseline = None
        proofs = []
        for prime in PRIMES:
            tree = stored_chain(key, prime)
            assert tree["prime"] == prime and tree["weight"] == list(weight)
            row, = tree["blocks"]
            assert row["length"] == length
            if regenerate:
                rebuilt = run_json(native / "g2_dotted_chain",
                    (*weight[:3], length, length, tree["seed"], workers, prime, *weight[3:]))
                assert rebuilt == tree
            incoming, outgoing, dimension = (row[k] for k in ("incoming", "outgoing", "dimension"))
            joined = incoming
            selected = []
            if recipes.exists():
                span = json.loads((HERE / f"joined-{key}-{prime}.json").read_bytes())
                assert span["weight"] == packed["weight"] and span["prime"] == prime
                assert span["length"] == length and span["products"] == len(packed["recipes"])
                assert span["dimension"] == dimension and span["boundary"] == incoming
                selected = span["selected"]
                joined = span["augmented"]
                assert len(set(selected)) == joined-incoming == len(selected)
                if regenerate:
                    proof = certify(executable, table, weight, length, prime, span["seed"], workers, native)
                    assert all(proof[k] == v for k,v in span.items())
                    proofs.append(dict(prime=prime, boundary_determinant=proof["boundary_minor"]["det"],
                                       joined_determinant=proof["joined_minor"]["det"],
                                       joined_minor_sha256=proof["joined_minor"]["sha256"]))
            assert joined + outgoing == dimension
            ranks = (dimension, incoming, outgoing, joined, len(selected))
            assert baseline is None or baseline == ranks
            baseline = ranks
        result.append(dict(weight=list(weight), length=length, dimension=baseline[0],
                           incoming=baseline[1], outgoing=baseline[2], joined=baseline[3],
                           products_added=baseline[4], minors=proofs))
    target = check_return(repo, endpoint)
    assert target == json.loads((HERE / "return.json").read_bytes())
    assert len(result) == 5 and target["return_coefficient"] == 6
    return dict(source_blocks=result, return_coefficient=6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--endpoint", required=True, type=Path)
    parser.add_argument("--charged", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("enable assertions and select a positive worker count")
    with TemporaryDirectory() as scratch:
        value = verify(args.repo.resolve(), args.endpoint.resolve(), args.charged.resolve(),
                       args.workers, args.rebuild, Path(scratch))
    print(json.dumps(value, indent=2))
