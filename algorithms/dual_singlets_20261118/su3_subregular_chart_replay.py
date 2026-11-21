"""Replay a scalar-chart first correction by integer polynomial multiplication."""
import argparse
from collections import Counter
import hashlib
import json
from math import lcm
from pathlib import Path

from su3_subregular_differential import qmonomial
from su3_subregular_image import ARCHIVES, ODD, write
from su3_subregular_root_coupled import d_variables, read_coefficients
from su3_subregular_root_linear import decode, root_degree
from su3_subregular_shards import digest


def assemble(folder, image):
    terms = {}
    for path in sorted(folder.glob("*.json.gz")):
        record = decode(path)
        spectator = "" if path.name == "unit.json.gz" else path.name.split(".")[0]
        if "solution" in record:
            if record["spectator"] != spectator or record["full_image_file_sha256"] != image:
                raise ValueError(f"Wrong component or image in {path}")
            record = record["solution"]
        if record["kind"] != "primitive":
            raise ValueError(f"Expected a primitive in {path}")
        part = read_coefficients(record["coefficients"])
        if any(d_variables(w).hex() != spectator for w in part) or terms.keys() & part.keys():
            raise ValueError(f"Overlapping or mislabeled components in {path}")
        terms.update(part)
    if not terms:
        raise ValueError("No correction coefficients found")
    return [[w.hex(), str(c)] for w, c in sorted(terms.items())]


def verify(image, recipe, roots, lowest_d, reference, *, also_through_roots=None, residual_path=None):
    if also_through_roots is not None and roots != 4:
        raise ValueError("Checking both filtrations requires the complete chart image")
    data = decode(image)
    level = data["level"]
    _, size, denominator, original_hash = ARCHIVES[level]
    sign = 1 if level == 32 else -1
    if (data["sign"], data["denominator"], data["full_image_words"], data["full_image_sha256"],
        data["degree"], data["root_limit"]) != (sign, denominator, size, original_hash, 15, 4):
        raise ValueError("Unexpected chart input or source convention")
    rows = data["coefficients"]
    leading = decode(reference / f"j{level}_through9.json.gz")
    if [r for r in rows if len(r[0]) <= 18] != leading["coefficients"]:
        raise ValueError("The original through-degree-nine expansion disagrees")
    if level == 34:
        receipt = json.loads((reference/"normal/a8_expansion.json").read_bytes())["receipt"]
        for degree, info in receipt["homogeneous"].items():
            part = [r for r in rows if len(r[0]) == 2*int(degree)]
            if len(part) != info["terms"] or digest(part) != info["payload_sha256"]:
                raise ValueError(f"A8 coefficient mismatch at degree {degree}")

    primitive = read_coefficients(recipe)
    scale = denominator
    for coefficient in primitive.values():
        scale = lcm(scale, coefficient.denominator)
    residual = {}
    degrees = Counter()
    field_degrees = Counter()
    previous = ""
    for label, n in rows:
        word = bytes.fromhex(label)
        if (label <= previous or not isinstance(n, int) or n == 0 or not 1 <= len(word) <= 15
                or sum(k >= 64 for k in word) != 4
                or list(word) != sorted(word) or any(k >= 128 or (k >= 64 and k % 16 == 0) for k in word)
                or any(i == j and ODD[i] for i,j in zip(word,word[1:]))):
            raise ValueError("Invalid colour monomial or coefficient")
        previous = label
        count = root_degree(word)
        degrees[len(word)] += 1
        field_degrees[count] += 1
        if count <= roots:
            residual[word] = sign*n*(scale//denominator)
    if dict(degrees) != {int(k):v for k,v in data["nonzero_by_degree"].items() if v}:
        raise ValueError("Polynomial degree counts disagree")
    selected = len(residual)
    for word, coefficient in primitive.items():
        if (list(word) != sorted(word) or any(k >= 128 or (k >= 64 and k%16 == 0) for k in word)
                or any(i == j and ODD[i] for i,j in zip(word,word[1:]))
                or sum(k >= 64 for k in word) != 3
                or root_degree(word) > roots or (len(d_variables(word)) < lowest_d and
                    (also_through_roots is None or root_degree(word) > also_through_roots))):
            raise ValueError("Primitive outside the selected filtration")
        integer = coefficient.numerator*(scale//coefficient.denominator)
        for target, n in qmonomial(word).items():
            if roots == 4 or root_degree(target) <= roots:
                value = residual.get(target, 0) + integer*n
                if value:
                    residual[target] = value
                else:
                    residual.pop(target, None)
    if any(len(d_variables(w)) >= lowest_d for w in residual):
        raise ValueError("The rational supercharge equation has a nonzero retained coefficient")
    result = dict(level=level, source_words=size, source_sha256=original_hash, sign_in_Q1=sign,
                full_chart_terms=len(rows), trace_groups=data["trace_groups"], degree_counts=dict(sorted(degrees.items())),
                fg_degree_counts=dict(sorted(field_degrees.items())), retained_fg_degree=roots,
                first_equation_input_terms=selected, primitive_terms=len(primitive), common_multiplier=scale,
                through_d_count=lowest_d, residual_terms_below_d_count=len(residual), retained_residual=0)
    if also_through_roots is not None:
        counts = Counter((len(d_variables(w)), root_degree(w)) for w in residual)
        if any(r <= also_through_roots for _, r in counts):
            raise ValueError("A coefficient survives at low off-diagonal ghost degree")
        result.update(also_through_fg_degree=also_through_roots,
                      residual_by_d_and_fg=[dict(d_degree=d, fg_degree=r, terms=n)
                                            for (d, r), n in sorted(counts.items())])
    if residual_path is not None:
        write(residual_path, dict(denominator=scale, coefficients=[[w.hex(), n] for w, n in sorted(residual.items())]))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image", type=Path, required=True)
    p.add_argument("--primitive", type=Path, required=True)
    p.add_argument("--components", type=Path, help="Assemble signed coefficients from a solve directory")
    p.add_argument("--roots", type=int, choices=range(5), required=True)
    p.add_argument("--through-d", type=int, choices=range(4), default=0)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("Enable assertions for the source replays")
    reference = Path(__file__).resolve().parent/"taylor_hooks/su3_scalar_levi/subregular"
    if args.components:
        recipe = assemble(args.components, hashlib.sha256(args.image.read_bytes()).hexdigest())
    else:
        recipe = decode(args.primitive)
    result = verify(args.image, recipe, args.roots, args.through_d, reference)
    if args.components:
        write(args.primitive, recipe)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))
