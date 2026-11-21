"""Compare an SO(7) BMN descendant with protected products at weight (6,3,3)."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations, product
import json
from math import prod
from pathlib import Path

from flint import fmpq_mat

from dual633_pairing import integral, load_source, outer_flavour, target_frame
from gaugino_primary_order import charge
import mixed_incoming as m
from so7_hair_first import lower, product as trace_product
from so8_level24_descendant import raising_images
from sp3_omega_product_replay import odd_jet, odd_slot, trace_factor
from sp3_omega_products import canonical, cartan
from su2_levi_images import native
from su3_levi_branching import sign, spectators
from trace_descent import dimensions, evaluate, settings


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "so7_cartan_target_20261118/representatives.json.gz"


def diagonal(polynomial):
    """Collect signed eigenvalue pairs in the same Weyl orbit coordinates as the product recipes."""
    out = Counter()
    for word, num in polynomial.items():
        if any(len(trace) % 2 for trace in word):
            continue
        jets = [tuple((a,b,sum(1 << j for j in theta)) for a,b,theta in tr) for tr in word]
        for planes in product(range(3), repeat=len(word)):
            sequence = tuple((i,l) for i,tr in zip(planes,jets) for l in tr)
            sign, ordered = canonical(sequence, lambda x: odd_jet(x[1]))
            if not sign:
                continue
            slots = tuple(tuple(l for j,l in ordered if j==i) for i in range(3))
            phase, key = canonical(slots, odd_slot)
            if phase:
                out[key] += sign * phase * 2**len(word) * num
    return {w:c for w,c in out.items() if c}


def export_trace(polynomial, output):
    nums, den = integral(polynomial)
    alphabet = sorted({l for w in nums for tr in w for l in tr})
    traces = sorted({tr for w in nums for tr in w})
    identifiers = {tr:i for i,tr in enumerate(traces)}
    with output.open("w") as stream:
        print(len(alphabet),len(traces),len(nums),file=stream)
        for l in alphabet:
            print(*native(l)[0],file=stream)
        for tr in traces:
            print(len(tr),*(alphabet.index(l) for l in tr),file=stream)
        for w,num in sorted(nums.items()):
            print(num*prod(native(l)[1] for tr in w for l in tr),len(w),*(identifiers[tr] for tr in w),file=stream)
    return len(nums), den, hashlib.sha256(output.read_bytes()).hexdigest()


def build(path, archive, output, residual_output):
    settings()
    source = load_source(path)
    source_q, _ = integral(m.q(source))
    frame = target_frame((3,3,1,1,1),9,771253)
    samples = evaluate([source_q],frame,7,"SO")[0]
    assert not any(samples)
    phi = 0, 0, (0,)
    quadratic = {((phi, phi),): 1}
    original = lower(source)
    descendant = trace_product(outer_flavour(original,0),quadratic)
    assert charge(descendant) == ((0,0,6,3,3),10)
    assert all(not any(row.values()) for row in raising_images(descendant))
    variation = m.q(descendant)
    propagated = trace_product(outer_flavour(lower(m.q(source)),0), quadratic)
    assert variation == {w:-v for w,v in propagated.items()}
    terms, den, digest = export_trace(descendant, output)
    restricted = diagonal(descendant)
    data = json.loads(gzip.decompress(archive.read_bytes()))
    assert (data["weight"],data["length"]) == ([0,0,6,3,3],10)
    columns = [tuple(tuple(map(tuple,s)) for s in col) for col in data["cartan_columns"]]
    recipes = [[(tuple(map(tuple,w)),c) for w,c in row] for row in data["highest_weight_recipes"]]
    product_rows = []
    for recipe in recipes:
        row = Counter()
        for factors, n in recipe:
            for w,c in cartan(factors).items():
                row[w] += n*c
        product_rows.append(row)
    square = fmpq_mat([[r.get(w,0) for w in columns] for r in product_rows])
    assert int(square.det()) % data["prime"] == data["minor"]["determinant"] != 0
    coordinates = square.transpose().solve(fmpq_mat([[str(restricted.get(w,0))] for w in columns]))
    coeff = [Fraction(str(coordinates[i,0])) for i in range(len(columns))]
    remainder = Counter(restricted)
    for row,n in zip(product_rows,coeff,strict=True):
        if n:
            for w,c in row.items():
                remainder[w] -= n*c
    residual = Counter(descendant)
    for row,mult in zip(recipes,coeff,strict=True):
        if mult:
            for factors,n in row:
                factor_product = {():1}
                for g in factors:
                    factor_product = trace_product(factor_product,trace_factor(g))
                for w,c in factor_product.items():
                    residual[w] -= mult*n*c
    residual = {w:c for w,c in residual.items() if c}
    assert diagonal(residual) == {w:c for w,c in remainder.items() if c}
    assert m.q(residual) == variation
    assert all(not any(row.values()) for row in raising_images(residual))
    corrected_terms, corrected_den, corrected_digest = export_trace(residual,residual_output)
    return dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), parent_closure_coordinates=len(samples),
                lowered_parent_terms=len(original), descendant_terms=terms, common_denominator=den,
                propagated_classical_variation_terms=len(variation), Cartan_orbit_coefficients=len(restricted),
                native_input_sha256=digest, product_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                product_coefficients=[[i,str(n)] for i,n in enumerate(coeff) if n],
                residual_Cartan_orbit_coefficients=sum(bool(c) for c in remainder.values()),
                residual_terms=corrected_terms, residual_denominator=corrected_den, residual_native_sha256=corrected_digest)


def minor(record, field, height, width, prime):
    mnr = record[field]
    rank = mnr["rank"]
    assert 0 < mnr["det"] < prime and 0 <= rank <= min(height,width)
    for indices, bound in ((mnr["rows"],height),(mnr["columns"],width)):
        assert len(set(indices))==len(indices)==rank and all(0<=i<bound for i in indices)
    return rank


def cartan_character(flavour=(6,3,3), length=10):
    """Three signed eigenvalue slots with graded permutation and flavour highest-weight projection."""
    shift = (2,1,0)
    weights = [(sign(p),tuple(x+d-shift[j] for x,d,j in zip(flavour,shift,p,strict=True)))
               for p in permutations(range(3))]
    cap = tuple(map(max,zip(*(w for _,w in weights))))
    terms = [(weight[:3]+(weight[-1],), (sum(weight[:3])+weight[-1])%2, n)
             for weight,n in spectators((*cap,0,0), length).items() if weight[-1]%2==0]

    @lru_cache(None)
    def trace(w, cycles):
        if not cycles:
            return int(not any(w))
        size = cycles[0]
        value = 0
        for degree,odd,n in terms:
            rest = tuple(x-size*y for x,y in zip(w,degree,strict=True))
            if min(rest)>=0:
                value += (-1)**((size-1)*odd)*n*trace(rest,cycles[1:])
        return value

    rows = []
    for phase,weight in weights:
        point = (*weight,length)
        numerator = trace(point,(1,1,1))+3*trace(point,(1,2))+2*trace(point,(3,))
        dimension, remainder = divmod(numerator,6)
        assert remainder==0
        rows.append(dict(weight=weight,sign=phase,dimension=dimension))
    return dict(weights=rows, highest_weight_dimension=sum(r["sign"]*r["dimension"] for r in rows))


def resolve(folder):
    result = {}
    for name, tag, group, colours, offset in (("so7","SO7","SO",7,99401),("sp3","Sp3","Sp",6,99787)):
        size = {y:sum(dimensions((0,0,6,3,3),y,colours,group).values()) for y in range(8,13)}
        ranks = {}
        for y in range(8,12):
            r = json.loads((folder/f"{name}_y{y}.json").read_bytes())
            assert (r["mode"],r["gauge"],r["sector"],r["length"],r["prime"],r["seed"],r["dimensions"]) == (
                "complete_classical_image",tag,[6,3,3,0,0],y,1000033,offset+7919*y,[size[y],size[y+1]])
            ranks[y] = minor(r,"minor",size[y+1],size[y],1000033)
        assert ranks[8]+ranks[9]==size[9] and ranks[10]+ranks[11]==size[11]
        result[tag] = dict(dimensions=size, exact_ranks=ranks, length_ten_cohomology=size[10]-ranks[9]-ranks[10])
    ranks = {}
    for label in ("candidate","residual"):
        r = json.loads((folder/f"{label}_so.json").read_bytes())
        assert (r["gauge"],r["mode"],r["sector"],r["length"],r["prime"],r["seed"]) == (
            "SO7","explicit_classical_image",[6,3,3,0,0],10,2147483647,317083)
        previous,centre,following = r["cochains"]
        assert [previous,centre,following]==[result["SO7"]["dimensions"][y] for y in (9,10,11)]
        incoming = minor(r,"incoming",centre,previous,r["prime"])
        outgoing = minor(r,"outgoing",following,centre,r["prime"])
        assert (incoming,outgoing)==(result["SO7"]["exact_ranks"][9],result["SO7"]["exact_ranks"][10])
        assert len(r["image_values_scaled"])==centre and all(0<=x<r["prime"] for x in r["image_values_scaled"])
        ranks[label] = minor(r,"augmented",centre,previous+1,r["prime"])
    assert ranks == dict(candidate=incoming+1,residual=incoming)
    descendant = json.loads((folder/"descendant.json").read_bytes())
    assert descendant["source_sha256"]==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert descendant["Cartan_orbit_coefficients"]>0 and descendant["residual_Cartan_orbit_coefficients"]==0
    cartan = json.loads((folder/"replay.json").read_bytes())
    assert cartan["cartan_rank"]==result["Sp3"]["length_ten_cohomology"]==result["SO7"]["length_ten_cohomology"]-1
    character = cartan_character()
    assert character["highest_weight_dimension"] == cartan["cartan_rank"]
    assert int(cartan["determinant"]) % 1000033 == cartan["determinant_residues"]["1000033"] != 0
    assert descendant["product_archive_sha256"]==cartan["recipe_sha256"]==hashlib.sha256((folder/"products.json.gz").read_bytes()).hexdigest()
    return dict(classical=result, SO7_incoming=incoming, SO7_augmented_mod_2147483647=ranks,
                common_products=cartan["cartan_rank"], full_Cartan_character=character,
                rational_Cartan_kernel={g:r["length_ten_cohomology"]-cartan["cartan_rank"] for g,r in result.items()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path)
    parser.add_argument("--residual-export", type=Path)
    parser.add_argument("--products", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if (not __debug__ or bool(args.export) == bool(args.folder)
            or bool(args.products)!=bool(args.export) or bool(args.residual_export)!=bool(args.export)):
        parser.error("use --export or --folder with enabled assertions")
    result = build(SOURCE,args.products,args.export,args.residual_export) if args.export else resolve(args.folder)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
