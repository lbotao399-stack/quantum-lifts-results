"""Recover the exterior multiplet of the SO(7) level-26 first supercharge image."""
import argparse
from collections import Counter
from functools import cache
import gzip
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import subprocess
import sys

from flint import fmpz_mat

HERE = Path(__file__).resolve().parent
ALG = HERE.parent
sys.path.insert(0, str(ALG))
import mixed_incoming as m
from dual633_pairing import load_source
from fortuity.bmn.words import w_gen
from fortuity.full.letters import Q_letter, is_ferm
from gaugino_primary_order import charge
from so7_hair_first import lower, product as trace_product, translated
from q1d_exact import pair_terms
from so7_level26.verify import calculate as verify_parent
from so8_level24_descendant import raising_images
from sp3_omega_product_replay import check_factor, product_coordinates
from sp3_omega_products import raising


def clean(row):
    return {k: v for k, v in row.items() if v}


@cache
def variation(letter, operator, direction):
    a, b, flavour = letter
    if operator == "E":
        old, new = direction
        if old in flavour and new not in flavour:
            sign = (-1)**sum(min(old, new) < k < max(old, new) for k in flavour)
            return sign, (a, b, tuple(sorted((set(flavour) - {old}) | {new})))
        return None
    sign = (-1)**sum(k > direction for k in flavour)
    if operator == "D" and direction not in flavour:
        return sign, (a, b, tuple(sorted((*flavour, direction))))
    degree = a + b + len(flavour)
    if operator == "S" and direction in flavour and degree > 1:
        return -(degree-1)*sign, (a, b, tuple(k for k in flavour if k != direction))
    return None


def pair_variation(row, operator, direction):
    result = Counter()
    for (x, y), coefficient in row.items():
        for slot, letter in enumerate((x, y)):
            term = variation(letter, operator, direction)
            if term:
                n, changed = term
                pair = (changed, y) if slot == 0 else (x, changed)
                result[pair] += coefficient*n*((-1)**is_ferm(x) if slot else 1)
    return result


def first(row):
    result = Counter()
    for pair, coefficient in row.items():
        for n, x, y in pair_terms(*pair):
            result[x, y] += coefficient*n
    return result


def check_operations():
    subsets = [s for k in range(4) for s in combinations(range(3), k)]
    letters = [(a,b,s) for a,b,s in product(range(2), range(2), subsets) if a or b or s]
    for letter, operator, direction in product(letters, ("D","S"), range(3)):
        tree = Counter({(x,y): n for n,x,y in Q_letter(letter)})
        residual = pair_variation(tree, operator, direction)
        term = variation(letter, operator, direction)
        if term:
            value, changed = term
            residual.update({(x,y): value*n for n,x,y in Q_letter(changed)})
        expected = Counter()
        if operator == "D":
            scalar = (0,0,(direction,))
            expected[letter,scalar] += 1
            expected[scalar,letter] -= 1
        assert clean(residual) == clean(expected)

    bmn = [l for l in letters if not any(l[:2])]
    for x,y,direction,operator in product(bmn, bmn, range(3), ("D","S")):
        row = {(x,y):1}
        residual = first(pair_variation(row, operator, direction))
        residual.update(pair_variation(first(row), operator, direction))
        assert not clean(residual)
    return dict(tree_checks=len(letters)*6, ordered_first_pair_checks_per_operator=len(bmn)**2*3)


def act(row, operator, direction):
    result = Counter()
    for word, coefficient in row.items():
        sign = 1
        for ti, trace in enumerate(word):
            for pos, letter in enumerate(trace):
                term = variation(letter, operator, direction)
                if term:
                    n, new = term
                    updated = trace[:pos]+(new,)+trace[pos+1:]
                    orientation, canonical = m.so_canon(word[:ti]+(updated,)+word[ti+1:])
                    if orientation:
                        result[canonical] += coefficient*sign*n*orientation
                if operator in ("D","S") and is_ferm(letter):
                    sign = -sign
    return clean(result)


def combination(*terms):
    result = Counter()
    for n, row in terms:
        result.update({w:n*k for w,k in row.items()})
    return clean(result)


def sources():
    data = ALG.parent/"so7_cartan_target_20261118/representatives.json.gz"
    assert hashlib.sha256(data.read_bytes()).hexdigest() == "eda182846ea9ac1f2ed989519f5b4a234e35ccf48f0ed6f120ce899b0ac37e20"
    primary = trace_product(lower(load_source(data)), translated(w_gen((0,),1)))
    assert len(primary) == 1668 and charge(primary) == ((0,0,5,4,4),10)
    lowered = act(primary, "E", (0,1))
    assert act(lowered,"E",(1,0)) == primary
    symmetric = act(primary, "D", 0)
    alternating = combination((1,act(primary,"D",1)),(-1,act(lowered,"D",0)))
    for row, w in ((symmetric,(6,4,4)),(alternating,(5,5,4))):
        assert charge(row) == ((0,0,*w),10)
        assert all(not clean(v) for v in raising_images(row))
    assert combination((1,act(symmetric,"S",0)),(8,primary),
                       (1,act(act(primary,"S",0),"D",0))) == {}
    assert combination((1,act(alternating,"S",1)),(10,primary),
                       (1,act(act(primary,"S",1),"D",1)),
                       (-1,act(act(lowered,"S",1),"D",0))) == {}
    return dict(parent_terms=len(primary), highest_weight_terms=[len(symmetric),len(alternating)],
                return_on_D_of_image=[-8,-10], return_on_Q1_of_descendants=[8,10])


def small_blocks():
    blocks = []
    for name, weight in (("444",(4,4,4)),("543",(5,4,3))):
        record = json.loads((HERE/f"native-{name}.json").read_bytes())
        assert record["sector"] == [*weight,0,0] and record["gauge"] == "SO7"
        assert (record["prime"],record["seed"]) == (1000003,772801)
        block, = record["blocks"]
        blocks.append([block[k] for k in ("dim","rank_in","rank_out","rho")])
    assert blocks == [[502,313,189,0],[2432,1546,873,13]]
    data = json.loads(gzip.decompress((HERE/"gravitons.json.gz").read_bytes()))
    assert data["group"] == "SO(7)" and data["weight"] == [0,0,5,4,3] and data["length"] == 10
    recipes = [[(tuple(map(tuple,w)),n) for w,n in row] for row in data["highest_weight_recipes"]]
    columns = {tuple(tuple(map(tuple,s)) for s in col):i for i,col in enumerate(data["cartan_columns"])}
    factors = {g for row in recipes for w,_ in row for g in w}
    for spec in factors:
        check_factor(spec)
    matrix = []
    for row in recipes:
        hw, coords = Counter(), Counter()
        for word,n in row:
            hw.update({v:n*k for v,k in raising(word).items()})
            coords.update({j:n*k for j,k in product_coordinates(word,columns).items()})
        assert not clean(hw)
        matrix.append([coords[j] for j in range(len(columns))])
    assert len(matrix) == len(columns) == 13
    determinant = int(fmpz_mat(matrix).det())
    assert determinant % 1000003 == 135801
    return dict(classical_blocks=blocks, product_determinant=str(determinant),
                product_residues=[determinant%p for p in (1000003,1000033)])


def exterior(axes):
    if len(axes) != len(set(axes)):
        return 0, ()
    return (-1)**sum(a>b for i,a in enumerate(axes) for b in axes[i+1:]), tuple(sorted(axes))


@cache
def lower_descendant(direction, axes, component):
    if not axes:
        return {}
    first_axis, rest = axes[0], axes[1:]
    result = Counter()
    if direction == first_axis:
        result[rest,component] -= 9+len(rest)
    if component == direction:
        result[rest,first_axis] += 1
    if direction in rest:
        sign, sorted_axes = exterior(tuple(first_axis if a==direction else a for a in rest))
        if sign:
            result[sorted_axes,component] += sign
    for (tail,j),n in lower_descendant(direction,rest,component).items():
        sign, sorted_axes = exterior((first_axis,)+tail)
        if sign:
            result[sorted_axes,j] -= n*sign
    return clean(result)


def returns():
    records = []
    for degree in range(4):
        states = list(product(combinations(range(3),degree),range(3)))
        matrix = []
        for tests, component in states:
            row = []
            for state in states:
                values = {state:1}
                for direction in tests:
                    following = Counter()
                    for (axes,i),n in values.items():
                        following.update({k:n*v for k,v in lower_descendant(direction,axes,i).items()})
                    values = clean(following)
                row.append(values.get(((),component),0))
            matrix.append(row)
        mat = fmpz_mat(matrix)
        assert mat.rank() == len(states)
        if degree == 3:
            assert all(mat[i,j] == (-880 if i==j else 0) for i,j in product(range(3),repeat=2))
        records.append(dict(degree=degree, dimension=len(states), determinant=str(mat.det()),
                            characteristic_coefficients=[int(v) for v in mat.charpoly()]))
    return records


def public_census(engine):
    raw = gzip.decompress((HERE/"public-j28.jsonl.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == "85772f819fe6aeac061271772dc0433ebd1d54f1f5316a180d17242c67acc099"
    rows = {}
    for line in raw.splitlines():
        record = json.loads(line)
        weight = tuple(record["lam"])
        key = weight,14-record["n"]
        assert sum(weight)==14 and record["group"]=="SO_odd" and record["N"]==3
        assert record["check"] and key not in rows
        rows[key] = record
    assert len(rows)==147
    dims = {}
    for a in range(15):
        for b in range(a+1):
            c = 14-a-b
            if not 0<=c<=b:
                continue
            command = [str(engine),"dims","SO7",*map(str,(a,b,c)), "0","0","1","14"]
            block = json.loads(subprocess.check_output(command))
            for y,dim in block["d_Y"].items():
                dims[(a,b,c),int(y)] = dim
    assert dims == {k:v["m"] for k,v in rows.items() if v["m"]}
    nonzero = []
    for (weight,y),dim in sorted(dims.items()):
        v = rows[weight,y]
        incoming_and_products = v["r_up"]+v["grav"]
        outgoing = rows.get((weight,y+1),{}).get("r_up",0)
        defect = dim-incoming_and_products-outgoing
        assert defect>=0
        if defect:
            nonzero.append([list(weight), y, dim, v["r_up"], v["grav"], outgoing, defect])
    assert nonzero == [[[5,5,4],10,10217,4057,15,6144,1],
                       [[6,4,4],10,13755,5412,37,8305,1]]
    return dict(nonempty_blocks=len(dims), defects=nonzero)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine",type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable assertions")
    parent = verify_parent(ALG/"so7_level26")
    result = dict(parent_rank_pair=[parent["classical_rank_over_Q"],parent["augmented_rank_over_Q"]],
                  operators=check_operations(), traces=sources(), lower_targets=small_blocks(),
                  returns=returns(), full_bmn_level28=public_census(args.engine.resolve()))
    print(json.dumps(result,indent=2))
