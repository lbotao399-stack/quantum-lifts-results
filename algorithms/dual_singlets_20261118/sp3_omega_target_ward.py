"""Exact Cartan and superconformal tests on the physical successors of omega."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations
import json
from pathlib import Path

from flint import fmpz_mat
import mixed_incoming as m
from so8_level24_descendant import substitute, substitution
from sp3_omega_actions import centre_basis, load_page
from sp3_omega_products import canonical, cartan, letter_odd, slot_odd
from sp3_omega_replay import pair
from terminal_commutators import audit_rules, sign, words
from trace_descent import dimensions


def clean(row):
    return {key:c for key,c in row.items() if c}


@lru_cache(maxsize=50000)
def cartan_substitution(slots, spin, flavour, increase):
    row = Counter()
    phase = 1
    for i,slot in enumerate(slots):
        for j,jet in enumerate(slot):
            x,y,mask = jet
            theta = tuple(k for k in range(3) if mask & (1<<k))
            scale, moved = substitution((x,y,theta), spin, flavour, increase)
            if scale:
                u,v,new_theta = moved
                new = (u,v,sum(1<<k for k in new_theta))
                inside, part = canonical(slot[:j]+(new,)+slot[j+1:],letter_odd)
                outside, key = canonical(slots[:i]+(part,)+slots[i+1:],slot_odd)
                row[key] += scale*phase*inside*outside
            if letter_odd(jet):
                phase = -phase
    return tuple(clean(row).items())


def act(row,spin,flavour,increase):
    result = Counter()
    for slots,n in row.items():
        for key,k in cartan_substitution(slots,spin,flavour,increase):
            result[key] += n*k
    return clean(result)


def returned(row, delta):
    result = Counter()
    for spin in (0,1):
        result.update(delta(delta(row,spin,0,True),spin,0,False))
    return clean(result)


def difference(a, b, scale):
    value = Counter(a)
    value.subtract({w:scale*c for w,c in b.items()})
    return clean(value)


def cartan_products(path):
    data = json.loads(gzip.decompress(path.read_bytes()))
    rows = []
    for recipe in data["highest_weight_recipes"]:
        value = Counter()
        for factors,coefficient in recipe:
            for key,n in cartan(tuple(map(tuple,factors))).items():
                value[key] += coefficient*n
        rows.append(clean(value))
    return rows


def rank(rows):
    keys = sorted({w for row in rows for w in row})
    if not keys:
        return 0
    return fmpz_mat([[r.get(w,0) for w in keys] for r in rows]).rank()


def products(root):
    first = cartan_products(root/"products_1.json.gz")
    third = cartan_products(root/"products_3.json.gz")
    assert rank(first)==len(first)==70 and rank(third)==len(third)==14
    raised = []
    for row in first:
        raised.append({(s,i,k):v for s in (0,1) for i in (1,2)
                       for k,v in act(row,s,i,True).items()})
    joint = rank(raised)
    assert joint==70
    spectrum = {}
    K = [returned(r,act) for r in third]
    for eigen in (0,9,16):
        spectrum[eigen] = len(third)-rank([difference(k,r,eigen) for k,r in zip(K,third)])
    assert sum(spectrum.values())==14
    # Test the annihilating polynomial before projection onto any chosen coordinates.
    for v in K:
        w = difference(returned(v,act),v,9)
        assert not difference(returned(w,act),w,16)
    return dict(grade_one_Cartan_dimension=70, joint_flavour_two_three_raising_rank=joint,
                grade_three_Cartan_dimension=14, grade_three_return_eigenspaces=spectrum)


def source(page_path, cycle_path, out, workers):
    page = load_page(page_path)
    cycle = json.loads(gzip.decompress(cycle_path.read_bytes()))
    assert cycle["page_sha256"]==hashlib.sha256(page_path.read_bytes()).hexdigest()
    rows = centre_basis(page)
    coefficients = Counter()
    for i,scale in cycle["cycle_basis_coefficients"]:
        for w,n in rows[i].items():
            coefficients[w] += Fraction(scale)*n
    row = clean(coefficients)
    assert len(row)==4264 and all(not any(i in theta for i in (1,2))
                                  for w in row for tr in w for _,_,theta in tr)
    second = returned(row,substitute)
    library = m.library(out,6,m.P)
    initial = pair(row,page,library,workers)
    k = pair(second,page,library,workers)
    assert initial==1 and k in (0,9,16)
    return dict(cycle_terms=len(row), returned_terms=len(second), prime=m.P,
                initial_pairing=initial,returned_pairing=k,
                possible_integer_eigenvalues=[0,9,16], return_eigenvalue=k)


def classical(root):
    archive = json.loads((root/"adjacent_ranks.json").read_bytes())
    grades = {v["data"]["sector"][-1]:v["data"] for v in archive.values()}
    bounds = []
    weight = (1,1,7,2,2)
    assert len(audit_rules())==12
    terms = []
    for ds in (0,1):
        for perm in permutations(range(3)):
            w = (1+ds,1-ds)+tuple(weight[i+2]+perm[i]-i for i in range(3))
            terms.append([w,(-1)**ds*sign(perm),len(words(w))])
    assert sum(s*n for _,s,n in terms)==0
    for t in range(4):
        data = grades[t]
        y = 13-t
        by_length = {r["Y"]:r for r in data["blocks"]}
        if t==1:
            n = sum(dimensions(weight,y,6,"Sp").values())
            before = by_length[y-1]["rank_out"]
            after = sum(dimensions(weight,y+1,6,"Sp").values())
        else:
            row = by_length[y]
            n,before,after = row["dim"],row["rank_in"],row["rank_out"]
        bounds.append(dict(grade=t,length=y,cochains=n,incoming_rank_lower=before,
                           outgoing_rank_lower=after,cohomology_upper=n-before-after))
    assert [b["cohomology_upper"] for b in bounds]==[0,70,0,14]
    return dict(stable_terminal_Weyl_terms=terms, bounds=bounds)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parent/"sp3_omega")
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--task",choices=("source","products"),required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.task=="source":
        record = source(a.root/"target.json.gz",a.root/"cycle.json.gz",a.out,a.workers)
    else:
        record = {**classical(a.root/"outgoing"),**products(a.root/"outgoing")}
    m.write(a.out/(a.task+".json"),record)
