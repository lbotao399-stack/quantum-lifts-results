"""Verify a rank-three Cartan relation between graviton products in SO(7) and USp(6)."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path

from flint import fmpq_mat
import mixed_incoming as m
from gaugino_primary_order import charge
from so7_bmn_633 import diagonal, export_trace, minor
from so7_hair_first import product
from so8_level24_descendant import raising_images
from sp3_omega_product_replay import check_factor, trace_factor
from sp3_omega_products import cartan
from su2_levi_images import check_convention

HERE=Path(__file__).resolve().parent


def recipes(raw):
    return [(tuple(map(tuple,factors)),n) for factors,n in raw]


@lru_cache(maxsize=1600)
def trace_monomial(factors):
    if not factors:
        return {():1}
    return product(trace_monomial(factors[:-1]),trace_factor(factors[-1]))


def expand(row, evaluate):
    out=Counter()
    for factors,n in row:
        for w,c in evaluate(factors).items():
            out[w]+=n*c
    return {w:c for w,c in out.items() if c}


def reconstruct():
    base=json.loads(gzip.decompress((HERE/"products.json.gz").read_bytes()))
    raw=(HERE/"product_recipe.json.gz").read_bytes()
    chosen=recipes(json.loads(gzip.decompress(raw)))
    reference=list(map(recipes,base["highest_weight_recipes"]))
    columns=[tuple(tuple(map(tuple,s)) for s in c) for c in base["cartan_columns"]]
    vectors=[expand(row,cartan) for row in reference]
    restriction=expand(chosen,cartan)
    basis=fmpq_mat([[str(v.get(c,0)) for c in columns] for v in vectors])
    assert int(basis.det())==int(json.loads((HERE/"replay.json").read_bytes())["determinant"])
    sol=basis.transpose().solve(fmpq_mat([[str(restriction.get(c,0))] for c in columns]))
    coefficients=[Fraction(str(sol[i,0])) for i in range(sol.nrows())]
    combined=Counter(dict(chosen))
    assert len(combined)==len(chosen)
    for mult,row in zip(coefficients,reference,strict=True):
        if mult:
            for factors,c in row:
                combined[factors]-=mult*c
    combined=[(k,v) for k,v in sorted(combined.items()) if v]
    factors={g for w,c in combined for g in w}
    for g in sorted(factors):
        check_factor(g)
    traces=expand(combined,trace_monomial)
    assert charge(traces)==((0,0,6,3,3),10)
    assert not m.q(traces) and all(not any(r.values()) for r in raising_images(traces))
    assert not expand(combined,cartan) and not diagonal(traces)
    rank_four=expand(combined,lambda row:cartan(row,rank=4))
    check_convention()
    result=dict(level=24,groups=["SO(7)","USp(6)"],weight=[6,3,3],length=10,
                recipe_sha256=hashlib.sha256(raw).hexdigest(),seed_recipe_terms=len(chosen),
                subtracted_recipes=[[i,str(c)] for i,c in enumerate(coefficients) if c],
                graviton_monomials=len(combined),verified_single_factors=len(factors),trace_terms=len(traces),
                formal_classical_variation_terms=0,
                rank_three_Cartan_terms=0,rank_four_Cartan_terms=len(rank_four))
    if rank_four:
        w,c=min(rank_four.items())
        result["first_rank_four_Cartan_coefficient"]=dict(orbit=w,coefficient=str(c))
    return traces,result


def run(native_path, rank_path):
    traces,result=reconstruct()
    n,den,digest=export_trace(traces,native_path)
    assert n==result["trace_terms"]
    result.update(common_denominator=den,native_trace_sha256=digest)
    if rank_path:
        record=json.loads(rank_path.read_bytes())
        assert (record["gauge"],record["prime"],record["sector"],record["length"],record["seed"])==(
            "SO7",1000033,[6,3,3,0,0],10,778117)
        y0,y1,y2=record["cochains"]
        assert (y0,y1,y2)==(1863,1837,772)
        ranks=[minor(record,"incoming",y1,y0,1000033),minor(record,"augmented",y1,y0+1,1000033)]
        assert ranks==[1151,1152] and minor(record,"outgoing",y2,y1,1000033)==665
        result.update(prime=1000033,augmented_ranks=ranks,augmented_minor_det=record["augmented"]["det"])
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--export",type=Path,required=True)
    p.add_argument("--rank",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:p.error("enable assertions")
    a.out.write_text(json.dumps(run(a.export,a.rank),indent=2)+"\n")
