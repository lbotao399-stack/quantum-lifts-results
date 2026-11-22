"""Rational highest-weight products with independent ordinary G2 Cartan restrictions."""
import argparse
from collections import Counter
from pathlib import Path
import sys


def run(args):
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    import mixed_incoming as m
    import numpy as np
    from mixed33_exact import apply, exact
    from sp3_omega_products import all_products, raising
    from su3_bmn_parity_products import labelled
    from sp3_omega_product_replay import trace_factor
    from flint import fmpq_mat

    families = all_products(6,tuple(args.weight),args.length)
    recipes, polynomials, echelon = [],[],{}
    all_recipes = []
    for family, words in sorted(families.items()):
        if any(g[0] not in (2,6) for g in family):
            continue
        images = [raising(w) for w in words]
        _, _, basis = exact(images)
        for combination in basis:
            assert not apply(combination,images)
            recipe = [(words[i],c) for i,c in sorted(combination.items())]
            all_recipes.append(recipe)
            row = Counter()
            for spec,c in recipe:
                for w,n in labelled(spec,True,'G2').items():
                    row[w] += c*n
            reduced = {k:int(v)%m.P for k,v in row.items() if v}
            while reduced:
                k=min(reduced)
                if k not in echelon:
                    iv=pow(reduced[k],-1,m.P)
                    echelon[k]={j:(v*iv)%m.P for j,v in reduced.items()}
                    recipes.append(recipe)
                    polynomials.append(row)
                    break
                c=reduced[k]
                for j,v in echelon[k].items():
                    n=(reduced.get(j,0)-c*v)%m.P
                    if n:
                        reduced[j]=n
                    else:
                        reduced.pop(j,None)
            if len(recipes)==args.bound and not args.all:
                break
        if len(recipes)==args.bound and not args.all:
            break
    columns=sorted(echelon)
    rational=fmpq_mat([[int(row[k]) for k in columns] for row in polynomials])
    det=rational.det()
    assert det and len(recipes)<=args.bound
    factors={g for comb in recipes for spec,c in comb for g in spec}
    assert all(not m.q(trace_factor(g)) for g in factors)
    result=dict(weight=args.weight,length=args.length,cartan_rank=len(recipes),
                determinant=str(det),columns=columns,recipes=recipes,factors=len(factors),
                highest_combinations_visited=len(all_recipes))
    m.write(args.out,result)
    if args.all:
        m.write(args.all,dict(weight=args.weight,length=args.length,recipes=all_recipes))
    print('weight',args.weight,'Cartan rank',len(recipes),'highest rows',len(all_recipes),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--weight',type=int,nargs=5,required=True)
    p.add_argument('--length',type=int,required=True)
    p.add_argument('--bound',type=int,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--all',type=Path)
    args=p.parse_args()
    if not __debug__:
        raise RuntimeError('enable Python assertions')
    run(args)
