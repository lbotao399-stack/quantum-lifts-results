"""Test whether a single target covector excludes incoming images on every grade."""
import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from dual_charge_bands import coordinates
from filtered_rank import matmul
from sp3_omega_actions import centre_basis, load_page
from trace_descent import action, highest_weights, settings

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"so7_cartan_target_20261118"))
from reduce_incoming import Echelon


def all_rows(weight, length):
    sources = highest_weights(weight, length)
    return [r for degree in sorted(sources) for r in sources[degree]]


def run(page_path, out, workers):
    settings()
    page = load_page(page_path)
    out.mkdir(exist_ok=True, parents=True)
    target = tuple(page["weight"])
    library = m.library(out, 6, m.P)
    letters = sorted({l for content in m.contents(target, 9) for l in content})
    probes = page["centre_frames"]["4"]["samples"]
    def values(rows):
        return coordinates(rows, 4, letters, probes, 6, "Sp", library, workers)
    parents = all_rows(target, 8)
    before = page["incoming"]
    selected = [m.q(parents[i]) for i in before["rows"]]
    boundary = np.vstack([values(selected[i:i+192]) for i in range(0,len(selected),192)])
    jj = before["columns"]
    free = sorted(set(range(boundary.shape[1])).difference(jj))
    square = nmod_mat(boundary[:,jj].tolist(), m.P)
    assert int(square.det()) == before["determinant"] != 0
    projection = np.asarray(square.solve(nmod_mat(boundary[:,free].tolist(), m.P)).tolist(), np.int64)
    def quotient(v):
        return (v[:,free] - matmul(v[:,jj], projection, m.P)) % m.P
    assert not np.any(quotient(boundary))
    # Check the quotient on every formal predecessor, including finite-colour identities.
    for i in range(0,len(parents),192):
        assert not np.any(quotient(values([m.q(r) for r in parents[i:i+192]])))
    cycle = {}
    for row, c in zip(centre_basis(page), page["representative"]["cycle_basis_coefficients"], strict=True):
        if c:
            for w,n in row.items():
                cycle[w] = (cycle.get(w,0)+c*n) % m.P
    target_row = quotient(values([{w:c for w,c in cycle.items() if c}]))[0]
    assert np.any(target_row)
    equations = Echelon(len(free), rhs=True)
    equations.add(np.asarray([list(target_row)+[1]], np.int64), [(-1,0)])
    result = dict(prime=m.P, target_page_sha256=hashlib.sha256(page_path.read_bytes()).hexdigest(),
                  target_dimension=boundary.shape[1], boundary_rank=len(jj), quotient_dimension=len(free),
                  all_classical_parents=len(parents), orders=[])
    for order in (1,2,3):
        t = 3-order
        weight = (t,t,8-t,3-t,3-t)
        length = 8+order
        contents = list(m.contents(weight,length))
        maximum = max((sum(sorted((len(l[2]) for l in content), reverse=True)[:order+1])
                       for content in contents), default=0)
        record = dict(order=order, weight=weight, content_profiles=len(contents),
                      maximum_selected_flavour_degree=maximum)
        result["orders"].append(record)
        if maximum < 3*order:
            record["identically_zero_by_flavour_degree"] = True
            continue
        source = all_rows(weight,length)
        record["formal_rows"] = len(source)
        for lo in range(0,len(source),192):
            actions = [action(r,order,6,"Sp") for r in source[lo:lo+192]]
            data = quotient(values(actions))
            equations.add(np.column_stack((data,np.zeros(len(data),dtype=np.int64))),
                          [(order,lo+i) for i in range(len(data))])
            record["checked_rows"] = lo+len(data)
            if equations.witness is not None:
                record["conflict"] = [[*k,c] for k,c in sorted(equations.witness.items())]
                result["target_only_solution"] = False
                m.write(out/"target_only.json.gz",result)
                print("combined images contain target at order",order,"row",lo+len(data),flush=True)
                return result
        record["system_rank"] = equations.rank
        print("through order",order,"rank",equations.rank,flush=True)
    coeff = equations.solution()
    functional = np.zeros(boundary.shape[1],dtype=np.int64)
    functional[free] = coeff
    functional[jj] = -(projection@coeff) % m.P
    assert len(free)*(m.P-1)**2 < 2**63
    assert not np.any(boundary@functional % m.P) and int(target_row@coeff % m.P)==1
    result.update(target_only_solution=True,target_functional=functional.tolist())
    m.write(out/"target_only.json.gz",result)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    assert __debug__
    run(args.page,args.out,args.workers)
