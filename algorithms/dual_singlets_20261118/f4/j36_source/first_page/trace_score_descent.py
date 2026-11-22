"""Eliminate the corrected quantum row against the five surviving trace coordinates."""
from functools import reduce
import gzip
import json
from math import gcd

from sympy import Matrix, Rational

from flag_boundaries import HERE
from trace_source_cycles import WEIGHT, monomials, quartic_pivots, reduce_tree


def calculate():
    data = json.loads(gzip.decompress((HERE/"trace_source_cycles.json.gz").read_bytes()))
    score = json.loads((HERE/"trace_source_score.json").read_bytes())
    src = [tuple(map(tuple,m)) for m in data["monomials"]]
    piv,index = quartic_pivots()
    columns = [reduce_tree(m,piv,index) for m in src]
    keys = sorted(set().union(*(v.keys() for v in columns)))
    split = Matrix(len(keys),len(src),lambda i,j: columns[j].get(keys[i],0))
    def expand(columns):
        out = Matrix.zeros(len(src),len(columns))
        for j,pairs in enumerate(columns):
            for i,v in pairs:
                out[i,j] = v
        return out
    singlet,cycle = [expand(data[k]) for k in ("singlet_columns","cycle_columns")]
    row = Matrix([list(map(Rational,score["monomial_scores"]))])
    terminal = split*singlet
    quantum = row*singlet
    assert [str(x) for x in quantum] == score["singlet_scores"]
    assert split*cycle == Matrix.zeros(len(keys),14) and row*cycle == Matrix.zeros(1,14)
    assert singlet.rank()==19 and cycle.rank()==14 and any(quantum)
    independent_rows = list(terminal.T.rref()[1])
    independent_columns = list(terminal.rref()[1])
    matrix = terminal.extract(independent_rows,independent_columns)
    assert matrix.rows == matrix.cols == 5 and matrix.det()
    functional = quantum[:,independent_columns]*matrix.inv()
    assert functional*terminal[independent_rows,:]==quantum
    numerators = [int(x.p) for x in functional]
    assert all(x.q==1 for x in functional)
    scale = reduce(gcd,map(abs,numerators))
    words = monomials(WEIGHT)
    return dict(singlets=19,classically_correctable_cycles=14,tree_rank=terminal.rank(),
                augmented_rank=terminal.col_join(quantum).rank(),
                singlet_basis_columns=independent_columns,
                terminal_quotient_trace_coordinates=[words[keys[i]] for i in independent_rows],
                minor=[[str(v) for v in r] for r in matrix.tolist()],minor_determinant=str(matrix.det()),
                eta_Q1_values_at_selected_columns=[str(quantum[i]) for i in independent_columns],
                eta_Q1_on_terminal_coordinates=list(map(str,functional)),
                common_factor=str(scale),primitive_row=[v//scale for v in numerators])


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
