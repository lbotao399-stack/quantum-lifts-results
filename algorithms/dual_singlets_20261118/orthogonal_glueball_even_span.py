"""Find an ordered-word minor of the even-colour relations and test their second pairing."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from multiprocessing import get_context

from flint import fmpq_mat
from orthogonal_glueball_even_colours import evaluate
from orthogonal_glueball_relation_span import determinant, probes, relation_dimension
from orthogonal_glueball_second_sources import ranks


def probe_row(row):
    relation = determinant(POWER, row['spec'])
    values = [relation.coefficient(w) for w in PROBES]
    relation.clear()
    return list(map(str,values))


def span(case, workers, attempts, word_count):
    global POWER, PROBES
    POWER = case['power']
    colours = case['colours']
    scores = [(tuple(map(tuple,w)),Fraction(c)) for w,c in case['dual']['scores']]
    PROBES = probes(POWER, case['dual']['scores'], word_count)
    dimension = relation_dimension(colours,POWER,2)
    assert dimension == case['flavour_relations']
    # A basis of all balanced-dotted, flavour-singlet relations covers the dotted singlet.
    selected = []
    values = []
    rows = case['rows'][:attempts]
    examined = 0
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for start in range(0,len(rows),workers):
            chunk = rows[start:start+workers]
            for row,v in zip(chunk,pool.map(probe_row,chunk)):
                examined += 1
                if fmpq_mat(values+[v]).rank() > len(values):
                    values.append(v)
                    selected.append(row)
                if len(values) == dimension:
                    break
            if len(values) == dimension:
                break
    rank = len(values)
    if not rank:
        raise ArithmeticError('no nonzero relations on probes')
    red,_ = fmpq_mat(values).rref()
    columns = [next(j for j in range(red.ncols()) if red[i,j]) for i in range(rank)]
    minor = [[row[j] for j in columns] for row in values]
    det = fmpq_mat(minor).det()
    assert det
    for row in selected:
        recomputed = evaluate((row['index'],row['spec'],colours,POWER,scores))
        assert recomputed == row
    scores_ranks = ranks(selected)
    return dict(colours=colours,power=POWER,balanced_flavour_relations=dimension,rank=rank,
                lower_full_weight_relations=[relation_dimension(colours,POWER,k,False) for k in (0,1)],
                examined=examined,probes=[PROBES[j] for j in columns],minor=minor,
                determinant=str(det),basis_rows=selected,**scores_ranks)
