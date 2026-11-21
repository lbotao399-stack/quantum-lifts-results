"""Return all nine-letter highest-weight traces in the Sp(3) omega block after four substitutions."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
from pathlib import Path

import mixed_incoming as m
from so8_level24_descendant import check_letter, substitute, substitution
from trace_descent import dimensions, highest_weights


STEPS = ((0,1),(1,1),(0,2),(1,2))


def test(row):
    result = row
    for spin, flavour in STEPS:
        result = substitute(result, spin, flavour, False)
    middle = len(result)
    for spin, flavour in reversed(STEPS):
        result = substitute(result, spin, flavour, True)
    assert result == {w: 72*n for w, n in row.items()}
    return middle, len(row)


def run(workers):
    weight, length = (3,3,5,0,0), 9
    rows = [r for group in highest_weights(weight,length).values() for r in group]
    alphabet = {l for row in rows for w in row for tr in w for l in tr}
    checks = 0
    for increase, steps in ((False,STEPS), (True,reversed(STEPS))):
        for spin, flavour in steps:
            for l in alphabet:
                check_letter(l, spin, flavour, increase)
                checks += 1
            alphabet |= {v for l in alphabet for n,v in [substitution(l,spin,flavour,increase)] if n}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        sizes = list(pool.map(test, rows, chunksize=8))
    digest=hashlib.sha256()
    for row in rows:
        digest.update(m.compact([[m.encode(w), n] for w,n in sorted(row.items())]))
    return dict(weight=weight,letters=length,formal_highest_weights=len(rows),
                formal_basis_sha256=digest.hexdigest(),classical_letter_checks=checks,
                symplectic_character=dimensions(weight,length,6,"Sp"),
                lowered_weight=[1,1,5,2,2],forward_steps=STEPS,
                returned_coefficient=72,zero_lowered_rows=sum(middle==0 for middle,_ in sizes),
                input_terms=sum(n for _,n in sizes),lowered_terms=sum(n for n,_ in sizes),
                incoming_derivative_grades_after_lowering={str(r):1-r for r in (1,2,3)})


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    result=run(args.workers)
    m.write(args.out,result)
    print(result,flush=True)
