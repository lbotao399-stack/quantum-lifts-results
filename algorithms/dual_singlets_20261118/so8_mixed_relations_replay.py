"""Replay the seven SO(8) mixed determinants by independent directed-loop expansion."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from flint import fmpz_mat
import mixed_incoming as m
from mixed33_exact import read
from rank_four_terminal import quotient
from replay_orthogonal import pairing


def evaluate(row):
    return {w:n for w,n in pairing(row["species"],row["slots"],9).items() if n}


def run(archive,workers):
    rec=read(archive)
    basis,columns=quotient(tuple(rec["weight"]),rec["formal"])
    lookup={w:i for i,w in enumerate(basis)}
    selected=rec["relations"]
    projections=[]
    with ProcessPoolExecutor(workers) as pool:
        for original,terms in zip(selected,pool.map(evaluate,selected)):
            raw={lookup[w]:n for w,n in terms.items()}
            assert raw==dict(original["raw"])
            projected=[sum(n*int(columns[i,j]) for i,n in raw.items()) for j in range(columns.ncols())]
            assert projected==original["singlet"]
            projections.append(projected)
    piv=rec["minor"]["columns"]
    det=int(fmpz_mat([[r[j] for j in piv] for r in projections]).det())
    assert det==rec["minor"]["determinant"] and det and len(selected)+rec["cartan"]==rec["formal"]
    return dict(templates=[r["index"] for r in selected],permutations_per_template=rec["terms_per_template"],
                raw_coefficients=sum(len(r["raw"]) for r in selected),determinant=det,
                classical_dimension=rec["formal"]-len(selected))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",type=Path,default=Path(__file__).with_name("rank_four")/"so8_mixed"/"orthogonal_paths.json")
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert __debug__ and a.workers>0
    out=run(a.archive,a.workers);m.write(a.out,out);print(out)
