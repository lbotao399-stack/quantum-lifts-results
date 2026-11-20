"""Finite-colour terminal relations in the level-30 mixed singlet of SO(9) and USp(8)."""
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path
import argparse,random,subprocess

from flint import fmpz_mat

import mixed_incoming as m
from orthogonal_relations import pack
from rank_four_terminal import bounds, make_case, quotient, relation


def make_orthogonal(rng, weight, index):
    labels=[i for i,n in enumerate(weight) for _ in range(n)]
    rng.shuffle(labels)
    slots=list(range(2*len(labels)))
    rng.shuffle(slots)
    if index%4==0:
        # A separate two-letter trace multiplies a ten-letter identity.
        rest=slots[:]
        letters=rng.sample(range(len(labels)),2)
        extra=[2*letters[0],2*letters[1],2*letters[0]+1,2*letters[1]+1]
        slots=[x for x in rest if x not in extra]+extra
    return labels,slots


def batch(exe, table, group, cases):
    cmd=[str(exe),str(table),group]
    payload=''.join(' '.join(str(n) for n in s+layout)+'\n' for s,layout in cases)
    output=subprocess.run(cmd,input=payload,text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output)==len(cases)
    rows=[]
    for line in output:
        pairs=list(map(int,line.split()))
        assert len(pairs)%2==0
        rows.append(dict(zip(pairs[::2],pairs[1::2])))
    return rows


def run(work,workers,orthogonal,symplectic,seed):
    work=work.resolve()
    work.mkdir(exist_ok=True,parents=True)
    screen=[bounds(t) for t in range(6)]
    case=screen[3];weight=case["weight"]
    basis,cov=quotient(weight,case["formal"])
    columns=[[int(x) for x in row] for row in cov.tolist()]
    table=work/"basis.txt"
    with table.open("w") as f:
        print(len(basis),file=f)
        for w in basis:
            print(len(w),*(pack(g) for g in w),file=f)
    exe=work/"brauer"
    subprocess.run(["c++","-std=c++17","-O3",str(Path(__file__).with_name("rank_four_brauer.cpp")),"-o",str(exe)],check=True)
    for group,num,block in (("Sp",symplectic,128),("SO",orthogonal,3)):
        rng=random.Random(seed)
        selected=[];reduced=[];attempted=0;nonzero=0
        for start in range(0,num,block*workers):
            n=min(num-start,block*workers)
            cs=[(make_case(weight,rng,i) if group=="Sp" else make_orthogonal(rng,weight,i))
                for i in range(start,start+n)]
            chunks=[cs[i:i+block] for i in range(0,len(cs),block)]
            with ThreadPoolExecutor(workers) as pool:
                records=[r for part in pool.map(lambda x:batch(exe,table,group,x),chunks) for r in part]
            for idx,((species,slots),raw) in enumerate(zip(cs,records)):
                attempted+=1
                if not raw:continue
                projected=[sum(c*columns[i][j] for i,c in raw.items()) for j in range(case["formal"])]
                if not any(projected):continue
                nonzero+=1
                if fmpz_mat(reduced+[projected]).rank()>len(reduced):
                    if group=="Sp":
                        actual=relation(species,slots)
                        assert raw=={i:actual[w] for i,w in enumerate(basis) if w in actual}
                    reduced.append(projected)
                    selected.append(dict(trial=start+idx,species=species,slots=slots,
                        raw=sorted(raw.items()),singlet=projected))
                    print(group,"independent",len(reduced),"at trial",start+idx,flush=True)
                    if len(reduced)==case["formal"]-case["cartan"]:
                        break
            minor={}
            if reduced:
                jj,det=next((c,int(sq.det())) for c in combinations(range(case["formal"]),len(reduced))
                       for sq in [fmpz_mat([[r[j] for j in c] for r in reduced])] if sq.det())
                minor=dict(columns=jj,determinant=det)
            m.write(work/(group+".json"),dict(level=30,terminal=screen,weight_dimension=len(basis),group=group,
                weight=weight,free_slots=10,terms_per_relation=945 if group=="Sp" else 3628800,
                calculated_through=start+n,attempted=attempted,nonzero_singlet_projections=nonzero,
                relations=selected,minor=minor,classical_lower=case["cartan"],classical_upper=case["formal"]-len(reduced)))
            if len(reduced)==case["formal"]-case["cartan"]:break
        print(group,"bounds",case["cartan"],case["formal"]-len(reduced),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--orthogonal",type=int,default=900)
    p.add_argument("--symplectic",type=int,default=25000)
    p.add_argument("--seed",type=int,default=1627924)
    assert __debug__
    run(**vars(p.parse_args()))
