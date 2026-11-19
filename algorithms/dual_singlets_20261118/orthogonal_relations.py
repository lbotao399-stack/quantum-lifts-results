"""Project orthogonal determinant relations at (2,2,2,2,2) to the spin/flavour singlet quotient."""
import argparse, gzip, hashlib, itertools, json, math, random, subprocess, tempfile
from collections import Counter
from functools import cache
from pathlib import Path
from flint import fmpz_mat, fmpq_mat
from terminal_commutators import words, parity
HERE=Path(__file__).resolve().parent
W=(2,2,2,2,2)

def norm(parts):
    key=lambda v:(sum(v),v)
    odd=[key(v) for v in parts if parity(v)]
    if len(set(odd))!=len(odd):return 0,()
    return (-1)**sum(v>w for i,v in enumerate(odd) for w in odd[i+1:]), tuple(sorted(parts,key=key))

def lower(word,i,j):
    out=Counter()
    for k,v in enumerate(word):
        if not v[i] or (j<2 and v[j]):continue
        t=list(v);t[i]-=1;t[j]+=1;s,x=norm(word[:k]+(tuple(t),)+word[k+1:]);out[x]+=s*v[i]
    return out

def coinvariants():
    basis=words(W);rows=[]
    for i,j in [(0,1),(2,3),(3,4)]:
        w=list(W);w[i]+=1;w[j]-=1
        for v in words(tuple(w)):
            z=lower(v,i,j);rows.append([z[b] for b in basis])
    m=fmpz_mat(rows); n,k=m.nullspace();assert len(basis)==224 and m.rank()==218 and k==6
    c=[]
    for j in range(k):
        vals=[int(n[i,j]) for i in range(len(basis))];d=math.gcd(*vals);c.append([x//d for x in vals]);
    cov=fmpz_mat(c).transpose();assert (m*cov).is_zero() and cov.rank()==6
    return basis,cov

def pack(g):return sum(v*8**i for i,v in enumerate(g))

def table(path,basis,cov):
    with path.open('w') as f:
        print(len(basis),file=f)
        for b,row in zip(basis,cov.tolist()):print(len(b),*[pack(g) for g in b],*[int(x) for x in row],file=f)

def make_case(rng):
    labs=list(range(5))*2; rng.shuffle(labs)
    layout=list(range(20));rng.shuffle(layout)
    return labs,layout

def run(exe,tab,cases):
    inp=''.join(' '.join(map(str,ls+lay))+'\n' for ls,lay in cases)
    p=subprocess.run([str(exe),str(tab)],input=inp,text=True,stdout=subprocess.PIPE,check=True)
    ans=[list(map(int,line.split())) for line in p.stdout.splitlines()];assert len(ans)==len(cases);return ans

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--attempts',type=int,default=180);ap.add_argument('--output',type=Path);ap.add_argument('--replay',type=Path);a=ap.parse_args();assert __debug__
    basis,cov=coinvariants(); print('singlet quotient',len(basis),cov.ncols(),flush=True)
    with tempfile.TemporaryDirectory() as td:
        t=Path(td);exe=t/'relations';tab=t/'columns';table(tab,basis,cov);subprocess.run(['c++','-O3','-std=c++17',str(HERE/'orthogonal_relations.cpp'),'-o',str(exe)],check=True)
        cases=[];selected=[];rows=[];coeff=[];rng=random.Random(202611190555);tried=0
        if a.replay:
            rec=json.loads(a.replay.read_bytes());cases=[(b['species'],b['slots']) for b in rec['independent_relations']]
        while tried<(len(cases) if a.replay else a.attempts) and len(rows)<5:
            batch=cases if a.replay else [make_case(rng) for _ in range(min(24,a.attempts-tried))]
            for (ls,lay),v in zip(batch,run(exe,tab,batch)):
                tried+=1
                raw=v[:-6];x=v[-6:];assert len(raw)==len(basis) and (fmpz_mat([raw])*cov).tolist()[0]==x
                if any(x) and fmpz_mat(rows+[x]).rank()>len(rows):
                    selected.append(dict(species=ls,slots=lay,formal_trace_coefficients=raw,coinvariant_coordinates=x));rows.append(x);print('rank',len(rows),'at',tried,flush=True)
                    if len(rows)==5:break
            if a.replay:break
        out=dict(weight=W,groups=['SO(7)'],formal_weight_space=len(basis),lowering_rank=218,formal_singlets=6,formal_basis_sha256=hashlib.sha256(tab.read_bytes()).hexdigest(),independent_relations=selected,rank=len(rows),singlet_upper_bound=6-len(rows));
        if a.output:a.output.write_text(json.dumps(out,indent=2)+'\n')
        if a.replay:assert json.loads(json.dumps(out))==rec
        print('upper bound',6-len(rows),flush=True)
if __name__=='__main__':main()
