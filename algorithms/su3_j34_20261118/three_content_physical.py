"""Construct the three-content correction directly in three-colour polynomial coordinates.
Only select parent words with images supported in the middle content on this projection;
the charge space outside these three contents remains unrestricted.
"""
import argparse, gzip, json, hashlib, os, subprocess
from fractions import Fraction
from pathlib import Path
from collections import Counter
from functools import reduce
from math import lcm
import replay_formal_blocks as s
import formal_trace_blocks as f
import projected_finite_correction as p
import finite_trace_identity as ft
import colour_outgoing as out

CONTENTS=tuple(bytes.fromhex(v) for v in ('0000000000010405050618','0000000000010505050710','0000000001040505060810'))
MID=CONTENTS[1]; PRIME=1000003
h=lambda x:hashlib.sha256(x).hexdigest()

def make(work,rhs):
    target=p.rhs_rows(rhs,set(CONTENTS));assert {p.content(v) for v in target}=={MID}
    parents=sorted(set().union(*(f.parents(w) for w in target)))
    isolated=out.incident(MID)-set.union(*(out.incident(c) for c in CONTENTS if c!=MID))
    image={}
    for w in parents:
        if p.content(w) not in isolated:continue
        y={v:c for v,c in f.q(w).items() if p.content(v) in CONTENTS}
        assert all(p.content(v)==MID for v in y)
        if y:image[w]=y
    words=sorted(set(target).union(*(z.keys() for z in image.values())));idx={w:i for i,w in enumerate(words)};den=lcm(*(k.denominator for k in target.values())); src=list(image)
    code=lambda y:sorted([idx[w],int(v)] for w,v in y.items());assert all((v*den).denominator==1 for v in target.values());b=code({w:v*den for w,v in target.items()})
    z=dict(rhs_sha256=s.INPUT,contents=[c.hex() for c in CONTENTS],middle=MID.hex(),all_input_terms=len(target),incident_parents=len(parents),isolated_contents=[c.hex() for c in sorted(isolated)],source=[s.hexes(w) for w in src],word=[s.hexes(w) for w in words],image=[code(image[w]) for w in src],rhs=b,den=den,odd=bytes(i for i,v in enumerate(ft.ODD) if v).hex());s.write(work/'recipe.json.gz',z)
    (work/'words.txt').write_text(str(len(words))+'\n'+''.join(f"{len(w)} {' '.join(s.hexes(w))}\n" for w in words));print('source, words, denominator:',len(src),len(words),den,flush=True);return z

def columns(work,z,exe,jobs):
    from concurrent.futures import ThreadPoolExecutor
    import shutil
    assert 1<=jobs<=64; exe=str(exe.resolve()); data=(work/'words.txt').read_text().splitlines(); n=int(data[0]); assert n==len(data)-1
    def chunk(j):
        i,k=n*j//jobs,n*(j+1)//jobs; path=work/('words-%02d.txt'%j); dest=work/('polys-%02d.txt'%j);path.write_text(str(k-i)+'\n'+'\n'.join(data[i+1:k+1])+'\n');subprocess.run([exe,str(path),z['odd'],str(dest)],check=True);return dest
    with ThreadPoolExecutor(jobs) as pool, (work/'polys.txt').open('wb') as dest:
        for path in pool.map(chunk,range(jobs)):
            with path.open('rb') as src:shutil.copyfileobj(src,dest)
            path.unlink()

def load_columns(work,z):
    import numpy as np
    from scipy.sparse import csr_matrix
    ind=[0];keys=[];data=[]
    with (work/'polys.txt').open() as file:
        for line in file:
            a=np.fromstring(line,sep=' ',dtype=np.int64);assert len(a)==1+2*a[0] and (a[1::2][:-1]<a[1::2][1:]).all() and (a[2::2]!=0).all()
            ind.append(ind[-1]+a[0]);keys.append(a[1::2].copy());data.append(a[2::2].copy())
    assert len(ind)==len(z['word'])+1
    unique,packed=np.unique(np.concatenate(keys),return_inverse=True);poly=csr_matrix((np.concatenate(data),packed,ind),shape=(len(z['word']),len(unique)),dtype=np.int64)
    def sparse(rows):
        a,b,c=[],[],[0]
        for row in rows:
            for i,v in row:a.append(i);b.append(v)
            c.append(len(a))
        return csr_matrix((np.asarray(b,dtype=np.int64),a,c),shape=(len(rows),len(z['word'])))
    L=sparse(z['image']+[z['rhs']]); bound=sum(abs(int(v)) for v in L.data)*max(abs(int(v)) for v in poly.data);assert bound<2**63
    A=(L@poly).tocsr();A.eliminate_zeros();A.sort_indices();A=A.tocsc();print('rows,monomials,nonzeros, RHS nonzeros:',A.shape,A.nnz,A.getrow(A.shape[0]-1).nnz,flush=True);return A,unique

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('mode',choices=['prepare','expand','solve','verify']);a.add_argument('--work',type=Path,required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--exe',type=Path);a.add_argument('--jobs',type=int,default=1);args=a.parse_args();assert __debug__;w=args.work;w.mkdir(exist_ok=True,parents=True)
    if args.mode=='prepare':make(w,args.rhs)
    if args.mode=='expand':columns(w,s.read(w/'recipe.json.gz'),args.exe,args.jobs)
    if args.mode in ('solve','verify'):
        z=s.read(w/'recipe.json.gz');A,keys=load_columns(w,z)
        from physical_image_solve import solve, verify
        if args.mode=='solve':solve(w,z,A,keys)
        else:verify(w,z,A,keys)
