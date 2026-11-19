"""Give one finite-field counterexample per failed content of a FIXED vector.

Convert each rational residual to F_1000003, evaluate all seven-generator exterior
coefficients, and keep a nonzero coordinate with its deterministic seed. This
certifies rational nonvanishing for that fixed vector. Zero samples are never
accepted as rational identities. Check the expander independently with explicit
integer matrices, path products and Grassmann permutations.
"""
import argparse,hashlib,random,subprocess
from collections import defaultdict,Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
import content_census as cen
import finite_trace_identity as ft
s,l=cen.s,cen.l;P=1000003

def assignment(seed):
    def draw():
        nonlocal seed;seed^=(seed<<13)&0xffffffff;seed^=seed>>17;seed^=(seed<<5)&0xffffffff;return seed%23-11
    M={}
    for k,odd in enumerate(l.f.ODD):
        for i in range(3):
            for j in range(3):
                for h in (range(7) if odd else [-1]):M[k,h,i,j]=-M[k,h,0,0]-M[k,h,1,1] if i==j==2 else draw()
    return M

def read_data(pl,exe,seed):
    odd=bytes(i for i,q in enumerate(l.f.ODD) if q).hex();data={}
    with pl.open('rb') as st:
        res=subprocess.run([str(exe.resolve()),str(seed),odd],stdin=st,check=True,capture_output=True,text=True)
    for v in res.stdout.splitlines():
        g,n,*a=v.split();n=int(n);masks=[k for k in range(128) if k.bit_count()==n];vs=list(map(int,a));assert g not in data and len(vs)==len(masks) and all(0<=v<P for v in vs);data[g]=dict(zip(masks,vs))
    return data

def controls(exe,seed,d):
    a=assignment(seed);ws=['000101','000105','020400','050007','000105070410','000105','050607090a0d10'];odd=bytes(i for i,q in enumerate(l.f.ODD) if q)
    control=d/'controls';words=[(bytes.fromhex(w),) for w in ws]+[(bytes.fromhex('000105'),bytes.fromhex('000105')), (bytes.fromhex('0001'),bytes.fromhex('0005'))];data=[]
    for i,w in enumerate(words):
        parity=sum(l.f.ODD[j] for j in b''.join(w));data.append(f'{i} {parity} 1\n1 {len(w)} '+ ' '.join(s.hexes(w))+'\n')
    control.write_text(''.join(data));values=read_data(control,exe,seed);tested=[]
    for i,w in enumerate(words):
        n=sum(l.f.ODD[j] for j in b''.join(w));ds=sorted(values[str(i)]);inds=sorted({ds[0],ds[len(ds)//2],ds[-1]})
        for mask in inds:
            ix=[k for k in range(7) if mask>>k&1];indices=list(range(n));M={(t,h):tuple(tuple(a[t,ix[h] if h!=-1 else -1,j,k] for k in range(3)) for j in range(3)) for t in set(b''.join(w)) for h in (indices if l.f.ODD[t] else [-1])};v=ft.exterior_value(w,M,n);assert v%P==values[str(i)][mask];tested.append([s.hexes(w),mask,v%P])
    assert any(v for _,_,v in tested);return tested

def run(rhs,vec,audit,exe,out,work,seeds):
    with TemporaryDirectory(dir=work) as t:
        d=Path(t);handles={}
        def emit(g,w,val):
            if g not in handles:handles[g]=(d/g.hex()).open('w')
            v=Fraction(val);assert v.denominator%P;num=(v.numerator*pow(v.denominator,-1,P))%P;handles[g].write(f'{num} {len(w)} '+ ' '.join(s.hexes(w))+'\n')
        try:c=cen.inventory(rhs,vec,emit)
        finally:
            for f in handles.values():f.close()
        sets=c['contents'];pl=d/'rows';import shutil
        with pl.open('wb') as dst:
            for g in sorted(handles):
                dst.write(f'{g.hex()} {sum(l.f.ODD[i] for i in g)} {sets[g.hex()][2]}\n'.encode());
                with (d/g.hex()).open('rb') as src:shutil.copyfileobj(src,dst)
        gold=s.read(audit);assert gold['ordered_source_sha256']==c['source_sha256'];replayed=gold['contents']+gold['holdouts'];assert len({v['content'] for v in replayed})==len(replayed)==70;zero={g for g,v in sets.items() if not v[-1]}|{r['content'] for r in replayed if r['nonzero_monomials']==0};known_nonzero={r['content'] for r in replayed if r['nonzero_monomials']};proof={};checks={}
        for seed in seeds:
            checks[str(seed)]=controls(exe,seed,d);prob=read_data(pl,exe,seed);assert set(prob)=={k.hex() for k in handles}
            for g,row in prob.items():
                for mask,value in row.items():
                    if value:assert g not in zero;proof.setdefault(g,[seed,mask,value]);break
        unresolved=sorted(set(sets)-zero-set(proof)-known_nonzero);res=dict(rhs_sha=c['rhs_sha256'],vector_sha=c['source_sha256'],words=c['source_words'],prime=P,generators=7,assignments=seeds,checks=checks,complete_contents=500,exact_zero=sorted(zero),mathematically_nonzero=sorted(set(proof)|known_nonzero),nonzero_witness=proof,nonzero_witness_order=['seed','exterior_mask','residue_mod_prime'],unresolved=unresolved,scope='First equation of the specified fixed representative. A nonzero evaluation is conclusive over Q; a zero evaluation is inconclusive.');s.write(out,res);print('fixed vector; zero, nonzero, unresolved:',len(zero),len(res['mathematically_nonzero']),len(unresolved),flush=True);return res

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--vector',type=Path,required=True);a.add_argument('--audit',type=Path,required=True);a.add_argument('--exe',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--seeds',type=int,nargs='+',default=[817234,419071]);b=a.parse_args();assert __debug__;b.work.mkdir(exist_ok=True,parents=True);run(b.rhs,b.vector,b.audit,b.exe,b.out,b.work,b.seeds)
