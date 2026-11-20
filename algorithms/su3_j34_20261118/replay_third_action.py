"""Recompute Q3 O0 with both path orientations and evaluate it and Q0 Q3 O0 at SU(3).

The three-colour specialization is over the integers and the exterior algebra
on fresh generators. A nonzero rational coefficient survives every trace identity;
a zero at this assignment is never a proof of a polynomial identity.
"""
import argparse,gzip,hashlib,json,random
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,combinations
from pathlib import Path
from third_action import Q3d_multi,hits,fw,alg,read_source,HERE,SOURCE_SHA,cartan,audit
from fortuity.full.letters import is_ferm
from verify import add_scaled

def load(path):
    h=hashlib.sha256();out={}
    with gzip.open(path,'rb') as stream:
        for line in stream:
            h.update(line);c,w=json.loads(line);key=tuple(tuple((a,b,tuple(s)) for a,b,s in tr) for tr in w);assert key not in out and Fraction(c);out[key]=Fraction(c)
    return out,h.hexdigest()

def double(arg):
    w,c=arg;return {v:c*k for v,k in Q3d_multi(w,3,'SU',both_paths=True).items()}

def basis_values(groups,seed):
    """Integers and ordinary 3-by-3 products in ordered exterior monomials."""
    rng=random.Random(seed);alphabet=sorted({L for g in groups for L in g});degree=max(sum(is_ferm(L) for L in g) for g in groups);M={}
    for L in alphabet:
        for k in (range(degree) if is_ferm(L) else [-1]):
            m=[rng.randrange(-4,5) for _ in range(8)];m.append(-m[0]-m[4]);M[L,k]=tuple(m)
    @lru_cache(96000)
    def tr(word,per):
        i=iter(per);A=(1,0,0,0,1,0,0,0,1)
        for L in word:
            B=M[L,next(i) if is_ferm(L) else -1];A=tuple(A[j]*B[k]+A[j+1]*B[k+3]+A[j+2]*B[k+6] for j in (0,3,6) for k in (0,1,2))
        return A[0]+A[4]+A[8]
    proof=[]
    for g,terms in sorted(groups.items()):
        n=sum(is_ferm(L) for L in g);perms=[(p,(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))) for p in permutations(range(n))];val=Fraction(0)
        for w,c in terms.items():
            sizes=[sum(map(is_ferm,t)) for t in w];v=0
            for p,sg in perms:
                offset=0;z=sg
                for t,ns in zip(w,sizes):z*=tr(t,p[offset:offset+ns]);offset+=ns
                v+=z
            val+=v*c
        proof.append(dict(content=g,terms=len(terms),exterior_degree=n,coefficient=str(val)))
    return dict(seed=seed,matrices=[[L,k,list(m)] for (L,k),m in sorted(M.items())],exterior_order='ascending 0 through degree minus one',contents=proof)

def grouped(row):
    z=defaultdict(dict)
    for w,c in row.items():z[tuple(sorted(L for t in w for L in t))][w]=c
    return z

def run(files,receipt,workers,out):
    receipt=json.loads(receipt.read_bytes());rec=receipt['receipt'];assert receipt['sacct'].split('|')[2]=='COMPLETED' and receipt['sacct'].endswith('|32|0:0')
    for n,(size,d) in rec['files'].items():p=files/n;assert p.stat().st_size==size and hashlib.sha256(p.read_bytes()).hexdigest()==d
    report=json.loads((files/'direct.json').read_bytes());assert report==rec['summary'] and report['source_sha256']==SOURCE_SHA
    stored,digest=load(files/'Q3.jsonl.gz');assert digest==report['output_uncompressed_sha256'];src=read_source(HERE/'O34_C_even.json.gz');selected=[];Q0source={}
    for w,c in src.items():
        add_scaled(Q0source,fw.Q_multi(w),c)
        # A nonzero local Q3 requires three copies of every theta flavour among
        # four selected letters. With one or no field-strength and only scalars/
        # fermions, check that necessary condition independently of the path code.
        fs=[L[2] for t in w for L in t]
        allowed=any(all(sum(i in g for g in quad)>=3 for i in range(3)) for quad in combinations([g for g in fs if len(g)>=2],4))
        if allowed:selected.append((w,c))
    assert len(selected)==419==report['surviving_theta_filter']
    total={};nonzero=0;before=0
    with ProcessPoolExecutor(workers) as pool:
        for row in pool.map(double,selected,chunksize=3):
            nonzero+=bool(row);before+=len(row);add_scaled(total,row,1)
    assert total==stored and nonzero==report['nonzero_individual_images'] and before==report['distinct_terms_before_weighted_sum'];audit(total,-1);assert alg.conjugate(total)==total
    ca=cartan(total);assert json.loads(json.dumps(ca))==json.loads((files/'Q3_Cartan.json').read_bytes())
    d={};by=grouped(total)
    for w,c in total.items():add_scaled(d,fw.Q_multi(w),c)
    seed=346712;pt=basis_values(by,seed);closed=None
    if d:
        gr=grouped(d);key=min(gr,key=lambda g:(len(gr[g]),g));closed=basis_values({key:gr[key]},seed)
        if not Fraction(closed['contents'][0]['coefficient']):closed=basis_values(gr,seed)
    nonzero_witness=sum(bool(Fraction(u['coefficient'])) for u in pt['contents']);record=dict(job=rec['job'],source_sha256=SOURCE_SHA,source_Q0_formal_terms=len(Q0source),source_words=len(src),independent_flavour_filter=len(selected),both_path_orientations_checked=len(selected),nonzero_individual_images=nonzero,ordered_image_sha256=digest,formal_trace_terms=len(total),distinct_contents=len(by),nonzero_at_three_colours=nonzero_witness,traceless_integer_specialization=pt,Cartan_formal_terms=ca['formal_diagonal_trace_terms'],Cartan_coefficients=len(ca['terms']),Q0_of_direct_Q3_formal_terms=len(d),Q0_of_direct_Q3_three_colour_test=closed,scope='Direct cochain only. Formal Q0 words on the finite-rank source are not a closure test. The direct image does not solve either lower descent equation.')
    out.write_text(json.dumps(record,indent=2)+'\n');print(len(total),nonzero_witness,len(d),None if closed is None else [(c['terms'],c['coefficient']) for c in closed['contents']][:6])
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);h=HERE/'direct_third';a.add_argument('--files',type=Path,default=h);a.add_argument('--receipt',type=Path,default=h/'a8_5666003.json');a.add_argument('--workers',type=int,default=1);a.add_argument('--out',type=Path,required=True);assert __debug__;run(**vars(a.parse_args()))
