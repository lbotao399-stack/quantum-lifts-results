#!/usr/bin/env python3
"""Standalone stdlib verifier of two SO15 classes and exact SO13 preimages.

--rebuild-lower additionally compiles the supplied C++17 contraction engine
and regenerates every integer row modulo the complete CRT certificate.
"""
import argparse,hashlib,json,subprocess,tempfile,time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb,factorial,isqrt,lcm,prod
from pathlib import Path
import minimal_scalar_certificates as minimal
import pfaffian_scalar_candidates as p
from so13_character_certificate import character,partitions,transpose,centralizer,predecessor_shapes

HERE=Path(__file__).resolve().parent

@lru_cache(None)
def lower_relation_trace(cycles):
    degree=sum(cycles)
    shapes=[transpose(tuple(2*n for n in reversed(parts))) for parts in partitions(degree) if len(parts)>13]
    if not shapes:return 0
    types=sorted(Counter(cycles).items());value=0
    for crossed in product(*(range(m+1) for _,m in types)):
        vc=[];coefficient=1
        for (k,m),a in zip(types,crossed):
            vc += [k]*(2*(m-a))+[2*k]*a;coefficient*=comb(m,a)
        value+=coefficient*sum(character(s,tuple(sorted(vc))) for s in shapes)
    value,remainder=divmod(value,2**len(cycles));assert not remainder
    return value

def lower_flavour(shape,distinguished=0):
    v=sum(Fraction(character(shape,c)*lower_relation_trace(tuple(sorted(c+(1,)*distinguished))),centralizer(c)) for c in partitions(sum(shape)))
    assert v.denominator==1 and v>=0
    return int(v)

def modular_rank(rows,prime):
    rows=[[v%prime for v in row] for row in rows];rank=0
    for j in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inv=pow(rows[rank][j],-1,prime);rows[rank]=[v*inv%prime for v in rows[rank]]
        for i in range(rank+1,len(rows)):
            c=rows[i][j]
            if c:rows[i]=[(v-c*w)%prime for v,w in zip(rows[i],rows[rank])]
        rank+=1
        if rank==len(rows):break
    return rank

def explicit_projection(shape,key,raw):
    probes=[]
    for line in (HERE/f'lower{key}_probe_expansions.dat').read_text().splitlines():
        if line.startswith('PROBE'):probes.append([]);continue
        c,w=line.split();probes[-1].append((int(c),tuple(tuple(map(int,t)) for t in w.rstrip('.').split('.'))))
    terminal=p.Coefficient(minimal.operators(shape)[1]);values=[]
    for probe in probes:
        value=Fraction()
        for c,word in probe:
            norm=Fraction(prod(factorial(n) for n in Counter(word).values()))
            for trace in word:
                orbit={v[i:]+v[:i] for v in [trace,trace[::-1]] for i in range(len(trace))}
                norm*=Fraction(2*len(trace),len(orbit))
            value+=c*norm*terminal(word)
        values.append(value)
    ratios=[Fraction(a,b) for a,b in zip(raw,values) if b]
    assert ratios and ratios[0] and len(set(ratios))==1
    ratio=ratios[0];assert all(ratio*b==a for a,b in zip(raw,values))
    return ratio,list(map(str,values))

def verify_lower(key):
    r=json.loads((HERE/f'SO13_{key}_EXACT_BOUNDARY.json').read_text());shape=tuple(r['shape'])
    source=sum(lower_flavour(s,1) for s in predecessor_shapes(shape));terminal=lower_flavour(shape)
    assert source==r['source_relation_dimension'] and terminal==r['terminal_relation_dimension']
    frame=r['frame_integer_rows'];target=r['all_singleton_Pf17_integer_row'];assert len(frame)==terminal
    cols=r['injective_terminal_probe_columns'];assert len(cols)==terminal
    assert modular_rank([[row[j] for j in cols] for row in frame],1000003)==terminal
    coefficients=list(map(Fraction,r['rational_boundary_coefficients']))
    denominator=lcm(*(v.denominator for v in coefficients))
    integers=[v.numerator*(denominator//v.denominator) for v in coefficients]
    assert all(sum(c*frame[i][j] for i,c in enumerate(integers))==denominator*target[j] for j in range(len(target)))
    assert modular_rank(frame[:len(coefficients)],1000003)==len(coefficients)
    primes=r['integer_CRT_primes'];assert prod(primes)>2**257
    assert all(all(n%d for d in range(2,isqrt(n)+1)) for n in primes+[r['extra_verification_prime']])
    assert max(abs(v) for row in frame+[target] for v in row)<2**256
    ratio,values=explicit_projection(shape,key,target[:2])
    expected={'773':Fraction(-4,3969),'953':Fraction(256,1091475)}[key];assert ratio==expected
    return r,dict(passed=True,shape=shape,source_dimension=source,terminal_dimension=terminal,
                  full_terminal_frame_rank=terminal,actual_Q0_rows_in_preimage=len(coefficients),
                  rational_probe_identities=len(target),raw_divided_by_explicit=str(ratio),
                  explicit_Pfaffian_probe_values=values,raw_probe_values=target[:2],
                  fortuitous_component_at_SO13='zero',scope='Zero connecting class; no assertion that the total finite cohomology class is zero.')

def rebuild(records):
    start=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='rank-birth-replay-',dir=HERE) as tmp:
        exe=Path(tmp)/'probe'
        subprocess.run(['c++','-std=c++17','-O3',str(HERE/'lower_rank_projection_probe.cpp'),'-o',str(exe)],check=True)
        def get(job):
            key,r,prime=job;count=len(r['all_singleton_Pf17_integer_row']);pairs=4 if key=='773' else 2
            ids=r['exact_Q0_image_basis_seed_indices'];bids=[0] if key=='773' else r['additional_actual_product_seed_indices']
            def command(n,mode,seeds):
                out=subprocess.run([str(exe),str(n),'3','8',str(count),'0',str(seeds),str(mode),str(prime),str(pairs)],text=True,capture_output=True,check=True)
                return [list(map(int,line.split()[1:])) for line in out.stdout.splitlines()]
            q=command(13,1,max(ids)+1);b=command(13,4,max(bids)+1);raw=command(16,0,1)[0]
            actual=[q[i] for i in ids]+[b[i] for i in bids]+[raw]
            expected=r['frame_integer_rows']+[r['all_singleton_Pf17_integer_row']]
            assert len(actual)==len(expected)
            assert all(v%prime==w for row,other in zip(expected,actual) for v,w in zip(row,other)),(key,prime)
            return len(actual)*count
        jobs=[(key,r,prime) for key,r in records.items() for prime in r['integer_CRT_primes']+[r['extra_verification_prime']]]
        with ThreadPoolExecutor(max_workers=4) as pool:entries=list(pool.map(get,jobs))
    return dict(passed=True,prime_case_replays=len(jobs),residue_entries_checked=sum(entries),
                proof='Integer integrality and the a priori 2^256 bound in RANK_BIRTH_PROOF.md, together with CRT modulus>2^257, certify every rebuilt integer.',seconds=time.monotonic()-start)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rebuild-lower',action='store_true');ap.add_argument('--out',type=Path,default=HERE/'RANK_BIRTH_VERIFIED.json');a=ap.parse_args()
    start=time.monotonic();upper=[];lower=[];records={}
    for key,shape in [('773',(7,7,3)),('953',(9,5,3))]:
        up=minimal.run(15,shape);upper.append(up)
        records[key],low=verify_lower(key);lower.append(low)
    result=dict(passed=True,field='Q',SO15=upper,SO13=lower,
                rank_support='For these fixed representatives the fortuitous component is supported exactly at SO(15) among odd orthogonal groups: zero connecting class for oddN<=13, nonzero at15, and not Q0-closed for oddN>=17.',
                integer_bound_proof='RANK_BIRTH_PROOF.md',seconds=time.monotonic()-start)
    if a.rebuild_lower:result['full_lower_rebuild']=rebuild(records)
    result['source_sha256']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['pfaffian_scalar_candidates.py','so13_character_certificate.py','so13_exact_certificate.py','minimal_scalar_certificates.py','lower_rank_projection_probe.cpp','verify_rank_birth.py']}
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'SO15_quotients':[r['exact_fortuitous_quotient_dimension'] for r in upper],'SO13_terminal_dimensions':[r['terminal_dimension'] for r in lower],'SO13_raw_to_explicit_ratios':[r['raw_divided_by_explicit'] for r in lower],'full_lower_rebuild':result.get('full_lower_rebuild'),'seconds':time.monotonic()-start},indent=2))

if __name__=='__main__':main()
