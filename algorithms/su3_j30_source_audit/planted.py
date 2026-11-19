"""Recover rational trace inputs with the production projection, Q0 maps and affine reducer."""
import argparse, hashlib, io, contextlib, json, pickle, sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import numpy as np

PRIMES=(1000003,1000033,1000037)
sha=lambda b:hashlib.sha256(b).hexdigest()

def mod(expr,p):
    out={w:int(c.numerator)*pow(int(c.denominator),-1,p)%p for w,c in expr.items()}
    return {w:c for w,c in out.items() if c}

def fresh_frame(maps,d,points,p):
    seed0=d['seed']+31000019;seed1=d['additional_seed']+31000019
    f=maps.Frame(3,points,d['k'],d['degrees'],d['alphabet'],seed=seed0,p=p)
    if d['additional_letters']:
        extra=maps.Frame(3,points,d['k'],d['degrees'],d['additional_letters'],seed=seed1,p=p); f.L.update(extra.L)
    return f

def run(tree,audit):
    assert __debug__
    runtime=tree/'runtime/work/sp3-native.nosync';sys.path.insert(0,str(runtime))
    import su3_second_correction_map_v1 as maps
    from sp_full_profiles import profile
    from sp_full_projected_basis import evaluate as dense
    from bounded_flint_rref import reduce_target
    from staged_sp2_first_correction import back_substitute
    from flint import nmod_mat
    plan=json.loads((runtime/'su3_second_correction_profiles_v1.json').read_bytes());i=8;source=plan['sources'][i];n=source['dimension'];assert n==30
    sector=maps.sector(tuple(plan['mu']),plan['source_length']);raw=(tree/'second_basis/shard0/sources08/independent_seeds.pkl').read_bytes()
    pin=next(r['seeds_sha256'] for r in audit['profiles'] if r['equation']==2 and r['index']==i); assert sha(raw)==pin
    words=pickle.loads(raw);assert len(words)==n
    coefficients=[Fraction((-1)**j*(j+2),2*j+5) for j in range(n)]
    exact={}; projected={p:[] for p in PRIMES};nterms=[]
    for word,c in zip(words,coefficients):
        initial=maps.algebra.parity_project({word:1},-1); rational=maps.algebra.projection.project(initial,sector);nterms.append(len(rational));assert rational
        for p in PRIMES:
            native=maps.project(initial,sector,p); assert native==mod(rational,p);projected[p].append(native)
        for w,x in rational.items():exact[w]=exact.get(w,0)+x*c
    exact={w:x for w,x in exact.items() if x};derivative=maps.algebra.differential(exact);assert derivative and not maps.algebra.differential(derivative)
    lookup={tuple((tuple(k),c) for k,c in plan['targets'][j]['profile']):j for j in source['targets']}
    pieces=defaultdict(dict)
    for w,c in derivative.items():pieces[lookup[profile(tuple(L for t in w for L in t))]][w]=c
    records=[]
    for p in PRIMES:
        planted=np.array([mod({0:c},p)[0] for c in coefficients],dtype=np.int64);stages=[]; blocks=[]; descs=[];target_records=[]
        for j in source['targets']:
            target=plan['targets'][j]; f,d=maps.target_frame(plan,j,p,points_override=4); descs.append(d)
            M=maps.evaluate_profile(projected[p],f,source['profile'],target['profile'],2).T.copy();b=dense([mod(pieces[j],p)],f)[0]; assert np.array_equal((M@planted)%p,b)
            with contextlib.redirect_stdout(io.StringIO()):stage=reduce_target(np.column_stack([M,b]).astype(np.int32),np.arange(n+1),stages,n,p)
            assert stage['consistent']; stages.append(stage); blocks.append(M);target_records.append(dict(target=j,training=len(b),nonzero_rhs=int(np.count_nonzero(b)),rank_increment=stage['rank']))
        M=np.concatenate(blocks);rank=int(nmod_mat(M.tolist(),p).rank());x,free=back_substitute(stages,n,p);delta=(x-planted)%p
        assert rank==sum(s['rank'] for s in stages)==n-len(free) and not np.any(x[free]) and not np.any((M@delta)%p)
        for row,j,d in zip(target_records,source['targets'],descs):
            f=fresh_frame(maps,d,7,p);target=plan['targets'][j];N=maps.evaluate_profile(projected[p],f,source['profile'],target['profile'],2).T; truth=dense([mod(pieces[j],p)],f)[0]
            assert np.array_equal((N@x)%p,truth) and not np.any(N@delta%p); row.update(held_out=len(truth),held_out_nonzero_rhs=int(np.count_nonzero(truth)),held_out_residual_nonzero=0)
        records.append(dict(prime=p,rank=rank,free=free,number_of_changed_source_coefficients=int(np.count_nonzero(delta)),particular=x.tolist(),targets=target_records))
    assert all(any(t['held_out_nonzero_rhs'] for t in r['targets']) for r in records)
    return dict(control='known rational input independent of modular solver',source_profile=i,seed_sha256=pin,source_dimension=n,source_coefficients=list(map(str,coefficients)),projected_trace_terms_per_recipe=nterms,rational_source_terms=len(exact),rational_Q0_terms=len(derivative),rational_Q0_squared_terms=0,tests=records,scope='All 30 full projected trace dictionaries in this profile at all three primes; rational Q0 applied before reduction. Does not check the physical right-hand side.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('third_tree',type=Path);p.add_argument('source_audit',type=Path);p.add_argument('output',type=Path);a=p.parse_args();r=run(a.third_tree.resolve(),json.loads(a.source_audit.read_bytes()));a.output.write_text(json.dumps(r,indent=2)+'\n');print([(x['prime'],x['rank'],x['free']) for x in r['tests']]);print('Rational trace terms',r['rational_source_terms'],r['rational_Q0_terms'])
