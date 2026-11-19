"""Compare native rows with rational trace differentiation and independent diagonal exterior products."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import numpy as np
import joint as j


def reduce_coefficients(expr):
    return {w:v for w,c in expr.items() if (v:=int(c.numerator)*pow(int(c.denominator),-1,j.P)%j.P)}


def wedge(left,right):
    out=defaultdict(int)
    for a,x in left.items():
        for b,y in right.items():
            if a&b:continue
            sign=(-1)**sum((b&((1<<i)-1)).bit_count() for i in range(a.bit_length()) if a>>i&1)
            out[a|b]+=sign*x*y
    return {m:c for m,c in out.items() if c}


def exact_diagonal(expr,frame):
    """Lift first two diagonal residues to integers and set the third to minus their sum."""
    masks=[sum(1<<v for v in s) for s in combinations(range(frame.k),frame.degrees[0])]
    result=[]
    for point in range(frame.P):
        @lru_cache(None)
        def letter(L,k):
            odd,r=frame.L[L]; a=r[point]; x=a[:,0,0];y=a[:,1,1];assert np.all((x+y+a[:,2,2])%j.P==0) and all(not a[:,h,l].any() for h in range(3) for l in range(3) if h!=l)
            return {(1<<b if odd else 0):(int(x[b]) if k==0 else int(y[b]) if k==1 else -int(x[b])-int(y[b])) for b in range(len(x))}
        @lru_cache(None)
        def trace(t):
            out=defaultdict(int)
            for k in range(3):
                cur={0:1}
                for L in t:cur=wedge(cur,letter(L,k))
                for m,c in cur.items():out[m]+=c
            return dict(out)
        total=defaultdict(Fraction)
        for w,c in expr.items():
            cur={0:1}
            for t in w:cur=wedge(cur,trace(t))
            for m,v in cur.items():total[m]+=c*v
        result.extend(total[m] for m in masks)
    return result


def control(root):
    assert j.inputs(root)==j.INPUTS; plan=j.load_plan();maps=j.runtime(root)
    from sp_full_profiles import profile
    from sp_full_projected_basis import evaluate as dense
    from native_su_parallel_v1 import evaluate
    sc=maps.sector(tuple(plan['mu']),9); rows=[]
    for i,k in ((2,34),(3,29),(30,5),(35,3)):
        s=plan['sources'][i];word=j.recipes(root,s)[k];init=maps.algebra.parity_project({word:1},-1);expr=maps.algebra.projection.project(init,sc);assert reduce_coefficients(expr)==maps.project(init,sc,j.P)
        d=maps.algebra.differential(expr);assert not maps.algebra.differential(d);pieces=defaultdict(dict);lookup={tuple((tuple(l),n) for l,n in plan['targets'][t]['profile']):t for t in s['targets']}
        for w,c in d.items():pieces[lookup[profile(tuple(L for t in w for L in t))]][w]=c
        fr=j.frames(maps,plan,i,2);values=j.run_rows(maps,plan,i,[word],fr,2);stats=[]
        for t in s['targets']:
            key=f'q{t:02d}';f=fr[key][0];rhs=dense([reduce_coefficients(pieces[t])],f)[0];assert np.array_equal(rhs,values[key][0]);desc=fr[key][1]
            # Exact cancellation after assigning diagonal matrices to the differentiated traces.
            diag=maps.Frame(3,1,f.k,f.degrees,list(f.L),seed=991+t,p=j.P,cartan=True)
            assert all(c==0 for c in exact_diagonal(pieces[t],diag))
            stats.append(dict(target=t,terms=len(pieces[t]),coordinates=len(rhs),nonzero=int(np.count_nonzero(rhs)),sha256=j.digest(rhs.tolist())))
        f=fr['c'][0];exact=exact_diagonal(expr,f);m=[int(c.numerator)*pow(int(c.denominator),-1,j.P)%j.P for c in exact];assert m==values['c'][0].tolist() and any(m)
        rows.append(dict(source=i,seed=k,fermions=s['fermions'],seed_file_sha256=s['recipe_sha256'],rational_projected_terms=len(expr),rational_Q0_terms=len(d),Q0_squared_terms=0,targets=stats,cartan_coordinates=len(m),cartan_nonzero=sum(bool(x) for x in m),first_exact_diagonal=str(exact[0]),first_diagonal_residue=m[0],cartan_sha256=j.digest(m)))
    assert all(any(t['nonzero'] for t in x['targets']) for x in rows)
    return dict(prime=j.P,plan_sha256=j.PLAN,inputs_sha256=j.INPUTS,driver_sha256=j.sha(j.__file__),control_sha256=j.sha(__file__),rows=rows,scope='Rational trace projection and differentiation before reduction; exact integer diagonal Grassmann multiplication; both compared with the joint worker on two evaluation points per frame. This is a map implementation control, not an obstruction.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();r=control(a.root);j.save(a.output,r);print([(v['source'],v['rational_Q0_terms'],v['first_diagonal_residue']) for v in r['rows']])
