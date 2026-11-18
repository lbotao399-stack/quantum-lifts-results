"""Full-letter SO(7)/USp(6) highest-weight spaces at the two specified physical charges."""
import argparse, json, sys
from collections import Counter
from functools import lru_cache
from itertools import permutations
from pathlib import Path
HERE=Path(__file__).resolve().parent
RUNTIME=HERE.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
sys.path.insert(0,str(RUNTIME))
from screen_gaugino_source_dimensions_v1 import content_patterns
from exact_dims import _multiplicity
from sp_dimensions import multiplicity_parts
from g2_full_classical_screen_v1 import sector

@lru_cache(None)
def dimensions(group,mu,length):
    dims=Counter();rho=(2,1,0)
    for spin,ss in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            w=(mu[0]+spin,mu[1]-spin)+tuple(mu[2+i]+rho[i]-rho[perm[i]] for i in range(3))
            sign=ss*(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            for b,f,num in content_patterns(w,length):
                m=_multiplicity(b,f,0,7) if group=='SO_odd' else multiplicity_parts(b,f,3)
                dims[sum(f)]+=sign*num*m
    assert all(v>=0 for v in dims.values())
    return {d:n for d,n in sorted(dims.items()) if n}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();rows=[]
    for h in (4,0,3,1,2):
        w=(h,h,4-h,4-h,4-h);L=12-h
        assert sector(w,L)==(0,0,8,8,8) and sector(w,L-1)==(1,1,8,8,8)
        for group in ('SO_odd','Sp'):
            ds={n:dimensions(group,w,n) for n in range(L-2,L+2)}
            rows.append(dict(group=group,mu=w,target_length=L,by_length_and_fermion_degree=ds));print(h,group,[sum(d.values()) for d in ds.values()],flush=True)
    args.output.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__': main()
