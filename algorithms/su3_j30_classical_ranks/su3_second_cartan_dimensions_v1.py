"""Exact Weyl-invariant Cartan cochain dimensions for the SU3 second target.

The SU3 Cartan representation of S3 has eigenvalues (1,1), (1,-1), and
(omega,omega^2) in the three conjugacy classes of sizes 1,3,2. Coefficient
sequences below are the exact symmetric/exterior powers of that representation.
These count Cartan cochains, not closed lifts or quantum obstructions.
"""
import sp_environment
from collections import Counter,defaultdict
from itertools import permutations
from functools import lru_cache
from math import comb
from pathlib import Path
import hashlib,json,time
from g2_full_classical_screen_v1 import contents
from fortuity.full.letters import is_ferm
from sp_full_profiles import profile

HERE=Path(__file__).resolve().parent


def factor(multiplicity,odd,kind):
    if odd:
        if multiplicity>2:return 0
        return ((1,2,1),(1,0,-1),(1,-1,1))[kind][multiplicity]
    if kind==0:return multiplicity+1
    if kind==1:return int(multiplicity%2==0)
    return (1,-1,0)[multiplicity%3]


@lru_cache(maxsize=None)
def invariant(bosons,fermions):
    total=0
    for kind,size in enumerate((1,3,2)):
        coefficient=1
        for multiplicity in bosons:coefficient*=factor(multiplicity,False,kind)
        for multiplicity in fermions:coefficient*=factor(multiplicity,True,kind)
        total+=size*coefficient
    assert total>=0 and total%6==0
    return total//6


def dimensions(mu,length):
    answer=defaultdict(int);rho=(2,1,0)
    for spin_shift,spin_sign in ((0,1),(1,-1)):
        for order in permutations(range(3)):
            sign=(-1)**sum(order[i]>order[j] for i in range(3) for j in range(i+1,3))
            shifted=(mu[0]+spin_shift,mu[1]-spin_shift)+tuple(mu[i+2]+rho[i]-rho[order[i]] for i in range(3))
            if min(shifted)<0:continue
            for content in contents(shifted,length):
                counts=Counter(content)
                bosons=tuple(sorted((n for L,n in counts.items() if not is_ferm(L)),reverse=True))
                fermions=tuple(sorted((n for L,n in counts.items() if is_ferm(L)),reverse=True))
                answer[sum(fermions),profile(content)]+=spin_sign*sign*invariant(bosons,fermions)
    assert all(n>=0 for n in answer.values())
    return {key:n for key,n in answer.items() if n}


def main():
    start=time.monotonic()
    for degree in range(15):
        expected=sum(2*quadratic+3*cubic==degree for quadratic in range(8) for cubic in range(6))
        assert invariant((degree,),())==expected
        assert invariant((),(degree,))==int(degree==0)
    mu=(2,2,5,2,2);length=9
    profiles=dimensions(mu,length);by_degree=defaultdict(int)
    for (degree,_),number in profiles.items():by_degree[degree]+=number
    rows=[dict(fermions=degree,profile=p,dimension=n) for (degree,p),n in sorted(profiles.items())]
    result=dict(status='exact_SU3_Cartan_Weyl_invariant_cochain_dimensions',mu=mu,length=length,
                dimension=sum(profiles.values()),by_fermion_degree=dict(by_degree),profiles=rows,
                one_species_controls=30,seconds=time.monotonic()-start,
                current_probe_coordinates_by_degree={0:2,2:56,4:140,6:56,8:2},
                closed_lift_surjectivity_proved=False,quantum_obstruction_computed=False,
                driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'su3_second_cartan_dimensions_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='profiles'}),flush=True)


if __name__=='__main__':main()
