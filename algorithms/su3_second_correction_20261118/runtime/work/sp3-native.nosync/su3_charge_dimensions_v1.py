"""SU(3) invariant dimensions resolved by charge conjugation.

In the disconnected component, use a real SO(3) torus with fundamental
eigenvalues (z,1,z^-1).  Ad(g) C on sl(3) has eigenvalues
  1,z,z^-1,-1,-z,-z^-1,-z^2,-z^-2.
The normalized twisted Weyl measure is (1-z^4)(1-z^-4)/2.
All coefficient extraction and cycle-index operations are integral.
"""
import sp_environment
from collections import defaultdict
from functools import lru_cache
from itertools import permutations
from math import factorial
from pathlib import Path
import hashlib
import json
import time
from lie_root_dimension_v1 import cycle_data, multiplicity_parts
from screen_gaugino_source_dimensions_v1 import content_patterns


@lru_cache(None)
def twisted_tensor_trace(cycles):
    character = {0: 1}
    eigenvalues = ((1,0),(1,1),(1,-1),(-1,0),(-1,1),(-1,-1),(-1,2),(-1,-2))
    for length in cycles:
        following = defaultdict(int)
        for weight, coefficient in character.items():
            for sign, exponent in eigenvalues:
                following[weight+length*exponent] += coefficient*sign**length
        character = {weight: value for weight,value in following.items() if value}
    numerator = 2*character.get(0,0)-character.get(4,0)-character.get(-4,0)
    assert numerator % 2 == 0
    return numerator//2


@lru_cache(None)
def parity_parts(bosons, fermions):
    terms, denominator = {(): 1}, 1
    for count, odd in [(n,False) for n in bosons if n]+[(n,True) for n in fermions if n]:
        denominator *= factorial(count)
        following = defaultdict(int)
        for previous, coefficient in terms.items():
            for cycle, weight in cycle_data(count,odd):
                following[tuple(sorted(previous+cycle,reverse=True))] += coefficient*weight
        terms = {cycle: weight for cycle,weight in following.items() if weight}
    numerator = sum(weight*twisted_tensor_trace(cycle) for cycle,weight in terms.items())
    assert numerator % denominator == 0
    trace = numerator//denominator
    dimension = multiplicity_parts(bosons,fermions,'SU3')
    assert abs(trace) <= dimension and (dimension-trace)%2 == 0
    return (dimension+trace)//2, (dimension-trace)//2


@lru_cache(None)
def dimensions(mu, length, parity):
    assert parity in (-1,1)
    result = defaultdict(int)
    rho = (2,1,0)
    for spin_shift, spin_sign in ((0,1),(1,-1)):
        for order in permutations(range(3)):
            sign = (-1)**sum(order[i]>order[j] for i in range(3) for j in range(i+1,3))
            shifted = (mu[0]+spin_shift,mu[1]-spin_shift)+tuple(
                mu[i+2]+rho[i]-rho[order[i]] for i in range(3))
            if min(shifted)<0:
                continue
            for bosons,fermions,count in content_patterns(shifted,length):
                result[sum(fermions)] += spin_sign*sign*count*parity_parts(bosons,fermions)[parity==-1]
    assert all(value>=0 for value in result.values()), (mu,length,parity,result)
    return {degree: value for degree,value in sorted(result.items()) if value}


def controls():
    symmetric, exterior = [0]*13, [0]*13
    for quadratic in range(7):
        for cubic in range(5):
            if 2*quadratic+3*cubic <= 12:
                symmetric[2*quadratic+3*cubic] += (-1)**cubic
    for degree, value in ((0,1),(3,1),(5,-1),(8,-1)):
        exterior[degree] = value
    for degree in range(13):
        even,odd = parity_parts((degree,),())
        assert even-odd == symmetric[degree], ('symmetric',degree,even,odd)
        even,odd = parity_parts((),(degree,))
        assert even-odd == exterior[degree], ('exterior',degree,even,odd)
    return dict(one_matrix_symmetric_degrees=13,one_matrix_exterior_degrees=13)


def main():
    import g2_full_classical_screen_v1 as unrefined
    start=time.monotonic()
    rows=[]
    for order in range(5):
        mu=(order,order,7-order,4-order,4-order)
        for role,offset in (('incoming',-1),('correction',0),('target',1),('outgoing',2)):
            length=10-order+offset
            positive=dimensions(mu,length,1)
            negative=dimensions(mu,length,-1)
            whole=unrefined.dimensions(mu,length,'SU3')
            assert {degree:positive.get(degree,0)+negative.get(degree,0)
                    for degree in set(positive)|set(negative)} == whole
            row=dict(order=order,role=role,mu=mu,length=length,
                positive=sum(positive.values()),negative=sum(negative.values()),
                total=sum(whole.values()),positive_by_fermions=positive,negative_by_fermions=negative)
            rows.append(row)
            print(json.dumps(row),flush=True)
    result=dict(status='exact_integer_parity_dimensions',controls=controls(),rows=rows,
        seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Cochain dimensions only; no cohomology or higher differential is inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
