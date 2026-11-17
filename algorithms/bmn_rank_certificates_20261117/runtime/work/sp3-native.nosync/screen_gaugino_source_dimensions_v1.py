"""Exact low-cost source-space screen for the G2 S^3 target and controls.

This computes cochain dimensions, not differential ranks. The three source
gradings follow Q_n: mu shifts by (n,n,-n,-n,-n) and Y by 1-n.
"""
import sp_environment
from lie_root_dimension_v1 import multiplicity_parts
from sp_dimensions import multiplicity_parts as sp_multiplicity
from itertools import product, permutations
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import signal
import time

HERE = Path(__file__).resolve().parent


@lru_cache(None)
def content_patterns(mu, letters):
    if min(mu) < 0 or letters < 0 or sum(mu) < letters:
        return ()
    # Since every admissible letter has positive total degree, no single
    # letter can exceed this bound in a word with the specified length.
    maximum = sum(mu)-letters+1
    alphabet = []
    for a,b in product(range(mu[0]+1), range(mu[1]+1)):
        for bits in product((0,1), repeat=3):
            weight = (a,b)+bits
            degree = sum(weight)
            if not 1 <= degree <= maximum or any(weight[i] > mu[i] for i in range(5)):
                continue
            alphabet.append((weight, (sum(bits)+1)%2))
    alphabet.sort(reverse=True)
    counts = defaultdict(int)
    def visit(start, remaining, left, bosons, fermions):
        if not left:
            if not any(remaining):
                counts[(tuple(sorted(bosons,reverse=True)),tuple(sorted(fermions,reverse=True)))] += 1
            return
        if sum(remaining) < left:
            return
        for index in range(start, len(alphabet)):
            weight, odd = alphabet[index]
            maximum_count = min([left]+[remaining[i]//weight[i] for i in range(5) if weight[i]])
            for count in range(1, maximum_count+1):
                new = tuple(remaining[i]-count*weight[i] for i in range(5))
                visit(index+1,new,left-count,bosons+(() if odd else (count,)),fermions+((count,) if odd else ()))
    visit(0,mu,letters,(),())
    return tuple((bos,fer,multiplicity) for (bos,fer),multiplicity in sorted(counts.items()))


def highest_dimension(mu, letters, group):
    total = 0
    rho = (2,1,0)
    for spin_shift, spin_sign in ((0,1),(1,-1)):
        for permutation in permutations(range(3)):
            sign = (-1)**sum(permutation[i]>permutation[k] for i in range(3) for k in range(i+1,3))
            shifted = (mu[0]+spin_shift,mu[1]-spin_shift)+tuple(mu[i+2]+rho[i]-rho[permutation[i]] for i in range(3))
            for bosons,fermions,count in content_patterns(shifted,letters):
                total += spin_sign*sign*count*multiplicity_parts(bosons,fermions,group)
    assert total >= 0,(mu,letters,group,total)
    return total


def controls():
    records = []
    for group, scalar_degrees, exterior_degrees in (('G2',(2,6),(3,11)),('SU3',(2,3),(3,5)),('Sp2',(2,4),(3,7))):
        scalar = [1]+[0]*9
        exterior = [1]+[0]*9
        for degree in scalar_degrees:
            for n in range(degree,10):scalar[n] += scalar[n-degree]
        for degree in exterior_degrees:
            for n in range(9,degree-1,-1):exterior[n] += exterior[n-degree]
        for n in range(10):
            a=multiplicity_parts((n,),(),group);b=multiplicity_parts((),(n,),group)
            assert (a,b)==(scalar[n],exterior[n]),(group,n,a,b)
            records.append(dict(group=group,degree=n,symmetric=a,exterior=b))
    mixed = 0
    for bosons,fermions in (((2,1),(1,)),((1,1),(2,)),((2,2),(1,1)),((1,1,1),(3,)),((3,),(2,2)),((2,),(3,3))):
        assert multiplicity_parts(bosons,fermions,'Sp2')==sp_multiplicity(bosons,fermions,2)
        mixed += 1
    return dict(symmetric_exterior_records=records,independent_mixed_Sp2_comparisons=mixed)


def main():
    start = time.monotonic()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('dimension screen cap')))
    signal.alarm(180)
    checked = controls()
    records = []
    for group in ('G2','SU3','Sp2'):
        for order in (1,2,3):
            mu=(3-order,3-order,order,order,order)
            source_letters=5+order
            row=dict(group=group,supercharge_order=order,source_mu=mu,source_letters=source_letters,
                source_highest_weight_cochains=highest_dimension(mu,source_letters,group),
                incoming_highest_weight_cochains=highest_dimension(mu,source_letters-1,group),
                outgoing_highest_weight_cochains=highest_dimension(mu,source_letters+1,group))
            records.append(row);print(json.dumps(row),flush=True)
    signal.alarm(0)
    result=dict(status='exact_integer_source_dimensions_only',controls=checked,records=records,seconds=time.monotonic()-start,
        sources={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                 ('lie_root_dimension_v1.py','screen_gaugino_source_dimensions_v1.py','sp_dimensions.py')},
        scope='No differential was evaluated. Positive cochain dimensions do not imply cohomology or any lift. The S^3 target is non-Coulomb in G2; its SU3 and Sp2 counterparts are classical-ring controls.')
    (HERE/'gaugino_source_dimension_screen_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('records','controls')}),flush=True)


if __name__=='__main__':main()
