#!/usr/bin/env python3
"""Independent exact Weyl-character count, with no state averaging operator."""
from native import *
from itertools import permutations
from fractions import Fraction

ROOTS=((1,1),(1,-1),(-1,1),(-1,-1),(2,0),(-2,0),(0,2),(0,-2))

@lru_cache(None)
def partitions(n,largest=None):
    if not n:return ((),)
    if largest is None:largest=n
    return tuple((first,)+rest for first in range(min(n,largest),0,-1)
                 for rest in partitions(n-first,first))

@lru_cache(None)
def denominator():
    result={(0,0):1}
    for a,b in ROOTS:
        new=Counter(result)
        for (u,v),n in result.items():new[u+a,v+b]-=n
        result={w:n for w,n in new.items() if n}
    assert result[0,0]==8
    return result

@lru_cache(None)
def invariant_cycles(cycles):
    character={(0,0):1}
    for length in cycles:
        new=Counter()
        for (u,v),n in character.items():
            new[u,v]+=2*n
            for a,b in ROOTS:new[u+length*a,v+length*b]+=n
        character=dict(new)
    value=sum(n*character.get((-a,-b),0) for (a,b),n in denominator().items())
    assert value%8==0
    return value//8

@lru_cache(None)
def invariant_content(bosons,fermions):
    terms={():Fraction(1)}
    for count,odd in [(n,False) for n in bosons]+[(n,True) for n in fermions]:
        new=defaultdict(Fraction)
        for parts in partitions(count):
            centralizer=math.prod(k**v*math.factorial(v) for k,v in Counter(parts).items())
            coefficient=Fraction((-1)**(count-len(parts)) if odd else 1,centralizer)
            for previous,n in terms.items():new[tuple(sorted(previous+parts,reverse=True))]+=n*coefficient
        terms={k:v for k,v in new.items() if v}
    value=sum(n*invariant_cycles(cycles) for cycles,n in terms.items())
    assert value.denominator==1 and value>=0
    return int(value)

def multiplicity(mu,length):
    result=Counter();rho=(2,1,0)
    for shift,spin_sign in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            moved=(mu[0]+shift,mu[1]-shift)+tuple(mu[i+2]+rho[i]-rho[perm[i]] for i in range(3))
            for letters in contents(moved,length):
                counts=Counter(letters)
                bosons=tuple(sorted((v for l,v in counts.items() if not is_ferm(l)),reverse=True))
                fermions=tuple(sorted((v for l,v in counts.items() if is_ferm(l)),reverse=True))
                result[sum(fermions)]+=sign*spin_sign*invariant_content(bosons,fermions)
    assert all(n>=0 for n in result.values())
    return {str(k):v for k,v in sorted(result.items()) if v}

def main():
    declaration=read(HERE/'declaration.json');records=[]
    for block in declaration['native_blocks']:
        length=block['length'];mu=tuple(block['weight_z_z_theta_theta_theta'])
        character=multiplicity(mu,length)
        metric=read(HERE/f'finite_metric_Y{length}.json.gz')
        metric_by_fermions=Counter(str(d) for d,i in metric['independent_basis_labels'])
        assert character==dict(metric_by_fermions),(character,metric_by_fermions)
        records.append({'length':length,'highest_weight_multiplicities_by_fermion_number':character,
                        'finite_highest_weight_dimension':sum(character.values()),
                        'all_flavour_triplet_states':3*sum(character.values())})
        print(records[-1],flush=True)
    result={'passed':True,'group':'USp(4)','adjoint_root_weights':[list(w) for w in ROOTS],
            'adjoint_zero_weight_multiplicity':2,'Weyl_group_order':8,
            'flavour_irreducible_representation_Dynkin_labels':[1,0],
            'dotted_spin':0,'left_spin':1,
            'charge_blocks':records,'total_highest_weight_multiplicity':sum(r['finite_highest_weight_dimension'] for r in records),
            'total_states_in_complete_flavour_irreducible_representations':sum(r['all_flavour_triplet_states'] for r in records),
            'scope':'Exact character coefficients count the complete irreducible-representation multiplicity space; no state averaging or extremal operator is used.'}
    write(HERE/'character_receipt.json',result)

if __name__=='__main__':main()
