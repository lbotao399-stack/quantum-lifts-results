"""Exact Sp(rank) full-letter SU(3) and dotted-SU(2) multiplicities.

The charge and letter conventions, and the Weyl-difference formula, are
unchanged from full_dimensions.py. Only the adjoint invariant multiplicity
is replaced by the already validated symplectic character formula.
"""
import sp_environment
from sp_dimensions import multiplicity_parts
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm
from itertools import permutations
from collections import Counter, defaultdict
from functools import lru_cache


@lru_cache(maxsize=None)
def highest_dimensions(sector, letters, rank):
    result=defaultdict(int)
    rho=(2,1,0)
    for derivative_shift, derivative_sign in ((0,1),(2,-1)):
        for permutation in permutations(range(3)):
            sign=(-1)**sum(permutation[i]>permutation[k] for i in range(3) for k in range(i+1,3))
            q=tuple(sector[i+2]+2*(rho[i]-rho[permutation[i]]) for i in range(3))
            shifted=(sector[0]+derivative_shift,sector[1]-derivative_shift)+q
            for content in fw.contents_of_sector(shifted,letters):
                if len(content)!=letters:
                    continue
                counts=Counter(content)
                bosons=tuple(sorted((n for letter,n in counts.items() if not is_ferm(letter)),reverse=True))
                fermions=tuple(sorted((n for letter,n in counts.items() if is_ferm(letter)),reverse=True))
                degree=sum(n for letter,n in counts.items() if is_ferm(letter))
                result[degree]+=sign*derivative_sign*multiplicity_parts(bosons,fermions,rank)
    assert all(value>=0 for value in result.values()),(sector,letters,rank,result)
    return {degree:value for degree,value in result.items() if value}
