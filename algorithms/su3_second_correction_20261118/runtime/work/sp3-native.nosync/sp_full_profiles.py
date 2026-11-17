"""Exact cochain dimensions resolved by field type and total jet order.

The profile is preserved by both highest-weight projectors. It is not
assumed to be preserved by Q0; its couplings must be retained in the map.
"""
import sp_environment
from sp_dimensions import multiplicity_parts
from sp_full_dimensions import highest_dimensions
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import permutations


def profile(content):
    return tuple(sorted(Counter((len(L[2]),L[0]+L[1]) for L in content).items()))


@lru_cache(None)
def dimensions(sector,letters,rank):
    result=defaultdict(int);rho=(2,1,0)
    for ds,ss in ((0,1),(2,-1)):
        for permutation in permutations(range(3)):
            sign=(-1)**sum(permutation[i]>permutation[k] for i in range(3) for k in range(i+1,3))
            q=tuple(sector[i+2]+2*(rho[i]-rho[permutation[i]]) for i in range(3))
            shifted=(sector[0]+ds,sector[1]-ds)+q
            for content in fw.contents_of_sector(shifted,letters):
                if len(content)!=letters:continue
                counts=Counter(content)
                bosons=tuple(sorted((n for L,n in counts.items() if not is_ferm(L)),reverse=True))
                fermions=tuple(sorted((n for L,n in counts.items() if is_ferm(L)),reverse=True))
                degree=sum(n for L,n in counts.items() if is_ferm(L))
                result[degree,profile(content)]+=sign*ss*multiplicity_parts(bosons,fermions,rank)
    assert all(n>=0 for n in result.values())
    result={key:n for key,n in result.items() if n};total=defaultdict(int)
    for (degree,_),n in result.items():total[degree]+=n
    assert dict(total)==highest_dimensions(sector,letters,rank)
    return result
