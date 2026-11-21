"""Apply the free dotted special supercharges to the bare SO(8) level-24 operator."""
import argparse
from collections import Counter
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path

from so8_instanton_exact_modes import run as exact_modes


# Each bit is an occupied SU(4) oscillator. The last three fields carry one
# dotted oscillator as well. The missing-flavour ordering is that of Q psi_i.
FIELDS = (3,5,9,13,11,7)
SIGNS = (1,1,1,1,-1,1)


@lru_cache(None)
def trace(word):
    odd = tuple(int(i.bit_count()==3) for i in word)
    degree = sum(odd)
    normal = {}
    for seq, flags, sign in [(word,odd,1),
            (word[::-1],odd[::-1],(-1)**(len(word)+degree*(degree-1)//2))]:
        for i in range(len(word)):
            parity = sum(flags[:i])
            rotated = seq[i:]+seq[:i]
            normal.setdefault(rotated,set()).add(sign*(-1)**(parity*(degree-parity)))
    first = min(normal)
    return (next(iter(normal[first])) if len(normal[first])==1 else 0), first


def canonical(terms):
    signed = [trace(t) for t in terms]
    sign = 1
    odd = []
    for s, t in signed:
        if not s:
            return 0, ()
        sign *= s
        if sum(v.bit_count()==3 for v in t)%2:
            odd.append(t)
    if len(set(odd)) != len(odd):
        return 0, ()
    inversions = sum(x>y for i,x in enumerate(odd) for y in odd[i+1:])
    return sign*(-1)**inversions, tuple(sorted(t for _,t in signed))


def initial(path):
    inputs = json.loads(gzip.decompress(path.read_bytes()))
    result = Counter()
    for n, multi in inputs:
        sign,t = canonical(tuple(tuple(FIELDS[i] for i in word) for word in multi))
        for word in multi:
            for i in word:
                sign *= SIGNS[i]
        result[t] += n*sign
    return {w:n for w,n in result.items() if n}


def special(source, flavour):
    result = Counter()
    for multi,n in source.items():
        sign = 1
        for i, word in enumerate(multi):
            for j, mask in enumerate(word):
                if mask.bit_count()!=3:
                    continue
                bit = 1<<flavour
                if mask&bit:
                    left = (mask&(bit-1)).bit_count()
                    changed = word[:j]+(mask^bit,)+word[j+1:]
                    phase,target = canonical(multi[:i]+(changed,)+multi[i+1:])
                    if phase:
                        result[target] += n*sign*(-1)**left*phase
                sign = -sign
    return {w:n for w,n in result.items() if n}


def run(path):
    mode = exact_modes(path)
    src = initial(path)
    assert all(sum(i.bit_count()==3 for t in w for i in t)==2 for w in src)
    first = [special(src,i) for i in range(4)]
    for i in range(4):
        assert not special(first[i],i)
    second = {}
    for i,j in combinations(range(4),2):
        row = special(first[i],j)
        reverse = special(first[j],i)
        assert row == {w:-n for w,n in reverse.items()}
        assert all(all(k.bit_count()==2 for t in w for k in t) for w in row)
        second[i,j]=row
    allowed = {():True, **{(i,):bool(first[i]) for i in range(4)},
               **{k:bool(v) for k,v in second.items()}}
    channels = []
    for channel in mode["standard_channels"]:
        retained=[]
        for assignment in channel["assignments"]:
            # Two zeta components are available at the left site for each flavour.
            # Any flavour of exact degree three must use its extra bare mode there.
            mandatory=tuple(i for i,k in enumerate(assignment["left_exact"]) if k==3)
            if allowed.get(mandatory,False):
                retained.append(assignment["left_exact"])
        channels.append(dict(left_extra=channel["left_extra_quartets"],
                             right_extra=channel["right_extra_quartets"],
                             action_quartic=channel["action_quartic_insertions"],
                             previous=len(channel["assignments"]),left_exact=retained))
    return dict(operator_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bare_terms=mode["terms"],canonical_terms=len(src),
                one_action_terms=[len(row) for row in first],
                two_action_terms={f"{i},{j}":len(row) for (i,j),row in second.items()},
                maximum_bare_dotted_modes=2,three_action_terms=0,channels=channels)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operator",type=Path,default=Path(__file__).with_name("so8_level24")/"even_bmn.json.gz")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    record=run(args.operator)
    args.out.write_text(json.dumps(record,indent=2)+"\n")
    print(record["one_action_terms"],record["two_action_terms"])
    print([len(row["left_exact"]) for row in record["channels"]])
