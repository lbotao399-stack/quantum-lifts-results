"""Verify one invariant first correction over the square-zero scalar scheme at level 32."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from math import lcm, prod
from pathlib import Path

import su3_square_zero_primitive as z
from su3_square_zero_image import expand_content


def parity(trace):
    return sum(z.ODD[k] for k in trace) % 2


def normalize(traces):
    factors = [z.canonical(t) for t in traces]
    sign, word = z.multi(tuple(t for _,t in factors))
    return sign * prod(n for n,_ in factors), word


def rotate(trace, position):
    before, after = trace[:position], trace[position:]
    return (-1)**(parity(before)*parity(after)), after+before


def gap_words(word, counts):
    """Expand a product of traces in words with a scalar--gaugino--scalar gap."""
    i, left, _, left_sign, single = z.commuting_gap(word)
    if not single:
        return {word:1}
    others = [(j,k,v) for j,t in enumerate(word) if i!=j
              for k,v in enumerate(t) if v in z.SCALARS]
    occurrences = Counter(v for _,_,v in others)
    same = left[0] in occurrences
    j,k,_ = next(item for item in others
                 if (item[2]==left[0] if same else occurrences[item[2]]>=2))
    right_sign, right = rotate(word[j], k)
    order = [i,j]+[h for h in range(len(word)) if h not in (i,j)]
    odd = [h for h in order if parity(word[h])]
    sign = left_sign*right_sign*(-1)**sum(h>v for l,h in enumerate(odd) for v in odd[l+1:])
    rest = tuple(word[h] for h in order[2:])
    original_sign, original = normalize((left,right)+rest)
    assert original==word and sign*original_sign==1
    ans = Counter()

    def insert(traces, scale):
        phase, term = normalize(traces+rest)
        if phase and not any(map(z.adjacent_scalars,term)):
            ans[term] += sign*scale*phase

    insert((left+right,),1)
    if same:
        counts["equal_scalar_splices"] += 1
    else:
        counts["polarized_splices"] += 1
        swapped_left, swapped_right = right[:1]+left[1:], left[:1]+right[1:]
        insert((swapped_left+swapped_right,),1)
        phase, swapped = normalize((swapped_left,swapped_right)+rest)
        if phase and not any(map(z.adjacent_scalars,swapped)):
            for term,n in gap_words(swapped, counts).items():
                ans[term] -= sign*phase*n
    ans = {w:c for w,c in ans.items() if c}
    assert all(not z.commuting_gap(w)[4] for w in ans)
    return ans


def polynomial_residual(difference, transpose):
    by_content = defaultdict(list)
    for word,c in difference.items():
        if transpose:
            sign = (-1)**sum(parity(t[:i])*int(z.ODD[k]) for t in word for i,k in enumerate(t))
            phase, word = normalize(tuple(t[::-1] for t in word))
            c *= sign*phase
        by_content[bytes(sorted(b"".join(word)))].append((word,c))
    return sum(len(expand_content(rows)) for rows in by_content.values())


def run(output, build):
    source = z.terms()
    counts = Counter()
    changed = Counter()
    for w,c in source.items():
        for u,n in gap_words(w,counts).items():
            changed[u] += c*n
    changed = {w:c for w,c in changed.items() if c}
    difference = Counter(source)
    difference.subtract(changed)
    difference = {w:c for w,c in difference.items() if c}
    checks = [polynomial_residual(difference, rev) for rev in (False,True)]
    assert checks == [0,0]
    if build:
        z.save(output,z.build(changed))
    archive = json.loads(gzip.decompress(output.read_bytes()))
    assert all(k<len(z.ALPH) for w,c in archive for t in w for k in bytes.fromhex(t))
    verification = z.replay(changed,output)
    assert not verification["residual_terms"]
    return dict(level=32, relation="all nine ordered matrix products of undifferentiated scalars vanish",
                first_image_words=verification["full_first_image_terms"], words_without_adjacent_scalars=len(source),
                singleton_gap_inputs=sum(z.commuting_gap(w)[4] for w in source), splices=dict(counts),
                rewritten_words=len(changed), trace_difference_terms=len(difference),
                independent_component_residuals=dict(common_image=checks[0],common_kernel=checks[1]),
                invariant_correction_terms=len(archive), common_denominator=lcm(*(Fraction(c).denominator for _,c in archive)),
                correction_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                ordinary_trace_variation_terms=verification["Q0_terms_after_scalar_and_projector_relations"],
                residual_after_trace_rewrite=verification["residual_terms"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--correction",type=Path,required=True)
    p.add_argument("--build",action="store_true")
    p.add_argument("--out",type=Path,required=True)
    a = p.parse_args()
    if not __debug__:
        p.error("enable assertions for replay")
    result = run(a.correction,a.build)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)
