"""Compare the integral Q0 images of *every* certified contraction with the pinned BMN engine."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import scalar_neighbors as production

HERE = Path(__file__).resolve().parent
UPSTREAM = HERE.parent / 'su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/bmn/words.py'
SHA = 'f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e'

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def audit(witness):
    assert digest(UPSTREAM) == SHA
    spec = importlib.util.spec_from_file_location('pinned_bmn', UPSTREAM)
    pinned = importlib.util.module_from_spec(spec); spec.loader.exec_module(pinned)
    cert = json.loads(gzip.decompress(witness.read_bytes()))
    @lru_cache(maxsize=65536)
    def reflect(word):
        # The pinned differential already quotients by *graded* cyclicity. Add
        # transpose conjugation only on its output, where at most one odd survives.
        assert word == pinned.canon_word(word)[1] and pinned.nferm(word) <= 1
        s, rev = pinned.canon_word(word[::-1]); s *= (-1)**len(word)
        if rev == word and s == -1: return 0, ()
        return (s, rev) if rev < word else (1, word)
    def canonical(raw):
        out = defaultdict(int); sign = 1; words = []
        for word in raw:
            s, w = pinned.canon_word(word)
            if not s: return {}
            sign *= s; words.append(w)
        pinned.add_term(out, tuple(words), sign); return dict(out)
    def reference(source):
        out = Counter()
        for multi,c in pinned.Q_expr(source).items():
            words = []; sign = 1
            for word in multi:
                s,w = reflect(word); sign *= s; words.append(w)
                if not sign: break
            if sign:
                s,w = pinned.canon_multi(words); out[w] += c*sign*s
        return {w:c for w,c in out.items() if c}
    def actual(raw):
        out = Counter(); phase = 1
        for i,w in enumerate(raw):
            if any(l >= 3 for l in w):
                other = [production.necklace(t) for j,t in enumerate(raw) if j != i]
                for u,c in production.dtrace(w):
                    k = phase*c; terms = [u]
                    for s,t in other: k *= s; terms.append(t)
                    out[tuple(sorted(terms))] += k
            if pinned.word_parity(w): phase = -phase
        return {w:c for w,c in out.items() if c}
    controls = {'field_strength':((6,0,1,2),), 'two_fermions_in_one_trace':((3,0,4,1),),
                'two_odd_traces':((3,0),(4,1)), 'swapped_odd_traces':((4,1),(3,0)),
                'repeated_odd_trace':((3,0),(3,0))}
    control_counts = {}; results = {}
    for label,raw in controls.items():
        src=canonical(raw); ref=reference(src); assert ref==actual(raw)
        assert not pinned.Q_expr(pinned.Q_expr(src))
        assert bool(ref) == (label != 'repeated_odd_trace')
        results[label]=ref; control_counts[label]=len(ref)
    assert results['two_odd_traces'] == {w:-c for w,c in results['swapped_odd_traces'].items()}
    recipes={production.decode(r) for m in cert['maps'].values() for r in m['source_contractions']}
    counts=Counter(); permap={k:len(m['source_contractions']) for k,m in cert['maps'].items()}
    for kind,layout in sorted(recipes):
        src=Counter()
        for letters,weight in production.assignments(kind):
            raw=tuple(tuple(letters[i] for i in t) for t in layout)
            for w,c in canonical(raw).items(): src[w] += c*weight
        ref=reference({w:c for w,c in src.items() if c}); same=production.expression((kind,layout))
        assert ref == same and ref, (kind,layout)
        if kind=='pp':
            tag='two_fermions_separate' if not any(7 in t and 9 in t for t in layout) else 'two_fermions_together'
            counts[tag]+=1
        else: counts[kind]+=1
    assert all(counts[k] for k in ('two_fermions_separate','two_fermions_together','f','p'))
    return dict(prime_used_for_sign_comparison=None, upstream_commit='4256f207986029540060fb44dd494746828f65ce', upstream_words_sha256=SHA,
                minor_file_sha256=digest(witness), checksum_basis='integral trace expressions before evaluation',
                distinct_recipes=len(recipes),by_source_type=dict(sorted(counts.items())),minor_rows=permap,
                graded_sign_controls=control_counts, pinned_Q0_squared_zero_on_controls=True, unequal_expressions=0)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--witness',type=Path,default=HERE/'scalar_minors.json.gz');p.add_argument('--out',type=Path);args=p.parse_args()
    report=json.dumps(audit(args.witness),indent=2)+'\n'
    if args.out:args.out.write_text(report)
    print(report,end='')
