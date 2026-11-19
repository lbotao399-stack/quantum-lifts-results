"""Reconstruct every parent of the saved outgoing words, then evaluate with upstream Q_multi.

No imports from the witness producer. Only selected parents are kept from the
full rational input stream. The recorded contributor lists are not read.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip, hashlib, json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
VENDOR = HERE.parent/'su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor'
sys.path.insert(0, str(VENDOR))
from fortuity.full.letters import letters_up_to
from fortuity.full.words import canon_word, canon_multi, Q_multi

ALPH = tuple(sorted(letters_up_to(1, 1)))
LETTER = {v: i for i, v in enumerate(ALPH)}
INPUT = '9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'

def read(path):
    with gzip.open(path, 'rb') as stream: return json.load(stream)

def write(path, data):
    raw = (json.dumps(data, sort_keys=True, separators=(',', ':'))+'\n').encode()
    path.write_bytes(gzip.compress(raw, mtime=0) if path.suffix == '.gz' else raw)

def packed(word): return tuple(bytes(LETTER[x] for x in t) for t in word)
def expanded(word): return tuple(tuple(ALPH[i] for i in t) for t in word)
def unpack_hex(word): return tuple(bytes.fromhex(t) for t in word)
def hexes(word): return [t.hex() for t in word]

@lru_cache(maxsize=16384)
def trace_degrees(t):
    ls = [ALPH[x] for x in t]
    return (sum(x[0] for x in ls), sum(x[1] for x in ls), *(sum(i in x[2] for x in ls) for i in range(3)))

def block(word): return tuple(sorted(trace_degrees(t) for t in word))

def parents(word):
    # Reverse Taylor-label addition in Psi^2, not the producer's inverse split map.
    found = set()
    for k, trace in enumerate(word):
        if len(trace) < 3: continue  # its putative predecessor is a traceless one-letter trace
        for i in range(len(trace)):
            rotated = trace[i:] + trace[:i]
            (z,p,s), (zz,pp,ss) = rotated[:2]
            if set(s).intersection(ss): continue
            merged = (z+zz, p+pp, tuple(sorted(s+ss)))
            if merged not in LETTER: continue
            sig, parent = canon_word((merged,)+rotated[2:])
            if sig:
                sig, candidate = canon_multi(word[:k]+(parent,)+word[k+1:])
                if sig: found.add(packed(candidate))
    return found

def stream_rows(path, digest):
    with gzip.open(path, 'rb') as f:
        for line in f:
            digest.update(line)
            coefficient, raw = json.loads(line)
            word = tuple(bytes(LETTER[z,p,tuple(s)] for z,p,s in tr) for tr in raw)
            yield word, coefficient

def run(rhs, witnesses):
    saved = read(witnesses)
    assert saved['letter_table'] == [[z,p,list(s)] for z,p,s in ALPH]
    targets = [unpack_hex(b['witness']['outgoing']) for b in saved['blocks']]
    assert len(targets) == len(set(targets)) == 733
    reverse = defaultdict(list)
    for i, w in enumerate(targets):
        for u in parents(expanded(w)):
            assert block(w) == block(u)
            reverse[u].append(i)
    contributions = [[] for _ in targets]; seen = set(); block_sizes = Counter(); digest = hashlib.sha256()
    contents = Counter()
    for u, coefficient in stream_rows(rhs, digest):
        block_sizes[block(u)] += 1; contents[bytes(sorted(b''.join(u)))] += 1
        if u not in reverse: continue
        assert u not in seen; seen.add(u)
        image = {packed(v): n for v,n in Q_multi(expanded(u)).items()}
        for i in reverse[u]:
            n = image.get(targets[i],0)
            if n: contributions[i].append((hexes(u), coefficient, n))
    assert digest.hexdigest() == saved['input_stream_sha256'] == INPUT
    assert sum(block_sizes.values()) == 517468 and len(block_sizes) == 733
    rows = []
    for b,w,terms in zip(saved['blocks'], targets, contributions):
        exact = sum((Fraction(c)*n for _,c,n in terms), Fraction())
        assert block(w) == tuple(map(tuple,b['degrees'])) and block_sizes[block(w)] == b['rhs_terms']
        assert exact and exact == Fraction(b['witness']['coefficient'])
        rows.append(dict(outgoing=hexes(w),coefficient=str(exact),parents=sorted(terms)))
    return dict(rhs_sha256=digest.hexdigest(),rhs_terms=sum(block_sizes.values()),replayed=len(rows),
                possible_parents=len(reverse),present_parents=len(seen),nonzero_parent_incidences=sum(map(len, contributions)),
                target_exact_letter_contents={k.hex():v for k,v in sorted(contents.items())}, coefficients=rows)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rhs',type=Path,default=HERE/'negative_Q1.jsonl.gz');p.add_argument('--witnesses',type=Path,default=HERE/'finite_colour/formal_blocks.json.gz');p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); assert __debug__; result=run(a.rhs,a.witnesses);write(a.out,result);print({k:v for k,v in result.items() if k not in ('target_exact_letter_contents','coefficients')})
