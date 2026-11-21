"""Encode the pinned parent in its original, unnormalized charge-conjugation sector."""
import argparse
import hashlib
import json
import pickle
import struct
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE = HERE.parent / 'su3_independent_prime_20261118'
RUNTIME = PACKAGE / 'runtime/work/sp3-native.nosync'
P = 1_000_033
PARENT = '63b1f6abf6765cebdb3a27b8fc30d30274f2281272616b193079d943ad729438'
PLAN = '58a440c904dda13ed2f0a34ee3e26543edbd607049dd2774e8af6dadaef65ae6'

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def alphabet():
    from itertools import combinations
    return sorted((a,b,t) for a in range(5) for b in range(5-a) for n in range(4)
                  for t in combinations(range(3),n) if n or a+b)

LETTERS = alphabet()
IDS = {l:i for i,l in enumerate(LETTERS,1)}

def write_alphabet(path):
    assert len(LETTERS)==119
    path.write_text(''.join(f'{a} {b} {sum(1<<i for i in t)}\n' for a,b,t in LETTERS))

def key(w): return b''.join(bytes(IDS[l] for l in tr)+b'\0' for tr in w)

def unpack(k):
    assert k.endswith(b'\0')
    return tuple(tuple(LETTERS[i-1] for i in tr) for tr in k[:-1].split(b'\0'))

def write_records(file, expr):
    with file.open('wb') as f:
        for k,c in sorted(expr.items()):
            if c%P:
                assert 3<=len(k)<=48 and k[-1]==0
                f.write(struct.pack('<BI',len(k),c%P)+k)

def records(file):
    with file.open('rb') as f:
        while head:=f.read(5):
            n,c=struct.unpack('<BI',head);k=f.read(n)
            assert len(k)==n and 3<=n<=48 and 0<c<P and k[-1]==0
            yield k,c

def prepare(folder):
    assert __debug__
    sys.path.insert(0,str(RUNTIME));import su_trace_adapter_v1 as alg
    from fortuity.full.letters import letters_up_to
    assert LETTERS==sorted(l for l in letters_up_to(4,4) if sum(l[:2])<=4)
    plan_file=HERE.parent/'su3_j30_classical_ranks/joint_plan.json'
    parent=ROOT/'data/SU3_j30_parent_p1000033/parent.npy'
    assert sha(plan_file)==PLAN and sha(parent)==PARENT
    plan=json.loads(plan_file.read_bytes());xx=np.load(parent,allow_pickle=False)
    assert xx.shape==(63983,) and xx.dtype==np.dtype('<u4') and 0<xx.min()<=xx.max()<P
    folder.mkdir(parents=True,exist_ok=True);write_alphabet(folder/'alphabet.txt')
    accumulated=defaultdict(int);offset=0;pins=[]
    for i,s in enumerate(plan['sources']):
        f=PACKAGE/'second_basis/shard0'/f'targets{i:02d}/independent_seeds.pkl'
        pin=sha(f);assert i==s['index'] and pin==s['recipe_sha256'];seeds=pickle.loads(f.read_bytes());assert len(seeds)==s['dimension']
        for x,w in zip(xx[offset:offset+len(seeds)],seeds,strict=True):
            for z,v in alg.parity_project({w:1},-1).items():accumulated[key(z)]+=int(x)*v
        offset+=len(seeds);pins.append(pin)
    assert offset==len(xx);out=folder/'parent.input.bin';write_records(out,accumulated)
    receipt=dict(prime=P,parent_sha256=PARENT,plan_sha256=PLAN,recipe_sha256=pins,
                 parent_coefficients=offset,alphabet_sha256=sha(folder/'alphabet.txt'),input_sha256=sha(out),
                 before_highest_weight=sum(1 for _ in records(out)))
    (folder/'input.json').write_text(json.dumps(receipt,indent=2)+'\n');print(receipt['before_highest_weight'],receipt['input_sha256'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path);prepare(p.parse_args().out)
