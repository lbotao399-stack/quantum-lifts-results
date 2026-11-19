"""Reconstitute and replay the rational 70-content vector without the production matrices."""
import argparse,hashlib,json
from fractions import Fraction
from math import gcd,lcm
from pathlib import Path
import replay_scalar_layers as r
s=r.s

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def unpack(shared,exact,audit,output):
    q=s.read(exact);z=s.read(shared);test=s.read(audit);D=int(z['denominator']);assert D>0 and D.bit_length()==q['vector_denominator_bits'] and gcd(D,*(int(v) for w,v in z['words']))==1
    words=z['words'];out=[[w,str(Fraction(int(v),D))] for w,v in words];ordered=[s.unpack_hex(x[0]) for x in words];assert ordered==sorted(set(ordered)) and len(out)==test['source_terms']==q['vector_words']==2724;assert digest(out)==q['vector_sha']==test['ordered_source_sha256'];assert q['integer_residual_coordinates']==0 and q['target_coordinates']==5475 and q['determinant_mod_prime']%q['prime'];s.write(output,out);return q

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--directory',type=Path,required=True);a.add_argument('--out-vector',type=Path,required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--exe',type=Path);a.add_argument('--work',type=Path);a.add_argument('--jobs',type=int,default=2);p=a.parse_args();assert __debug__;d=p.directory;q=unpack(d/'vector_shared.json.gz',d/'exact.json',d/'audit.json',p.out_vector)
    if p.exe:
        v=r.replay(p.rhs,p.out_vector,'joint',p.exe,p.work,p.jobs);assert v==s.read(d/'audit.json');print('Identical three-colour replay:',len(v['contents']),sum(x['monomials_collected'] for x in v['contents']),sum(x['nonzero_monomials'] for x in v['contents']),digest(v))
