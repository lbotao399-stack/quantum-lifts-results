"""Check the production generators on every letter, in the free associative algebra over Z."""
import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path
from inputs import LETTERS, IDS, P, RUNTIME, HERE, sha


def run(native, alphabet):
    assert __debug__
    sys.path.insert(0, str(RUNTIME)); import sp_environment
    from fortuity.full import letters
    qpin = '0b3561e55f22d5135f2a04ba2996f1652e92af135556a7388fd423974072c8c5'
    assert sha(letters.__file__) == qpin and sha(alphabet) == '8671cea568fe94fa8c0675b0fc5006a3b2f3394f8307aa780195829f89e05356'
    q = {IDS[l]:[(c, (IDS[x], IDS[y])) for c,x,y in letters.Q_letter(l)] for l in LETTERS}
    odd = {IDS[l]:letters.is_ferm(l) for l in LETTERS}
    tbl = subprocess.check_output([native,alphabet]); rows = [list(map(int,line.split())) for line in tbl.splitlines()]
    assert len(rows)==119 and all(len(a)==28 and a[0]==i for i,a in enumerate(rows,1))
    names = ('S0','S1','S2','Tp0','Tp1','E12','F21','E02','F20','E01','F10','spin_raise','spin_lower','spin_exchange')
    checks={}
    for k,name in enumerate(names):
        swap = name=='spin_exchange'; degree = int(k<5)
        op={r[0]:(1,r[-1]) if swap else (r[2+2*k],r[1+2*k]) for r in rows}
        for i,(c,j) in op.items():
            assert (j in odd and (odd[i]^odd[j])==degree) if c else j==0
        left_nonzero=right_nonzero=compared=integer_max=0
        for i in q:
            ca,u=op[i]; left=Counter();right=Counter()
            if ca:
                for b,w in q[u]:left[w]+=ca*b
            for b,(x,y) in q[i]:
                cx,vx=op[x];cy,vy=op[y]
                if swap:right[vx,vy]+=b
                else:
                    if cx:right[vx,y]+=b*cx
                    if cy:right[x,vy]+=(-1 if degree and odd[x] else 1)*b*cy
            want={w:(-b if degree else b) for w,b in right.items() if b};got={w:b for w,b in left.items() if b}
            assert got==want, (name,i,got,want)
            left_nonzero+=bool(got);right_nonzero+=bool(want);compared+=len(got);integer_max=max(integer_max,*map(abs,got.values())) if got else integer_max
        checks[name]=dict(letters=len(q),nonzero_comparisons=left_nonzero,nonzero_ordered_word_coefficients=compared,max_abs_coefficient=integer_max)
    return dict(field_of_letter_check='integers',parent_field=P,alphabet=119,ordered_product=True,letter_differential_sha256=qpin,
        table_sha256=hashlib.sha256(tbl).hexdigest(),lowering_sha256=sha(HERE/'lower.cpp'),native_projector_sha256=sha(RUNTIME/'su_full_projection_native_v1.cpp'),
        all_letter_identities=sum(v['letters'] for v in checks.values()),checks=checks,
        parent_verifier='algorithms/su3_j30_classical_ranks/verify_parent.py',
        conclusion='Both composite lowerings anticommute with Q0 on traces; the saved modular descendants are closed at the parent prime. No boundary verdict follows.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('binary',type=Path);p.add_argument('alphabet',type=Path);p.add_argument('--output',type=Path);a=p.parse_args();r=run(a.binary.resolve(),a.alphabet.resolve());txt=json.dumps(r,indent=2)+'\n';
    if a.output:a.output.write_text(txt)
    print(txt)
