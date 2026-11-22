"""Check the highest components of the six finite-colour relations and the R20 prefactors."""
import argparse
from collections import defaultdict
import json
from math import lcm
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from cartan_multiplier import add, mul, prefactors, trace
from syzygy import TERMS, load_vendor, preimage20, relation20

PARENT = Path(__file__).resolve().parents[1]


def restrict(expression):
    masks = (1,2,4,6,5,3,7)
    terms = []
    for word,c in expression.items():
        row = {():1}
        for tr in word:
            row = mul(row,trace(tuple(masks[i] for i in tr)))
            c *= (-1)**(tr.count(4)+tr.count(6))
        terms.append((c,row))
    return add(*terms)


def check_prefactors(dmk):
    pieces = [[] for _ in range(3)]
    for _,coefficient,upper,lower in TERMS[-3:]:
        for k in range(3):
            hi,lo = [0]*5,[2]*3
            hi[int(upper)] = lo[int(lower)] = k
            row = dmk.mul(dmk.v2(hi[1],lo[1]),dmk.v3(hi[2],hi[3],lo[2]))
            pieces[hi[4]].append((coefficient,row))
    assert [restrict(dmk.add(*part)) for part in pieces] == prefactors()


def colour_expand(exe,folder,words,polynomial):
    groups = defaultdict(dict)
    for multi,c in polynomial.items():
        groups[tuple(sorted(i for tr in multi for i in tr))][multi] = c
    assert groups
    distinct = nonzero = 0
    for content,rows in sorted(groups.items()):
        scalar = max((i for i in set(content) if not words.is_ferm(i)),key=content.count)
        den = lcm(*(v.denominator for v in rows.values()))
        lines = [f'{scalar} 030405 {len(rows)}']
        for multi,c in rows.items():
            encoded = ' '.join(bytes(tr).hex() for tr in multi)
            lines.append(f'{int(c*den)} {len(multi)} {encoded}')
        input_file,output_file = folder/'words.txt', folder/'result.json'
        input_file.write_text('\n'.join(lines)+'\n')
        subprocess.run([str(exe),str(input_file),str(output_file),str(len(content))],check=True)
        packet = json.loads(output_file.read_bytes())
        distinct += packet['distinct_monomials']
        nonzero += packet['nonzero_monomials']
    return dict(contents=len(groups),formal_words=len(polynomial),distinct=distinct,nonzero=nonzero)


def verify(vendor):
    dmk,words = load_vendor(vendor)
    check_prefactors(dmk)
    highest = [
        ('12[0,3]',dmk.tensor('r12_03')(2,2,2),dmk.tensor('R12_03')(2,2,2)),
        ('12[1,1]',dmk.tensor('r12_11')(0,2),dmk.tensor('R12_11')(0,2)),
        ('12[2,2]',dmk.tensor('r12_22')(0,0,2,2),dmk.tensor('R12_22')(0,0,2,2)),
        ('14',dmk.tensor('r14_13')(0,2,2,2),dmk.tensor('R14_13')(0,2,2,2)),
        ('16',dmk.tensor('r16_23')(0,0,2,2,2),dmk.tensor('R16_23')(0,0,2,2,2)),
        ('20',preimage20(dmk,0),relation20(dmk,0)),
    ]
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        exe = folder/'colour'
        subprocess.run(['c++','-O3','-std=c++17',str(PARENT/'source_closure.cpp'),'-o',str(exe)],check=True)
        result = {}
        for label,primitive,relation in highest:
            residual = dmk.add((1,words.Q_expr(primitive)),(-1,relation))
            data = colour_expand(exe,folder,words,residual)
            control = colour_expand(exe,folder,words,relation)
            assert data['nonzero'] == 0 and control['nonzero'] > 0,(label,data,control)
            result[label] = dict(residual=data,positive_control=control)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vendor',type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('enable assertions')
    print(json.dumps(verify(args.vendor),indent=2))
