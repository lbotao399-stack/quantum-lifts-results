"""Check the full grade-one incoming partner of the level-30 SU(3) F1 seed."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from su3_f1_grade_one_ranks import sectors
from trace_descent import dimensions


def verify():
    directory = HERE/'7659195'
    receipt = json.loads((directory/'receipt.json').read_bytes())
    controls = json.loads((HERE/'grade-one-second-prime.json').read_bytes())
    assert len(controls) == 2
    rows = []
    for i, weight in enumerate(sectors(7)):
        path = directory/f'f7_y8_sector{i}.json'
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == receipt['sha256'][path.name]
        r = json.loads(raw)
        assert r['mode']=='charge_conjugation' and r['gauge']=='SU3' and r['job']=='7659195'
        assert tuple(r['sector'])==weight and r['length']==8 and r['prime']==2147483647
        assert r['seed']==830021+113*8+977*7+1021*i
        charge = (*weight[3:], *weight[:3])
        even, odd = [sum(dimensions(charge,8,3,g).values()) for g in ('SU+','SU-')]
        assert even+odd == r['cochains'] and odd == r['odd_in']+r['odd_out']
        if i >= 3:
            control = controls[i-3]
            assert control['prime'] == 1000003 and control['seed'] == 7310218+i
            for field in ('mode','gauge','sector','length','cochains',
                          'even_in','odd_in','even_out','odd_out'):
                assert control[field] == r[field]
        rows.append((odd,r['odd_in'],r['odd_out']))
    assert rows[0] == (12931,2855,10076)
    return rows


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('run with assertions enabled')
    print(verify())
