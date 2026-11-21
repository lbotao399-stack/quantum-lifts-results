"""Check grade-zero and grade-two physical partners of the two SU(3) F1 seeds."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from b2_schur_complex import partner
from rank_two_bmn_product_replay import highest, witness
from rank_two_taylor_products import record, replay
from su3_bmn_parity_products import labelled
from trace_descent import dimensions
from two_flavour_hooks import schur_colour_dimensions
from unitary_hook_parity import resolve, weights


ROOT = Path(__file__).resolve().parent/'taylor_hooks/su3_f1_partners'
PUBLIC = {'bmn_hw_cohomology_j30.json': 'b6ffb2269078652ac417d21e0f3d02c25784cfc4',
          'hw_cohomology_j32_records/su3_c_8_4_4.jsonl': '7405507c4b0f60784863a24da7b971a83b08d181'}


def public_upper(directory):
    inputs = []
    for name, digest in PUBLIC.items():
        raw = (directory/name).read_bytes()
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest() == digest
        inputs.append(raw)
    level30, = [c for c in json.loads(inputs[0])['columns'] if c['lam'] == [7,4,4]]
    assert level30['all_checks'] and all(r['check'] for r in level30['sectors'])
    result = {'30': {15-r['twoJ']: r['h'] for r in level30['sectors'] if r['h']}}
    rows = {r['n']:r for line in inputs[1].splitlines() if line for r in [json.loads(line)]}
    assert set(rows) == set(range(9))
    values = {}
    for n, row in rows.items():
        assert row['check'] and row['N'] == 3 and row['lam'] == [8,4,4]
        assert row['m'] == sum(schur_colour_dimensions(partner((8,4,4),0),16-n,3))
        h = row['m']-row['r_up']-(rows[n-1]['r_up'] if n else 0)
        assert h >= 0
        if h:
            values[16-n] = h
    result['32'] = values
    return json.loads(json.dumps(result))


def checked_products(prefix, expected_weight, expected_length):
    data = json.loads((ROOT/'bmn'/(prefix+'_products.json')).read_bytes())
    assert data['length'] == expected_length
    rows = []
    for recipe in data['inputs']:
        recipe = tuple(map(tuple, recipe))
        assert sum(sum(g[:5]) for g in recipe) == expected_length
        rows.append(highest(labelled(recipe), expected_weight))
    assert len(rows) == data['rational_product_rank'] and json.loads(json.dumps(witness(rows))) == data['minor']
    return len(rows)


def read_receipt(job):
    folder = ROOT/job
    manifest = json.loads((folder/'receipt.json').read_bytes())
    assert manifest['job'] == job
    for name, digest in manifest['files'].items():
        assert name == Path(name).name and hashlib.sha256((folder/name).read_bytes()).hexdigest() == digest
    return set(manifest['files'])


def grade_two(level, tail, length, full_rebuild):
    rows, product_counts = [], []
    native_folder, product_folder = ROOT/'7658864', ROOT/'5692604'
    for i, weight in enumerate(weights(5,tail,True)):
        row = json.loads((native_folder/f'j{level}_sector_{i}.json').read_bytes())
        assert row['sector'] == list(weight) and row['prime'] == 2147483647
        assert row['seed'] == 611209+level*61+i*1009 and row['length'] == length and row['gauge'] == 'SU3'
        for parity, group in [('even','SU+'),('odd','SU-')]:
            dim = sum(dimensions((*weight[3:],*weight[:3]),length,3,group).values())
            row[parity+'_dim'] = dim
            row[parity+'_h'] = dim-row[parity+'_in']-row[parity+'_out']
        assert row['even_dim']+row['odd_dim'] == row['cochains']
        rows.append(dict(row, row_length=5, tail=tail, sector_index=i))
        ranks = []
        for raised in (False,True):
            if raised and i == 4:
                ranks.append(0)
                continue
            case = f'su{level}'
            data = json.loads(gzip.decompress((product_folder/f'{case}_{i}_{int(raised)}.json.gz').read_bytes()))
            expected = (*weight[:3],weight[3]+int(raised),weight[4]-int(raised))
            assert (data['sector'],data['index'],data['length'],data['case'],data['prime'],data.get('raised_spin',False)) == (
                list(expected),i,length,case,2147483647,raised)
            assert replay(data) == data['det']
            if full_rebuild:
                assert json.loads(json.dumps(record(case,i,raised,length=length))) == data
            ranks.append(data['rank'])
        product_counts.append(ranks[0]-ranks[1])
        assert product_counts[-1] == row['even_h']+row['odd_h']
    result = resolve(rows)
    lower = sum((-1)**i*r for i,r in enumerate(product_counts))
    parity = 'odd' if level == 30 else 'even'
    assert result['parities'][parity]['h'] == lower and sum(v['h'] for v in result['parities'].values()) == lower
    return dict(length=length, partition=result['partition'], parities=result['parities'],
                protected_product_multiplicities=product_counts, rational_product_rank=lower)


def calculate(upstream, full_rebuild):
    assert len(read_receipt('7658864')) == 10 and len(read_receipt('5692604')) == 18
    upper = json.loads((ROOT/'bmn/public_upper.json').read_bytes())
    if upstream:
        assert public_upper(upstream) == upper
    checked = {}
    for level, first in [(30,7),(32,8)]:
        shape = (first,4,4)
        y = first+3
        bmn = {which: checked_products(f'F1_{first}_{which}',shape,length)
               for which, length in [('incoming',y-1),('outgoing',y+1)]}
        ub = {int(k):v for k,v in upper[str(level)].items()}
        assert bmn['outgoing'] == ub[y+1]
        if level == 30:
            assert bmn['incoming'] == ub[y-1]
        else:
            extra = {n:checked_products(f'F1_8_length{n}',shape,n) for n in (8,9,11)}
            assert extra == {8:1,9:3,11:7} and ub == {8:1,9:3,10:8,11:8,12:1}
            even_index = sum((-1)**n*schur_colour_dimensions(partner(shape,0),n,3)[0]
                             for n in range(max(shape),sum(shape)+1))
            assert even_index == 8
            assert even_index-extra[8]-bmn['outgoing']+(ub[9]-extra[9])+(ub[11]-extra[11]) == bmn['incoming']
        grade1 = {name:schur_colour_dimensions(partner(shape,1),n,3)[level==30]
                  for name,n in [('incoming_cochains',y-2),('outgoing_cochains',y)]}
        assert partner(shape,3) is None
        checked[str(level)] = dict(seed_length=y, charge_conjugation=(-1)**(y+1),
            BMN_product_dimensions=bmn, grade_one_selected_C=grade1,
            grade_two_incoming=grade_two(level,first-4,y-3,full_rebuild))
    return checked


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream',type=Path,help='Pinned fortuity-ai/data/su3 directory')
    parser.add_argument('--rebuild-products',action='store_true')
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error('Run with character and rank checks enabled')
    args.out.write_text(json.dumps(calculate(args.upstream,args.rebuild_products),indent=2)+'\n')
