"""Check the archived dimension arithmetic and direct highest-weight bounds."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    sources = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    for name, digest in sources.items():
        assert sha(ROOT / name) == digest, name
    hashes = json.loads((ROOT / 'evidence/SHA256SUMS.json').read_text())
    for name, digest in hashes.items():
        assert sha(ROOT / 'evidence' / name) == digest, name
    checked = []
    for group, filename in [('Sp3', 'sp3_parity_index_j30_j32_v1.json'),
                            ('G2', 'g2_parity_index_j28_j30_v1.json')]:
        path = ROOT / 'evidence/inputs' / filename
        data = json.loads(path.read_text())
        columns = next(level['columns'] for level in data['levels'] if level['level'] == 30)
        gate = json.loads((ROOT / 'evidence' / group / 'index_gate.json').read_text())
        assert gate['source_index_sha256'] == sha(path)
        for row in gate['records']:
            column = next(value for value in columns if value['weights'] == row['weights'])
            for sector in column['sectors']:
                assert sector['dimension'] == sum(sector['by_fermions'].values())
                assert sector['dimension'] >= 0
            index = sum((-1) ** ((sum(row['weights']) - sector['sector'][0]) % 2)
                        * sector['dimension'] for sector in column['sectors'])
            assert index == column['even_minus_odd_length_index'] == row['exact_parity_index']
            sectors = []
            for result in sorted((ROOT / 'evidence' / group).glob('s*/result.json')):
                value = json.loads(result.read_text())
                if value['sector'][1:] == row['weights']:
                    assert (sum(row['weights']) - value['sector'][0]) % 2 == 0
                    sectors.append(value)
            lower = sum(value['rank_lower_bound'] for value in sectors)
            assert len(sectors) == row['completed_product_sectors']
            assert lower == row['rational_Cartan_product_lower_bound']
            assert index - lower == row['kernel_parity_index_upper_bound']
            assert max(0, lower - index) == row['forced_odd_length_cohomology_lower_bound']
            checked.append(dict(group=group, **row))
    result = dict(status='passed_archived_dimension_arithmetic_and_direct_rank_bounds',
                  source_hashes=len(sources), payload_hashes=len(hashes), columns=len(checked),
                  negative_upper_bounds=sum(row['kernel_parity_index_upper_bound'] < 0 for row in checked),
                  scope='Checks the archived exact dimensions and certified minor lower bounds; does not recompute character integrals.')
    (ROOT / 'evidence/index_replay.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
