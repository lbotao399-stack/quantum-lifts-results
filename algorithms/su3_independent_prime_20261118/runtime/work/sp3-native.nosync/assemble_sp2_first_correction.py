"""Audit common frames and all structural edges before coupling the equations."""
import argparse
import json
from datetime import datetime, timezone

import numpy as np
from correction_profile_inventory import inventory
from integer_j30_boundary_system import HERE, PRIME, digest
from projected_checkpoint import atomic, json_bytes
from sp2_q1_profile_evaluation import target_configuration, save_array


def main(write_matrix=False):
    graph = json.loads((HERE / 'first_correction_profile_graph.json').read_text())
    source_inventory = inventory()
    sources = {(r['degree'], r['number']):dict(r) for r in source_inventory['records']}
    degree7 = HERE / 'sp2_j30_first_correction_map7'
    certificate7 = json.loads((degree7 / 'rational_rank.json').read_text())
    assert certificate7['status'] == 'passed'
    assert certificate7['rational_source_sha256'] == digest(degree7 / 'rational_source7.pkl')
    sources[7, 0] = dict(degree=7, number=0, dimension=552,
        folder=degree7.name, basis_file='rational_source7.pkl',
        rational_basis_sha256=certificate7['rational_source_sha256'],
        certificate_sha256=digest(degree7 / 'rational_rank.json'))
    source_records, target_records, edges = [], [], []
    offset = 0
    for r in graph['source_profiles']:
        key = r['degree'], r['number']
        source_records.append(dict(**r, offset=offset, certificate=sources.get(key)))
        offset += r['dimension']
    assert offset == 30148
    columns = 0
    for target in graph['target_profiles']:
        key = target['degree'], target['number']
        config = target_configuration(*key)
        folder = HERE / 'sp2_j30_first_correction_rhs' / ('target%d_%d' % key)
        result_path = folder / 'result.json'
        record = dict(**target, offset=columns, columns=config['columns'],
                      homogeneous=None, complete=False)
        if result_path.exists():
            result = json.loads(result_path.read_text())
            manifest = json.loads((folder / 'manifest.json').read_text())
            assert all(manifest[k] == v for k,v in config.items())
            assert result['manifest_sha256'] == digest(folder / 'manifest.json')
            assert result['vector_sha256'] == digest(folder / 'vector.npy')
            values = np.load(folder / 'vector.npy', allow_pickle=False)
            assert values.shape == (config['columns'],)
            record.update(complete=True, homogeneous=result['terms'] == 0,
                          terms=result['terms'], vector_sha256=result['vector_sha256'],
                          manifest_sha256=result['manifest_sha256'])
        target_records.append(record)
        columns += config['columns']
    for record in graph['records']:
        source = tuple(record['source'][:2])
        for target_degree, target_number, target_dimension in record['targets']:
            config = target_configuration(target_degree, target_number)
            folder = HERE / 'sp2_j30_first_correction_blocks' / (
                'source%d_%d_target%d_%d' % (*source, target_degree, target_number))
            edge = dict(source=list(source), target=[target_degree, target_number],
                        source_dimension=record['source'][2], target_dimension=target_dimension,
                        complete=False, folder=folder.name)
            if (folder / 'result.json').exists():
                result = json.loads((folder / 'result.json').read_text())
                manifest = json.loads((folder / 'manifest.json').read_text())
                assert result['status'] == 'complete_evaluation'
                assert all(manifest[k] == v for k,v in config.items())
                assert manifest['rational_source_sha256'] == sources[source]['rational_basis_sha256']
                assert result['matrix_sha256'] == digest(folder / 'matrix.npy')
                assert result['manifest_sha256'] == digest(folder / 'manifest.json')
                values = np.load(folder / 'matrix.npy', mmap_mode='r', allow_pickle=False)
                assert values.shape == (record['source'][2], config['columns'])
                edge.update(complete=True, matrix_sha256=result['matrix_sha256'],
                            manifest_sha256=result['manifest_sha256'],
                            rank_lower_bound=result['rank_lower_bound'])
            edges.append(edge)
    assert len(edges) == 53 and len(target_records) == 13 and len(source_records) == 17
    complete = all(r['certificate'] for r in source_records) and all(r['complete'] for r in target_records + edges)
    result = dict(status='complete_inputs' if complete else 'incomplete_inputs',
                  updated_at=datetime.now(timezone.utc).isoformat(), prime=PRIME,
                  rational_source_dimension=sum(r['dimension'] for r in source_records if r['certificate']),
                  full_source_dimension=30148, full_target_dimension=33685,
                  evaluation_columns=columns, complete_map_blocks=sum(r['complete'] for r in edges),
                  required_map_blocks=len(edges), complete_rhs_profiles=sum(r['complete'] for r in target_records),
                  required_rhs_profiles=13, structural_zero_blocks=17*13-53,
                  homogeneous_target_dimension=5160, homogeneous_profiles=7,
                  sources=source_records, targets=target_records, edges=edges,
                  driver_sha256=digest(__import__('pathlib').Path(__file__)),
                  full_obstruction_computed=False,
                  scope='Audit of every common-frame block and source ordering. '
                        'Missing maps remain missing; homogeneous profiles are never discarded.')
    out = HERE / 'sp2_j30_first_correction_assembly'
    out.mkdir(exist_ok=True)
    if write_matrix:
        assert complete, 'Complete every source, map block and RHS profile first'
        matrix = np.lib.format.open_memmap(out / 'matrix.npy', mode='w+', dtype=np.int32,
                                          shape=(30148, columns))
        matrix[:] = 0
        source_offsets = {(r['degree'],r['number']):r['offset'] for r in source_records}
        target_offsets = {(r['degree'],r['number']):r['offset'] for r in target_records}
        for edge in edges:
            block = np.load(HERE / 'sp2_j30_first_correction_blocks' / edge['folder'] / 'matrix.npy',
                            mmap_mode='r', allow_pickle=False)
            row = source_offsets[tuple(edge['source'])]
            column = target_offsets[tuple(edge['target'])]
            matrix[row:row+block.shape[0], column:column+block.shape[1]] = block
        matrix.flush()
        rhs = np.concatenate([np.load(HERE / 'sp2_j30_first_correction_rhs' /
            ('target%d_%d' % (r['degree'],r['number'])) / 'vector.npy', allow_pickle=False)
            for r in target_records])
        save_array(out / 'rhs.npy', rhs)
        result.update(matrix_sha256=digest(out / 'matrix.npy'), rhs_sha256=digest(out / 'rhs.npy'))
    atomic(out / 'result.json', json_bytes(result))
    print(json.dumps({k:v for k,v in result.items() if k not in ('sources','targets','edges')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-matrix', action='store_true')
    main(parser.parse_args().write_matrix)
