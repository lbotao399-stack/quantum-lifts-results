"""Evaluate every stored Q1 coefficient in the existing coupled-map frames.

Immutable contribution chunks make interruption/resumption exact. The input
stores 2 Q1; the inverse of two is applied once, before accumulation.
"""
import argparse
import io
import json
import struct
import time
from math import ceil, comb
from pathlib import Path

import numpy as np
import native_trace_frame as native
from fortuity.full.letters import letters_up_to
from integer_j30_boundary_system import HERE, PRIME, digest
from projected_checkpoint import atomic, json_bytes
from sp_full_projected_basis import Frame


def target_configuration(degree, number):
    graph = json.loads((HERE / 'first_correction_profile_graph.json').read_text())
    target = next(r for r in graph['target_profiles']
                  if (r['degree'], r['number']) == (degree, number))
    alphabet = sorted(L for L in letters_up_to(2, 2) if L[0] + L[1] <= 2)
    points = ceil((target['dimension'] + 32) / comb(9, degree))
    config = dict(target=[degree, number], target_dimension=target['dimension'],
                  target_profile=target['profile'], prime=PRIME,
                  exterior_generators=9, points=points,
                  columns=points * comb(9, degree),
                  seed=735103 + 7919 * degree + 104729 * number,
                  alphabet=alphabet)
    return json.loads(json.dumps(config))


def make_frame(config):
    alphabet = [(dx, dy, tuple(indices)) for dx, dy, indices in config['alphabet']]
    return Frame(2, config['points'], config['exterior_generators'],
                 [config['target'][0]], alphabet, seed=config['seed'],
                 p=PRIME, group='Sp')


def read_terms(stream, count, alphabet):
    words, coefficients = [], []
    for _ in range(count):
        length_bytes = stream.read(4)
        assert len(length_bytes) == 4
        length, = struct.unpack('<I', length_bytes)
        encoded = stream.read(length)
        coefficient_bytes = stream.read(4)
        assert len(encoded) == length and encoded[-1] == 0
        assert len(coefficient_bytes) == 4
        coefficient, = struct.unpack('<I', coefficient_bytes)
        assert 0 < coefficient < PRIME
        words.append(tuple(tuple(alphabet[i - 1] for i in trace)
                           for trace in encoded[:-1].split(b'\0')))
        coefficients.append(coefficient * pow(2, -1, PRIME) % PRIME)
    return words, np.asarray(coefficients, dtype=np.int64)


def contribution(words, coefficients, frame):
    assert len(words) * (PRIME - 1) ** 2 < 2 ** 63
    values = native.evaluate_words(words, frame)
    return np.asarray(coefficients @ values, dtype=np.int64) % PRIME


def save_array(path, values):
    stream = io.BytesIO()
    np.save(stream, values, allow_pickle=False)
    atomic(path, stream.getvalue())


def run(degree, number, seconds, chunk_size=512):
    start = time.monotonic()
    config = target_configuration(degree, number)
    source = HERE / 'sp2_j30_one_field_strength/direct_Q1_native'
    partitions = source / 'profile_partitions'
    census = json.loads((partitions / 'result.json').read_text())
    assert census['status'] == 'passed' and census['stored_multiple_of_Q1'] == 2
    assert census['profile_graph_sha256'] == digest(HERE / 'first_correction_profile_graph.json')
    record = next(r for r in census['records']
                  if (r['degree'], r['profile_number']) == (degree, number))
    path = partitions / record['file']
    assert record['sha256'] == digest(path)
    assert census['alphabet_sha256'] == digest(source / 'alphabet.txt')
    alphabet_data = [list(map(int, line.split()))
                     for line in (source / 'alphabet.txt').read_text().splitlines()]
    assert alphabet_data[0][0] == PRIME
    alphabet = [(dx, dy, tuple(i for i in range(3) if mask & (1 << i)))
                for dx, dy, mask in alphabet_data[1:]]
    validation = json.loads((HERE / 'native_trace_frame_validation.json').read_text())
    assert validation['status'] == 'passed'
    for name, value in validation['source_sha256'].items():
        assert digest(HERE / name) == value
    for marker in (HERE / 'sp2_j30_first_correction_blocks').glob('*_target%d_%d/manifest.json' % (degree, number)):
        old = json.loads(marker.read_text())
        assert all(old[key] == value for key, value in config.items()), marker
    out = HERE / 'sp2_j30_first_correction_rhs' / ('target%d_%d' % (degree, number))
    out.mkdir(parents=True, exist_ok=True)
    chunks = out / 'chunks'
    chunks.mkdir(exist_ok=True)
    manifest = dict(**config, input_sha256=record['sha256'], terms=record['terms'],
                    input_multiple=2, output_multiple=1, chunk_size=chunk_size,
                    alphabet_sha256=census['alphabet_sha256'],
                    driver_sha256=digest(Path(__file__)),
                    native_validation_sha256=digest(HERE / 'native_trace_frame_validation.json'))
    marker = out / 'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text()) == manifest
    else:
        atomic(marker, json_bytes(manifest))
    marker_hash = digest(marker)
    total = np.zeros(config['columns'], dtype=np.int64)
    frame = None
    last = time.monotonic()
    with path.open('rb') as stream:
        source_terms, terms = struct.unpack('<QQ', stream.read(16))
        assert source_terms == census['source_terms'] and terms == record['terms']
        for offset in range(0, terms, chunk_size):
            count = min(chunk_size, terms - offset)
            array_path = chunks / ('%08d.npy' % offset)
            receipt_path = array_path.with_suffix('.json')
            begin = stream.tell()
            if receipt_path.exists():
                receipt = json.loads(receipt_path.read_text())
                assert receipt['manifest_sha256'] == marker_hash
                assert receipt['offset'] == offset and receipt['count'] == count
                assert receipt['begin'] == begin
                assert receipt['sha256'] == digest(array_path)
                values = np.load(array_path, allow_pickle=False)
                assert values.shape == total.shape and np.all((0 <= values) & (values < PRIME))
                stream.seek(receipt['end'])
            else:
                if frame is None:
                    frame = make_frame(config)
                words, coefficients = read_terms(stream, count, alphabet)
                values = contribution(words, coefficients, frame)
                save_array(array_path, values)
                atomic(receipt_path, json_bytes(dict(manifest_sha256=marker_hash,
                    offset=offset, count=count, begin=begin, end=stream.tell(),
                    sha256=digest(array_path))))
            total = (total + values) % PRIME
            progress = dict(status='checkpoint', target=[degree, number],
                            completed_terms=offset + count, terms=terms,
                            seconds=time.monotonic() - start)
            atomic(out / 'progress.json', json_bytes(progress))
            if time.monotonic() - last > 20:
                print(json.dumps(progress), flush=True)
                last = time.monotonic()
            if time.monotonic() - start > seconds and offset + count < terms:
                return progress
        assert not stream.read(1)
    save_array(out / 'vector.npy', total)
    result = dict(status='complete_evaluation', target=[degree, number],
                  terms=record['terms'], coordinates=len(total),
                  nonzero_coordinates=int(np.count_nonzero(total)), prime=PRIME,
                  manifest_sha256=marker_hash, vector_sha256=digest(out / 'vector.npy'),
                  input_multiple=2, output_multiple=1, seconds=time.monotonic() - start,
                  scope='Complete Q1 right-hand side in a shared target frame. '
                        'It does not determine the coupled cohomology obstruction.')
    atomic(out / 'result.json', json_bytes(result))
    atomic(out / 'progress.json', json_bytes(result))
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--target-degree', type=int, required=True)
    parser.add_argument('--target-profile', type=int, required=True)
    parser.add_argument('--seconds', type=int, default=700)
    parser.add_argument('--chunk-size', type=int, default=512)
    args = parser.parse_args()
    run(args.target_degree, args.target_profile, args.seconds, args.chunk_size)
