"""Compare Dot's published and Codex's ordered Q1 maps on the SO(11) chain."""
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

from terminal_primitive import GHOSTS, HERE, LETTERS, PHI, Q1d_multi, traces
from fortuity.full.letters import Q_letter
from q1d_exact import pair_terms

sys.path.insert(0, str(HERE.parent))
from verify_dot_coordinates import read_bundle


def check_classical_letters():
    # Dot uses an oriented bracket for each second jet. Codex uses divided
    # derivatives. The double identical ghost has coefficient two in Dot.
    pairs = ['BC', 'AC', 'AB', '--', '-+', '++', '-A', '-B', '-C', '+A', '+B', '+C']
    for labels in pairs:
        x, y = (LETTERS[c] for c in labels)
        jet = (x[0]+y[0], x[1]+y[1], tuple(sorted(x[2]+y[2])))
        classical = {(a, b): c for c, a, b in Q_letter(jet)}
        dot = Counter({(x, y): 1})
        both_odd = not x[2] and not y[2]
        dot[y, x] += 1 if both_odd else -1
        assert dot == {k: v * (2 if x == y else 1) for k, v in classical.items()}
    return len(pairs)


def check_ordered_pairs():
    jets = dict(LETTERS)
    for i, symbol in enumerate('abc'):
        jet = (0, 0, tuple(j for j in range(3) if i != j))
        jets[symbol] = jet
        sign = -((-1)**i)
        for pair in ((jet, PHI[i]), (PHI[i], jet)):
            assert {(x,y): c for c,x,y in pair_terms(*pair)} == {
                GHOSTS: sign, GHOSTS[::-1]: -sign}
        assert all(not pair_terms(jet, phi) for j, phi in enumerate(PHI) if j != i)
    return jets


def probes():
    samples = []
    for psi, scalar in zip('abc', 'ABC'):
        samples.append([psi+scalar])
        for x, y in product('+-AB', repeat=2):
            samples.extend(([psi+x+scalar+y], [x+psi+y+scalar], [psi+x, scalar+y]))
        samples.extend(([psi+scalar+scalar, 'ABC'], [psi+'AB', 'C'+scalar],
                        ['+A', psi+scalar, '-B'], ['A-', psi+'+', 'B'+scalar]))
    return samples


def main(folder):
    classical_pairs = check_classical_letters()
    labels = check_ordered_pairs()
    samples = probes()
    data = read_bundle(folder, 'so11-j36')
    engine = data['so11-j36-inert-bank-streaming.mjs'].decode().split('// b4=', 1)[0]
    node = """
const fs = require('fs'), vm = require('vm');
const input = JSON.parse(fs.readFileSync(0, 'utf8'));
const emit = vm.runInNewContext(input.engine + '\\nsoQ1emit4');
const output = input.samples.map(word => {
    const row = [];
    emit(word, 11, (terms, value) => row.push([terms, value]));
    return row;
});
process.stdout.write(JSON.stringify(output));
"""
    raw = subprocess.check_output(['node', '-e', node],
                                  input=json.dumps(dict(engine=engine, samples=samples)).encode())
    def decode(word):
        return tuple(tuple(labels[c] for c in t) for t in word)

    nonzero = 0
    for sample, emitted in zip(samples, json.loads(raw), strict=True):
        dot = Counter()
        for word, value in emitted:
            sign, key = traces.so_canon(decode(word))
            if sign:
                dot[key] += sign * value
        codex = Counter()
        sign, key = traces.so_canon(decode(sample))
        if sign:
            for word, value in Q1d_multi(key, 11, 'SO_odd').items():
                orientation, key_out = traces.so_canon(word)
                if orientation:
                    codex[key_out] += sign * orientation * value
        assert all(dot[key] == -2 * codex[key] for key in dot.keys() | codex.keys()), sample
        nonzero += any(codex.values())
        if sample == ['aA']:
            _, target = traces.so_canon((GHOSTS[::-1],))
            assert (dot[target], codex[target]) == (72, -36)
    assert nonzero
    converted = (-2)**3 * 3072000
    assert 9 * converted + 25600 * 8640 == 0
    return dict(classical_pairs=classical_pairs, ordered_matter_pairs=6,
                trace_inputs=len(samples), nonzero_trace_images=nonzero,
                dot_raw_per_codex=-2, dot_normalized_per_codex='-1/2',
                dot_third_in_codex_coupling=converted, codex_third=8640)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle_directory', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    print(json.dumps(main(args.bundle_directory)))
