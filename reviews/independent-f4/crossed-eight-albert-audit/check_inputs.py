"""Compare the inputs and all three six-cubic pivot recipes with Dot's pinned archive."""
from argparse import ArgumentParser
from ast import literal_eval
from fractions import Fraction
from hashlib import sha256
import json
from itertools import product
from pathlib import Path
import re
import subprocess
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ARCHIVE_HASH = 'b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972'
BASE = 'f4-j36-eight-Albert-certificate/'
BANK = BASE + 'f4-j36-eight-cubic/'
SEAL = BASE + 'f4-j36-nested1-seal/'


def native_logs(packet, side, prime):
    names = [name for name in packet.namelist() if name.startswith(BANK+'native-B-')
             and name.endswith(f'-{prime}.log')]
    found = {}
    for name in names:
        mark = Path(name).name.split('-')[2]
        if side == 'left' and mark != 'left':
            continue
        plans = {}
        for line in packet.read(name).decode().splitlines():
            if line.startswith('plan '):
                _, idx, remaining = line.split(maxsplit=2)
                _, bits, tree = literal_eval(remaining.rsplit(' ',1)[0])
                plans[int(idx)] = bits, tree
            hit = re.match(r'graph (\d+) poles \[(\d+), (\d+)\] weighted raw (\d+)', line)
            if hit:
                idx, u, v, weight = map(int, hit.groups())
                if side == 'right' and mark == 'left' and idx >= 12:
                    continue
                bits, tree = plans[idx]
                record = dict(poles=[u,v], weight=weight, bits=bits, tree=tree)
                record = json.loads(json.dumps(record))
                assert str(idx) not in found or found[str(idx)] == record
                found[str(idx)] = record
    return found


def six_cubic_pivots(packet, chart, data):
    source = packet.read(SEAL + 'B-exact-typed-replay.js').decode()
    original = 'const AX=[1,0,3,0,1,3,0,2,3,1,2,3]'
    assert source.count(original) == 1
    ends = [(0,1,2),(0,2,3),(2,0,1),(1,0,3)]
    middles = [(0,1,2),(0,1,3),(1,2,3),(0,0,1)]
    layouts = list(product(ends, middles, ends, middles))
    files = {0:'nested0-calibration-bank.json',1:'graph_bank.json'}
    graph_ids = {}
    for i, graph in enumerate(data['graphs']):
        assert graph['cubics'] == 6 and graph['cp'] == graph['cm']
        graph_ids[json.dumps([graph['E'],graph['d'],graph['cp']])] = i
    assert len(graph_ids) == len(data['graphs'])
    for pos, row in enumerate(chart['rows']):
        if row[0] != 'nested':
            continue
        recipe = row[1]
        assert recipe in (0,1,3) and recipe // 256 == 0
        axes = [axis for group in layouts[recipe] for axis in group]
        script = source.replace(original, 'const AX='+json.dumps(axes))
        built = json.loads(subprocess.check_output(['node','-e',script]))
        assert all(not proof['altResidual'] for proof in built['ward'].values())
        if recipe in files:
            assert built['bank'] == json.loads(packet.read(SEAL+files[recipe]))
        actual = []
        for terms in built['bank'][:2]:
            actual.append([[graph_ids[json.dumps([g['E'],g['d'],g['ch']])],g['w']] for g in terms])
        assert actual == data['rows'][str(pos)]


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('archive', type=Path)
    flags = cli.parse_args()
    assert sha256(flags.archive.read_bytes()).hexdigest() == ARCHIVE_HASH
    inputs = json.loads((HERE/'audit-data.json').read_bytes())
    colour = json.loads((HERE/'colour.json').read_bytes())
    pivots = json.loads((HERE/'section-graphs.json').read_bytes())
    with ZipFile(flags.archive) as packet:
        assert packet.testzip() is None
        read = lambda key: json.loads(packet.read(key))
        basis = read(BASE+'fixtures/color-basis-native.json')
        cubic = read(SEAL+'cubic.json')
        assert all(value == (cubic['cubic'] if key == 'cubic' else basis[key]) for key,value in colour.items())
        assert Fraction(cubic['k']) == Fraction(28,3)
        home = read(BANK+'bank-B-reduced.json')
        for side, graphs in inputs['graphs'].items():
            assert not home[side]['ordinary'] and len(home[side]['banks']) == len(graphs)
            for g, original in zip(graphs,home[side]['banks'],strict=True):
                assert all(g[k] == original[k] for k in ('E','d','cp','cm','w','cubics'))
        for p in (1009,1013):
            coord = inputs['coordinates'][str(p)]
            assert coord['coefficients'] == read(BANK+f'coordinates-{p}.json')['coordinates'][1]
            assert coord['outgoing'] == [row[1] for row in read(BANK+f'outgoing-{p}.json')]
            assert coord['chart_sha256'] == sha256(packet.read(SEAL+f'selected-{p}.u16.gz')).hexdigest()
            for side, observations in inputs['observations'][str(p)].items():
                assert native_logs(packet,side,p) == observations
        traces = read(SEAL+'trace_C4_ordered.json')['rows']
        scores = [-104166480711910256640000,-17361080118651709440000,-109953507418127493120000,0]
        column = [sum(Fraction(v)*s for v,s in zip(row,scores,strict=True)) for row in traces]
        chart = read(SEAL+'selected-1009.json')
        pivot = chart['rows'].index(['nested',0])
        assert len(column) == 703 and column[pivot] == 0
        column[pivot] = -Fraction(598031033628600,7)
        assert {str(i):str(v) for i,v in enumerate(column) if v} == inputs['section']
        six_cubic_pivots(packet,chart,pivots)
    print('Pinned colour tensors, 48 graph entries, individual logs, classical rows, and three regenerated pivot recipes agree.')


if __name__=='__main__':
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    main()
