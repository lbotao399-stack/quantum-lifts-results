"""Rebuild Dot's SU(3) endpoints and the three-scalar-jet Cartan probe."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import gzip
import importlib.util
import io
import json
from math import prod
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import zipfile

PIN = '8ba864d4eb9413131922b75aa653ab119d10a04a'
BASE = 'nonperturbative-lift/contributions/'
REPO = Path(__file__).resolve().parents[3]
HASHES = {
    'replay.py':'a6bc8ff37a271b0e49f04fba58622a146e88aa490558a875c12b66dadfa801de',
    'R007_recipe.json':'c82ec0745328b8482c12dec843bb25d25345f9f7dffa6776720a37971acba4df',
    'su30_0.json.gz':'81f936300745c652c4279bb769365b59df01fdcb321860d657b7b6d3e1188ebc',
}
A, B, C = (0,0,(0,)), (0,1,(1,)), (1,0,(2,))


def blob(path):
    return subprocess.check_output(['git','-C',str(REPO),'show',PIN+':'+BASE+path])


def endpoint(engine, inputs):
    products = json.loads(gzip.decompress(inputs['su30_0.json.gz']))['selected_recipes']
    recipe = json.loads(inputs['R007_recipe.json'])['source_linear_combination']
    terms = Counter()
    for index,c in recipe:
        part = {():Fraction(c)}
        for factor in products[index]:
            part = engine.multiply(part,engine.factor(tuple(factor)))
        terms.update(part)
    target = engine.clean(terms)
    for spin,flavour in [(0,1),(1,1),(0,2),(1,2)]:
        target = engine.act(target,('T',spin,flavour))
    assert len(target) == 5080
    return target


def verify_source():
    raw_zip = blob('dot-r007-exact-input/phase34-public-q0-q1-replay.zip')
    with zipfile.ZipFile(io.BytesIO(raw_zip)) as archive:
        raw = next(archive.read(p) for p in archive.namelist() if p.endswith('/dmk_o330_hw.json'))
    assert sha256(raw).hexdigest() == '8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34'
    terms = json.loads(raw)['terms']
    for _,word in terms:
        letters = sum(word,[])
        assert len(letters) == 10 and all(0 <= i < 6 for i in letters)
        assert sum(i < 3 for i in letters) == 5
    assert len(terms) == 6139
    return len(terms)


def verify_target(engine, target):
    abelian_line = Counter()
    profiles = Counter()
    selected = {}
    labels = {A:'A',B:'B',C:'C'}
    for w,c in target.items():
        letters = sum(w,())
        assert sorted(map(len,w)) == [2,2,2,3]
        sign, monomial = engine.ordered(letters, engine.odd)
        if sign:
            abelian_line[monomial] += sign*c
        counts = tuple(sum(len(l[2]) == n for l in letters) for n in range(4))
        profiles[counts] += 1
        if set(letters) <= labels.keys():
            selected[tuple(''.join(labels[l] for l in tr) for tr in w)] = c
    assert len(abelian_line) == 566 and not any(abelian_line.values())
    assert selected == {
        ('AA','AAA','BB','CC'):Fraction(8,5),
        ('AA','AAA','BC','BC'):Fraction(-8,5),
        ('AAA','AB','AB','CC'):Fraction(32,15),
        ('AAA','AB','AC','BC'):Fraction(-64,15),
        ('AAA','AC','AC','BB'):Fraction(32,15),
    }
    H = {'A':(1,1,-2), 'B':(1,0,-1), 'C':(0,1,-1)}
    value = sum(c*prod(sum(prod(H[l][i] for l in tr) for i in range(3)) for tr in w)
                for w,c in selected.items())
    assert value == Fraction(-2016,5)
    degree = min(n0+2*n1+n2 for n0,n1,n2,n3 in profiles)
    assert degree == 10 and sum(profiles.values()) == 5080
    return dict(target_terms=5080, collinear_coefficients=566, nonzero_collinear=0,
                profile_counts=[[list(k),v] for k,v in sorted(profiles.items())],
                three_jet_traces=[[' '.join(w),str(c)] for w,c in sorted(selected.items())],
                cartan_pairing=str(value), minimum_target_mode_degree=degree)


def main():
    source_count = verify_source()
    inputs = {name:blob('dot-r007-t2-primary/'+name) for name in HASHES}
    assert all(sha256(data).hexdigest() == HASHES[name] for name,data in inputs.items())
    with TemporaryDirectory() as tmp:
        script = Path(tmp)/'replay.py'
        script.write_bytes(inputs['replay.py'])
        spec = importlib.util.spec_from_file_location('dot_cartan_candidate',script)
        engine = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(engine)
        result = verify_target(engine,endpoint(engine,inputs))
    result.update(source_terms=source_count, minimum_source_mode_degree=15, one_instanton_modes=24)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with assertions enabled')
    main()
