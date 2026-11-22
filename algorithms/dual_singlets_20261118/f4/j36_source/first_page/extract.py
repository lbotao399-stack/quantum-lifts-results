"""Extract the Lie algebra and seed tensors from Dot's pinned native archive."""
import argparse
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import tarfile
import zipfile

ARCHIVE = '0a06b306a12abbe77cadbcfab616c08c579d82086e9798aece27b8c83a0f7b52'
ROOT = 'f4-review-reproduction/'
DEPENDENCIES = ROOT + 'independent-audit/dependencies/'

def member(z, package, name):
    with tarfile.open(fileobj=io.BytesIO(z.read(DEPENDENCIES + package)), mode='r:gz') as t:
        return json.load(t.extractfile(name))

def extract(path):
    assert sha256(path.read_bytes()).hexdigest() == ARCHIVE
    with zipfile.ZipFile(path) as z:
        tables = json.loads(z.read(ROOT + 'collected-O0/color-basis-native.json'))
        base = 'f4-j30-S5-E2-exact.tar.gz'
        sixth = 'f4-j36-S6-E2-exact.tar.gz'
        active = 'f4-j36-native-active-projected-q-exact16-and-four-closed-cochains.tar.gz'
        simple = member(z, base, 'f4_j30_exact/c72-prv-basis.json')
        seeds = member(z, sixth, 'f4_j36_s6/f4-a5-prv-K.json')
        solution = member(z, sixth, 'f4_j36_s6/f4-j36-dphi-exact-solve.json')
        mirror = member(z, active, 'f4_j36_active_project_exact/f4-j30-w0-native.json')
    return dict(archive_sha256=ARCHIVE, brackets=tables['bracketPairs'], metric=tables['metric'],
                inverse=tables['inverseMetric'], raising=simple['E'], lowering=simple['F'],
                longest_weyl=mirror, seeds=[dict(i=x['i'],j=x['j'],basis=x['basis']) for x in seeds],
                coefficients=solution['coeff'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    blob = json.dumps(extract(args.archive),separators=(',',':')).encode()
    args.output.write_bytes(gzip.compress(blob,mtime=0))
