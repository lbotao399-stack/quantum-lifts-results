"""Test the prior rank-three graviton polynomial in the defining seven of G2."""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


def convert(directory, target):
    data = gzip.decompress((directory/'product.txt.gz').read_bytes())
    report = json.loads((directory/'product.json').read_bytes())
    assert sha256(data).hexdigest() == report['native_trace_sha256']
    assert report['common_denominator'] == 540 and report['rank_three_Cartan_terms']==0
    rows = iter([int(i) for i in s.split()] for s in data.splitlines())
    letters, traces, words = next(rows)
    assert (letters, traces, words)==(7,149,1212)
    alphabet = [next(rows) for _ in range(letters)]
    coefficients = []
    with target.open('w') as stream:
        print(traces, words, 1, file=stream)
        for _ in range(traces):
            n, *term = next(rows)
            assert len(term)==n
            print(n,*(v for j in term for v in alphabet[j]),file=stream)
        for _ in range(words):
            c, n, *term = next(rows)
            assert len(term)==n
            print(n,*term,file=stream)
            coefficients.append(c)
        print(words,*(f'{i} {c} 540' for i,c in enumerate(coefficients)),file=stream)
    assert next(rows,None) is None


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--engine',type=Path,required=True)
    p.add_argument('--prime',type=int,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--group',choices=('G2','SO7','Sp3'),default='G2')
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    if not __debug__ or args.workers<1:raise ValueError('positive worker count and assertions required')
    case=args.repo/'algorithms/dual_singlets_20261118/taylor_hooks/bmn_633'
    with TemporaryDirectory() as scratch:
        packet=Path(scratch)/'product.dat'
        convert(case,packet)
        cmd=[args.engine,6,3,3,0,0,10,args.seed,args.workers,args.prime,packet]
        if args.group!='G2':cmd.append(args.group)
        result=json.loads(subprocess.check_output(list(map(str,cmd)),timeout=620))
    assert result['products']==1 and result['prime']==args.prime
    assert result['augmented']-result['boundary']==len(result['selected'])
    assert set(result['selected'])<={0}
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
