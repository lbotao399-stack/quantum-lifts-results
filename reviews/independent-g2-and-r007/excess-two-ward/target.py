"""Rebuild the G2 Cartan restriction, classical ranks, and augmented product ranks."""
from argparse import ArgumentParser
from collections import Counter
import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE=Path(__file__).resolve().parent


def replay(repo, filtered, workers):
    sys.path.insert(0,str(filtered))
    sys.path.insert(0,str(repo/'algorithms/dual_singlets_20261118'))
    from build_g2_full_chain import build
    from su3_bmn_parity_products import labelled
    import g2_product_words
    data=json.loads(gzip.decompress((HERE/'g2-d2-333-rows.json.gz').read_bytes()))
    cartan=[]
    for recipe in data['recipes']:
        terms=Counter()
        for monomial,coef in recipe:
            for word,value in labelled(tuple(map(tuple,monomial)),True,'G2').items():
                terms[word]+=coef*value
        cartan.append(sum(v!=0 for v in terms.values()))
    assert cartan==[0]*8
    primes=(2147483647,1000033)
    records=[]
    with TemporaryDirectory() as scratch:
        native=Path(scratch)
        build(repo,native)
        exe=native/'span'
        subprocess.run(['c++','-O3','-std=c++17','-pthread','-I'+str(native),
                        str(filtered/'g2_product_span.cpp'),'-o',str(exe)],check=True)
        table=native/'products.dat'
        from argparse import Namespace
        g2_product_words.run(Namespace(repo=repo,recipes=HERE/'g2-d2-333-rows.json.gz',out=table))
        assert table.read_bytes()==(HERE/'g2-d2-333-rows.dat').read_bytes()
        for p,chain_seed,product_seed in zip(primes,(761291,602177),(998677,775531)):
            def command(program,args):
                return json.loads(subprocess.check_output([str(program),*map(str,args)]))
            ranks=command(native/'g2_dotted_chain',(3,0,0,6,9,chain_seed,workers,p,3,3))
            span=command(exe,(3,0,0,3,3,8,product_seed,workers,p,table))
            assert [(r['length'],r['dimension'],r['incoming'],r['outgoing']) for r in ranks['blocks']]==[
                (6,235,72,151),(7,252,151,101),(8,115,101,13),(9,13,13,0)]
            assert [span[k] for k in ('dimension','boundary','augmented','selected')]==[115,101,102,[0]]
            records.append(dict(chain=ranks,augmented=span))
    return dict(cartan_supports=cartan,checks=records)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--filtered',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('enable assertions and select a positive worker count')
    result=replay(args.repo.resolve(),args.filtered.resolve(),args.workers)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
