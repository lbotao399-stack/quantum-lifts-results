"""Compare the rationally reconstructed native coordinates and two points to Dot's packet."""
import argparse
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
import gzip
import json
import numpy as np
from integer_section import graph_row,sample,trace_row

HERE=Path(__file__).resolve().parent
ROOT='f4-j36-eight-Albert-certificate/'
PIVOTS=ROOT+'f4-j36-nested1-seal/'
GRAPH=ROOT+'f4-j36-eight-cubic/'
DIGEST='b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972'


def check(archive,dense):
    assert sha256(archive.read_bytes()).hexdigest()==DIGEST
    packed=json.loads(gzip.decompress((HERE/'section-coordinates.json.gz').read_bytes()))
    basis=json.loads((HERE/'selected-1009.json').read_bytes())['rows']
    scale={'reference':1,'four':7056,'nested':592704}
    native=[Fraction(int(n)*scale[tag],int(packed['denominator'])*49787136)
            for (tag,idx),n in zip(basis,packed['numerators'],strict=True)]
    assert 49787136==28*252*84**2
    with ZipFile(archive) as z:
        for name,sub in [('cubic.json',PIVOTS),('trace_C4_ordered.json',PIVOTS),('selected-1009.json',PIVOTS)]:
            assert json.loads(z.read(sub+name))==json.loads((HERE/name).read_bytes())
        for p in (1009,1013):
            coordinates=json.loads(z.read(GRAPH+f'coordinates-{p}.json'))['coordinates'][1]
            assert [q.numerator*pow(q.denominator,-1,p)%p for q in native]==coordinates
            print(len(native),'native coordinates match modulo',p)
        if dense:
            for point in (0,1):
                x=sample(point,True)
                with np.errstate(all='ignore'):
                    t=trace_row(x);g,d=graph_row(x)
                for p in (1009,1013):
                    m=np.frombuffer(gzip.decompress(z.read(PIVOTS+f'selected-{p}.u16.gz')),dtype='<u2').reshape(703,703)
                    b=json.loads(z.read(GRAPH+f'outgoing-{p}.json'))[point][1]
                    assert [a%p for a in t]+[v*pow(s,-1,p)%p for v,s in zip(g,d)]==m[:,point].tolist()+[b]
                    print(point,p,'all 704 unscaled evaluations match')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path)
    ap.add_argument('--dense',action='store_true');a=ap.parse_args();check(a.archive,a.dense)
