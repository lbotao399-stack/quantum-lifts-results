"""Solve the crossed-source classical section over Q and verify every integer equation."""
from fractions import Fraction
from functools import reduce
from math import gcd,lcm
import argparse
import gzip
import json
from pathlib import Path

from flint import fmpz_mat, nmod_mat

HERE=Path(__file__).resolve().parent
ANCHORS=(-104166480711910256640000,-17361080118651709440000,-109953507418127493120000,0)
NESTED=Fraction(-598031033628600,7)
ORDERED=Fraction(-299015516814300,49)


def check(folder, output, saved):
    rows=[]
    for f in sorted(folder.glob('points-*.json.gz')):
        rows.extend(json.loads(gzip.decompress(f.read_bytes())))
    assert [i for i,r in rows]==list(range(720)) and all(len(v)==704 for i,v in rows)
    data=[v for _,v in rows]
    # Rank and pivot points are determined over the prime field. The same nonzero
    # integer determinant then makes the square interpolation unique over Q.
    p=1009
    rref,rank=nmod_mat([v[:-1] for v in data],p).transpose().rref()
    assert rank==703,rank
    indices=[next(i for i in range(len(data)) if rref[k,i]) for k in range(rank)]
    mat=[data[i][:-1] for i in indices]
    det=int(nmod_mat(mat,p).det())
    assert det
    tags=json.loads((folder/'selected-1009.json').read_bytes())['rows']
    # The primitive native corner scales are 28,252 and 84 for low, high and exchange.
    # The trace rows are integral; four/nested/crossed contain zero/one/two exchanges.
    scale={'reference':1,'four':28*252,'nested':28*252*84}
    target_scale=28*252*84**2
    scales=[scale[tag] for tag,_ in tags]
    if saved is None:
        rational=fmpz_mat(mat).solve(fmpz_mat([[data[i][-1]] for i in indices]))
        co=[Fraction(str(v)) for v in rational.entries()]
        denominator=lcm(*(v.denominator for v in co))
        nums=[int(v*denominator) for v in co]
        encoding=dict(denominator=str(denominator),numerators=list(map(str,nums)))
        (output/'section-coordinates.json.gz').write_bytes(gzip.compress(json.dumps(encoding,separators=(',',':')).encode(),compresslevel=9,mtime=0))
    else:
        encoding=json.loads(gzip.decompress(saved.read_bytes()))
        denominator=int(encoding['denominator']);nums=list(map(int,encoding['numerators']))
    assert len(nums)==703 and denominator>0 and reduce(gcd,nums,denominator)==1
    assert all(sum(a*b for a,b in zip(v[:-1],nums))==v[-1]*denominator for v in data)
    c4=json.loads((folder/'trace_C4_ordered.json').read_bytes())['rows']
    sc=[sum(Fraction(x)*y for x,y in zip(row,ANCHORS)) for row in c4]
    sc[tags.index(['nested',0])]=NESTED
    physical=[Fraction(n*s,denominator*target_scale) for n,s in zip(nums,scales)]
    result=sum(a*b for a,b in zip(physical,sc));remainder=ORDERED-result
    def mod(q,p):return q.numerator*pow(q.denominator,-1,p)%p
    assert [mod(result,p) for p in (1009,1013)]==[44,92]
    cert=dict(evaluation_points=len(data),rank=rank,selected_points=indices,scaled_minor_mod_1009=det,
              scaled_minor_mod_1013=int(nmod_mat(mat,1013).det()),exact_equations=720,
              nonzero_section_coordinates=sum(bool(q) for q in nums),
              normalized_crossed_score=str(ORDERED),rational_classical_section=str(result),
              corrected_rational_pairing=str(remainder),common_coordinate_denominator_digits=len(str(denominator)))
    (output/'section-result.json').write_text(json.dumps(cert,indent=2)+'\n')
    print({k:v for k,v in cert.items() if k!='selected_points'},flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--inputs',type=Path,default=HERE)
    ap.add_argument('--out',type=Path,default=HERE);ap.add_argument('--verify',type=Path)
    a=ap.parse_args();check(a.inputs,a.out,a.verify)
