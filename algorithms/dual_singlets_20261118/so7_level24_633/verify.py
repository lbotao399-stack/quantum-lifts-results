"""Check both integer target frames and the finite-field comparison at doubled charges (2,2,12,6,6)."""
import gzip
import hashlib
import json
from pathlib import Path

from flint import fmpz_mat, nmod_mat


def run(folder):
    for location in ("a8","local"):
        page=json.loads((folder/(location+"_pairing.json")).read_bytes())
        raw=(folder/(location+"_coordinates.json.gz")).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==page["coordinates_sha256"]
        coordinates=json.loads(gzip.decompress(raw))
        boundaries=coordinates["classical_images"]
        target=coordinates["candidate"]
        numerator=coordinates["first_numerator"]
        denominator=coordinates["first_denominator"]
        assert (len(boundaries),len(target),denominator)==(476,250,126)
        residual=[a+denominator*b for a,b in zip(numerator,target,strict=True)]
        for extension,rank in (([],209),([target],210),([numerator],210),
                               ([target,numerator],210),([residual],209)):
            assert fmpz_mat(boundaries+extension).rank()==rank
        augmented=boundaries+[target]
        for z in page["finite_field_controls"]:
            minor=[[augmented[i][j] for j in z["columns"]] for i in z["rows"]]
            assert z["determinant"]==int(nmod_mat(minor,z["prime"]).det())!=0
        print(location,page["seed"],"rational ranks",209,210,210,210,209)
    for group in ("so","sp"):
        outgoing,incoming=[json.loads((folder/(group+str(n)+".json")).read_bytes()) for n in (1,2)]
        for band in (outgoing,incoming):
            assert band["flavours"]==[6,3,3] and band["prime"]==1000003
            intervals=band["interval_ranks_and_pairings"]
            ranks={(a,b):r for a,b,r,_ in intervals}
            for a,b,r,multiplicity in intervals:
                assert multiplicity==r-ranks.get((a+1,b),0)-ranks.get((a,b-1),0)+ranks.get((a+1,b-1),0)
            assert band["ranks_by_order"]==[sum(n for a,b,r,n in intervals if b-a==k) for k in range(4)]
        middle=sum(incoming["target_dimensions"])
        assert middle-sum(outgoing["ranks_by_order"])-sum(incoming["ranks_by_order"])==150
        print(group,middle,"minus ranks",sum(incoming["ranks_by_order"]),sum(outgoing["ranks_by_order"]),"= 150")


if __name__=="__main__":
    assert __debug__
    run(Path(__file__).resolve().parent)
