"""Replay the saved F4 determinants, G2 edge ranks and symplectic index bound."""
import gzip
import hashlib
import json
from pathlib import Path
import struct
from generator_bounds import compare
from verify import cochains, read, receipt

HERE=Path(__file__).resolve().parent


def check_minor(matrix, record, prime):
    rows,cols=record["rows"],record["columns"]
    assert len(rows)==len(set(rows))==len(cols)==len(set(cols))==record["rank"]
    square=[[matrix[i][j]%prime for j in cols] for i in rows]
    raw=b"".join(struct.pack("<I",v) for row in square for v in row)
    assert hashlib.sha256(raw).hexdigest()==record["minor_sha256"]
    determinant=1
    for j in range(len(square)):
        pivot=next(i for i in range(j,len(square)) if square[i][j])
        if pivot!=j:
            square[j],square[pivot]=square[pivot],square[j]
            determinant=-determinant
        value=square[j][j]
        determinant=determinant*value%prime
        inverse=pow(value,-1,prime)
        for i in range(j+1,len(square)):
            scale=square[i][j]*inverse%prime
            if scale:
                for k in range(j+1,len(square)):
                    square[i][k]=(square[i][k]-scale*square[j][k])%prime
    assert determinant==record["determinant"] and determinant
    return len(square)


def checked_level24():
    folder=HERE/"f4_scalar"
    previous={}
    ranks={}
    for length in (10,11):
        d=json.loads(gzip.decompress((folder/f"level24_length{length}.json.gz").read_bytes()))
        previous[length]=d["previous_evaluations"]
        check_minor(previous[length],d["previous"],d["prime"])
        ranks[length]=check_minor(previous[length]+d["inserted_values"],d["augmented"],d["prime"])
    pair=json.loads(gzip.decompress((folder/"level24_two_fermions.json.gz").read_bytes()))
    old=[row[:pair["columns"]] for row in previous[10]]
    check_minor(old,pair["previous"],pair["prime"])
    ranks[10]=max(ranks[10],check_minor(old+pair["inserted_values"],pair["augmented"],pair["prime"]))
    assert ranks=={10:80,11:61} and 141-ranks[10]-ranks[11]==0
    return ranks


def verify_stage():
    folder=HERE/"f4_scalar"
    chain=read(folder/"chain_333.json")["rows"]
    levels=cochains(9,"F4")[1]
    assert [r["dimension"] for r in chain]==[levels[y][(3,3,3)][0] for y in range(4,10)]
    assert [r["cohomology_upper"] for r in chain]==[0,0,2,0,0,1]
    for prime in (1000003,1000033):
        d=read(folder/f"scalar_{prime}.json")
        matrix=d["augmented_evaluations"]
        assert matrix[-1]==d["polynomial_evaluations"]
        assert check_minor(matrix[:-1],d["boundary_minor"],prime)==10
        assert check_minor(matrix,d["augmented_minor"],prime)==11
        b=read(folder/f"channels_{prime}.json")
        assert len(b["lower_scalar_factors"])==9
        for r in b["lower_scalar_factors"]:
            assert r["incoming"]["rank"]==r["dimension"]==levels[r["degree"]].get(tuple(r["weight"]),[0])[0]
        assert [sum(r["dimensions_by_fermions"].values())-r["incoming"]["rank"]-r["outgoing"]["rank"]
                for r in b["source_blocks"]]==[0,0,2,0]
        (a,c),(e,f)=b["exact_current_product_minor"]["rows"]
        assert a*f-c*e==b["exact_current_product_minor"]["determinant"]==-26873856
    ranks=checked_level24()
    for job in ("5691136","5691142","5691420","5691422","7658467"):
        receipt(job)
    exceptional=cochains(18,"G2")[1]
    for job in ("5691136","5691142","5691420","5691422"):
        for f in (HERE/job).glob("chain_6_*.json"):
            data=read(f)
            for row in data["blocks"]:
                assert row["dimension"]==exceptional[row["length"]][tuple(data["weight"])][0]
                assert row["dimension"]==row["incoming"]+row["outgoing"]
    c2=next(x for x in read(HERE/"5690901/symplectic.json")["rows"] if x["weight"]==[8,8,8])
    totals=[sum(v[0] for l,v in c2["cochains_by_length"].items() if int(l)%2==parity)
            for parity in (0,1)]
    assert totals==[878313,878315] and c2["primary_euler"]==[-2]
    for group in ("C2","G2"):
        p=HERE/"rank_two_indices"
        tested=compare(read(p/f"{group}_ordinary_exact.json"),read(p/f"{group}_primary_exact.json"))
        assert json.loads(json.dumps(tested))==read(p/f"{group}_bound_results.json")
    for group,shape,expected in (("F4",(3,3,3),-1),("B3",(18,18,17),-6),("C3",(18,18,17),-59)):
        table=read(HERE/"7658467"/(group+".json"))
        assert int(table["modulus"])>2*int(table["coefficient_bound"])
        assert next(row[3] for row in table["rows"] if tuple(row[:3])==shape)==expected
    return dict(F4_scalar_dimension=1,F4_possible_unprotected_incoming_sources=0,
                F4_level24_adjacent_ranks=ranks,USp4_level48_primary_index=-2)


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("Enable certificate assertions")
    print(json.dumps(verify_stage(),indent=2))
