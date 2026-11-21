"""Verify the protected product spans and both sides of the rational Taylor dimension bounds."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from rank_two_taylor_products import PRIME, SECTORS, record, replay
from taylor_hook_cohomology import resolve, weights


ROOT = Path(__file__).resolve().parent/"taylor_hooks"


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix==".gz" else raw)


def check_receipt(folder,accounting):
    record = read(folder/"receipt.json")
    job,state,_,_,exit_code = accounting.split("|")
    assert job==record["job"] and state=="COMPLETED" and exit_code=="0:0"
    for name,digest in record["files"].items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest


def calculate(base, replay_minors=False, rebuild_ranks=False):
    folder = base/"cartan_products"
    jobs = read(folder/"jobs.json")
    for stage in ("weights","raised_spin"):
        check_receipt(folder/stage,jobs[stage])
    symplectic = []
    for k in ("hwsp0","hwsp1","hwsprest"):
        stage = folder/"symplectic_classical"/k
        check_receipt(stage,jobs[k])
        symplectic.extend(read(p) for p in sorted(stage.glob("Sp2*.json")))
    symplectic.sort(key=lambda r:r["sector_index"])
    blocks = {}
    for case,(group,tail,length) in SECTORS.items():
        sectors = weights(5,tail,True)
        if group=="Sp2":
            native = symplectic
        else:
            classical = base/("colour_parity/2147483647" if case=="su30" else "su3_level32/second_target")
            native = [read(classical/f"sector_{i}.json") for i in range(5)]
        rows = []
        for i,sector in enumerate(sectors):
            assert native[i]["prime"]==PRIME and native[i]["sector"]==list(sector)
            ranks = []
            for name in ("weights","raised_spin"):
                if name=="raised_spin" and not sector[4]:
                    ranks.append(0)
                    continue
                data = read(folder/name/(case+f"_{i}.json.gz"))
                expected = sector if name=="weights" else (*sector[:3],sector[3]+1,sector[4]-1)
                assert (data["sector"],data["prime"],data["length"],data["case"],data["index"]) == (
                    list(expected),PRIME,length,case,i)
                assert data.get("raised_spin",False)==(name=="raised_spin")
                assert 0<data["det"]<PRIME and len(data["selected_recipes"])==len(data["selected_monomials"])==data["rank"]
                if rebuild_ranks:
                    computed = record(case, i, name == "raised_spin")
                    assert json.loads(json.dumps(computed)) == data
                if replay_minors:
                    replay(data)
                ranks.append(data["rank"])
            upper = (native[i]["blocks"][0]["h"] if group=="Sp2"
                     else native[i]["even_h"]+native[i]["odd_h"])
            assert ranks[0]-ranks[1]==upper
            rows.append(dict(native_sector=sector,weight_image=ranks[0],raised_weight_image=ranks[1],
                             spin_highest_product_rank=upper,classical_upper=upper))
        target = sum((-1)**i*r["spin_highest_product_rank"] for i,r in enumerate(rows))
        if case=="sp36":
            dimensions = resolve(native)["blocks"][0]
        else:
            parity = "odd" if length%2 else "even"
            resolved=read(classical/"resolved.json")
            dimensions={s:sum((-1)**i*r[parity+"_"+s] for i,r in enumerate(native)) for s in ("dim","in","out","h")}
            assert resolved["parities"][parity]==dimensions
            assert dimensions["dim"]-dimensions["in"]-dimensions["out"]==target
        assert dimensions["h"]==target
        blocks[case]=dict(group=group,partition=[5,5]+[1]*tail,length=length,native_weights=rows,
                          classical=dimensions,Cartan_detected_products=target)
    return blocks


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",type=Path,default=ROOT)
    p.add_argument("--replay-minors",action="store_true")
    p.add_argument("--rebuild-ranks",action="store_true",
                   help="enumerate all products in each ordinary and raised weight and check the complete modular ranks")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    data=calculate(args.base,args.replay_minors,args.rebuild_ranks)
    args.out.write_text(json.dumps(data,indent=2)+"\n")
    print({k:r["Cartan_detected_products"] for k,r in data.items()},flush=True)
