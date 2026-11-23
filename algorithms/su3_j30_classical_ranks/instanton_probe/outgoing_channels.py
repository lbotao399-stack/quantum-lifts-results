"""Verify the odd R-007 outgoing channels and, optionally, regenerate their native minors."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ALG = HERE.parents[1] / "dual_singlets_20261118"
sys.path.insert(0, str(ALG))
from trace_descent import dimensions


def checked_minor(record, field, height, width):
    minor = record[field]
    rank = minor["rank"]
    assert 0 < minor["det"] < record["prime"]
    for key, bound in (("rows", height), ("columns", width)):
        chosen = minor[key]
        assert len(chosen) == len(set(chosen)) == rank
        assert all(0 <= i < bound for i in chosen)
    return rank


def verify(rebuild, threads):
    archive = json.loads(gzip.decompress((HERE / "outgoing_minors.json.gz").read_bytes()))
    result = {}
    with TemporaryDirectory() as tmp:
        binaries = {}
        for name, entry in archive.items():
            data = entry["data"]
            raw = (json.dumps(data, separators=(",", ":")) + "\n").encode()
            assert hashlib.sha256(raw).hexdigest() == entry["sha256"]
            t = data["sector"][-1]
            weight = (t,t,7-t,4-t,4-t)
            assert data["gauge"] == "SU3" and data["j"] == 30 and t in (1,3,4)
            assert data["sector"] == [*weight[2:],*weight[:2]]
            y = data["length"]
            prime, seed = data["prime"], data["seed"]
            if data["mode"] == "C_odd_projected_classical_image":
                for l, side in ((y,"source"),(y+1,"target")):
                    odd = sum(dimensions(weight,l,3,"SU-").values())
                    even = sum(dimensions(weight,l,3,"SU+").values())
                    assert (data["odd_"+side],data["full_"+side]) == (odd,odd+even)
                result[name] = checked_minor(data,"minor",data["odd_target"],data["odd_source"])
                driver = "odd_colour_minor"
                args = [*data["sector"],y,prime,seed]
            else:
                assert data["mode"] == "classical_Cartan_kernel" and y==7 and t==4
                for l, field in ((6,"parent"),(7,"cochains"),(8,"next")):
                    assert sum(sum(dimensions(weight,l,3,g).values()) for g in ("SU-","SU+")) == data[field]
                result[name] = [checked_minor(data,"incoming",data["cochains"],data["parent"]),
                                checked_minor(data,"outgoing",data["next"],data["cochains"]),
                                checked_minor(data,"outgoing_and_Cartan",data["stacked_rows"],data["cochains"])]
                driver = "su2_cartan_kernel"
                args = [*data["sector"],y,prime,seed,"SU3"]
            if rebuild=="full" or (rebuild=="small" and t==4):
                if driver not in binaries:
                    binaries[driver] = Path(tmp)/driver
                    source = ALG/"taylor_hooks"
                    subprocess.run(["c++","-O3","-std=c++17","-pthread","-I"+str(source/"native"),
                                    str(source/(driver+".cpp")),"-o",str(binaries[driver])],check=True)
                output = subprocess.check_output([str(binaries[driver]),*map(str,args)],
                                                 env=dict(os.environ, PAR_THREADS=str(threads)))
                assert output == raw
    assert result["n4r7oddin-final"] + result["n4r7oddout-final"] == 34430
    assert result["t3-7"] + result["t3-8"] == 26400
    assert result["t4-in"] == result["t4-in-second"] == 825
    assert result["t4-out"] == result["t4-out-second"] == 1037
    assert result["t4-5"]+825 == 1169 and 1037+result["t4-8"] == 1588
    assert result["t4-stacked"] == result["t4-stacked-second"] == [1698,2026,2060]
    assert 1896-825-1037 == 34 and 3758-1698-2060 == 0
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rebuild",choices=("small","full"))
    p.add_argument("--threads",type=int,default=4)
    a=p.parse_args()
    if not __debug__ or a.threads<1: p.error("enable assertions and choose a positive thread count")
    print(json.dumps(verify(a.rebuild,a.threads),indent=2))
