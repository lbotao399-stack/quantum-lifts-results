"""Classical spin singlets that can receive a differential from a BMN flavour singlet."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess


def run(engine, power, out, threads):
    out.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(engine.read_bytes()).hexdigest()
    groups = ("Sp2","Sp3","SO7")
    def compute(group):
        seed = 783721 + 151*power + groups.index(group)
        prime = 2147483647
        environment = dict(os.environ, PAR_THREADS=str(max(1,threads//len(groups))),
                           A_VERB="1", A_NORHO="1")
        arguments = [str(engine),"chain",group,"0","0","0",str(power),str(power),
                     "2",str(2*power),str(seed),str(prime)]
        with (out.parent/(group+"_native.log")).open("w") as log:
            raw = subprocess.check_output(arguments, env=environment, stderr=log)
        data = json.loads(raw)
        assert data["sector"]==[0,0,0,power,power] and data["gauge"]==group
        assert [b["Y"] for b in data["blocks"]]==list(range(2,2*power+1))
        assert all(b["dim"]-b["rank_in"]-b["rank_out"]==b["h"]>=0 for b in data["blocks"])
        for left,right in zip(data["blocks"],data["blocks"][1:]):
            assert left["rank_out"]==right["rank_in"]
        data.update(engine_sha256=digest, arguments=arguments[1:])
        result = out/(group+".json")
        stage = result.with_suffix(".part")
        stage.write_text(json.dumps(data,indent=2)+"\n")
        stage.replace(result)
        print(group, power, [(b["Y"],b["h"]) for b in data["blocks"]], flush=True)
        return result
    with ThreadPoolExecutor(len(groups)) as pool:
        return list(pool.map(compute,groups))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine",type=Path,required=True)
    p.add_argument("--power",type=int,choices=(5,6),required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--threads",type=int,default=6)
    a=p.parse_args()
    run(a.engine.resolve(),a.power,a.out,a.threads)
