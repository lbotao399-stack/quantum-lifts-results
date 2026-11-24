"""Rebuild B05 from trace enumeration, integer colour frames and exact maps."""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path

HERE=Path(__file__).resolve().parent


def run(work,receipt):
    work.mkdir(parents=True,exist_ok=True)
    assert not any(work.iterdir()),"Use an empty reconstruction directory."
    start=time.monotonic();stages=[]
    def call(script,*args):
        before=time.monotonic();command=[sys.executable,"-B",str(HERE/script),*map(str,args)]
        subprocess.run(command,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1"))
        stages.append(dict(script=script,seconds=time.monotonic()-before,exit_code=0))
    call("characters.py")
    for length in (6,7,8):call("prepare.py",length,"--output",work/f"stable{length}.pkl")
    for length in (6,7,8):call("finite_frame.py",work/f"stable{length}.pkl",work/f"frame{length}")
    for length in (6,7):call("maps.py",work,length)
    call("products.py",HERE/"certificate/proper_product_recipes.json",work/"products.pkl")
    call("extract.py",work)
    call("relation_test.py",work)
    call("export.py",work,work/"certificate")
    core=["basis.json.gz","maps.json.gz","frames.json.gz","generator.json.gz","rank_certificates.json",
          "proper_product_recipes.json","complete_graviton_recipes.json","summary.json","closure_relation.json"]
    hashes={}
    for filename in core:
        expected=(HERE/"certificate"/filename).read_bytes();actual=(work/"certificate"/filename).read_bytes()
        assert actual==expected,filename
        hashes[filename]=hashlib.sha256(actual).hexdigest()
    call("verify.py",work/"certificate","--skip-manifest","--receipt",work/"archived_check.json")
    call("independent_verify.py",work/"certificate","--receipt",work/"independent_check.json")
    result=dict(clean_source_rebuild=True,private_inputs=False,source_enumeration_rebuilt=True,
                all_finite_relations_reconstructed=True,all_native_polynomial_maps_reconstructed=True,
                independent_full_polynomial_residual_verified=True,bit_identical_core_artifacts=hashes,
                stages=stages,seconds=time.monotonic()-start)
    if receipt:receipt.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path);parser.add_argument("--receipt",type=Path)
    args=parser.parse_args();assert __debug__;run(args.work,args.receipt)
