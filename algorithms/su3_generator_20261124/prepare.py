"""Build the complete stable dotted-singlet source at the specified length."""
import argparse
import json
from pathlib import Path
import pickle
import time
from trace_algebra import highest_weights, differential, skeleton

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("length",type=int)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args(); start=time.monotonic()
    rows=highest_weights(args.length,progress=lambda value:print(json.dumps(value),flush=True))
    images=list(map(differential,rows))
    assert all(not differential(row) for row in images)
    result=dict(length=args.length,parity=-1,weight=[3,0,0,4,4],sources=rows,images=images)
    args.output.write_bytes(pickle.dumps(result,protocol=5))
    print(json.dumps(dict(length=args.length,singlets=len(rows),skeletons=len({skeleton(row) for row in rows}),source_terms=sum(map(len,rows)),image_terms=sum(map(len,images)),seconds=time.monotonic()-start)),flush=True)
