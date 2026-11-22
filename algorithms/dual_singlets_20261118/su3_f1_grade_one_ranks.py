"""Evaluate the five highest-weight sectors extracting a grade-one SU(3) F1 instanton partner."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def sectors(first):
    return [(first-1,3,3,1,1),(first,3,3,1,0),(first-1,4,3,1,0),
            (first,4,3,0,0),(first-1,4,4,0,0)]


def run(engine, output, threads, job, level30_incoming_only=False):
    output.mkdir(exist_ok=True)
    prime = 2147483647
    # Checkpoint each of the four smaller sectors before the full two-derivative sector.
    partners = [(7,8)] if level30_incoming_only else [(7,8),(8,9),(7,10),(8,11)]
    for first,length in partners:
        for i in (3,4,1,2,0):
            weight = sectors(first)[i]
            name = f'f{first}_y{length}_sector{i}.json'
            if (output/name).exists():
                continue
            seed = 830021 + 113*length + 977*first + 1021*i
            args = [str(engine.resolve()),'SU3',*map(str,weight),str(length),str(seed),str(prime)]
            row = json.loads(subprocess.check_output(args,env=dict(os.environ,PAR_THREADS=str(threads))))
            assert row['sector'] == list(weight) and row['length'] == length and row['prime'] == prime
            assert row['mode'] == 'charge_conjugation' and row['seed'] == seed
            row['job'] = job
            pending = (output/name).with_suffix('.part')
            pending.write_text(json.dumps(row,indent=2)+'\n')
            pending.replace(output/name)
            print(name,row['cochains'],row['even_in'],row['odd_in'],row['even_out'],row['odd_out'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('engine',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--threads',type=int,required=True)
    parser.add_argument('--job',required=True)
    parser.add_argument('--level30-incoming-only', action='store_true',
                        help='rebuild just the five archived level-30 incoming blocks')
    args = parser.parse_args()
    if not __debug__ or args.threads<1:
        parser.error('positive thread count and assertions required')
    run(args.engine,args.output,args.threads,args.job,args.level30_incoming_only)
