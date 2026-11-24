#!/usr/bin/env python3
"""Build complete exact native input spaces for the predeclared physical block."""
from native import *
import time, argparse

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--length',type=int);args=ap.parse_args()
    start=time.monotonic();records=[]
    for record in read(HERE/'declaration.json')['native_blocks']:
        length=record['length']
        if args.length and length!=args.length:continue
        mu=tuple(record['weight_z_z_theta_theta_theta'])
        for offset in (-1,0):
            file=HERE/f'frames_Y{length}_offset{offset}.json.gz'
            if file.exists():
                saved=read(file)
            else:
                rows,wordcount=basis(mu,length+offset)
                images={d:[q(r) for r in rs] for d,rs in rows.items()}
                assert all(not q(r) for rs in images.values() for r in rs)
                saved={'weight':mu,'length':length+offset,'offset':offset,
                       'formal_word_count':wordcount,'dimensions':{d:len(rs) for d,rs in rows.items()},
                       'basis':{d:export_rows(rs) for d,rs in rows.items()},
                       'images':{d:export_rows(rs) for d,rs in images.items()}}
                write(file,saved)
            out={k:saved[k] for k in ('weight','length','offset','formal_word_count','dimensions')}
            out['sha256']=hashlib.sha256(file.read_bytes()).hexdigest();records.append(out)
            print(out,'seconds',time.monotonic()-start,flush=True)
    write(HERE/f'frames_receipt{args.length or ""}.json',{'passed':True,'records':records,'seconds':time.monotonic()-start})

if __name__=='__main__':main()
