"""Write independently reproducible batches of exact integer evaluations."""
from concurrent.futures import ProcessPoolExecutor
from argparse import ArgumentParser
from multiprocessing import get_context
from pathlib import Path
import gzip,hashlib,json,os
import numpy as np
from integer_section import work,words,tensors

if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    root=parser.parse_args().out
    root.mkdir(parents=True,exist_ok=True)
    workers=int(os.environ.get('F4_INTEGER_WORKERS','32'))
    words();tensors()
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for start in range(0,720,60):
            filename=root/f'points-{start:03d}.json.gz'
            if filename.exists():
                saved=json.loads(gzip.decompress(filename.read_bytes()))
                assert list(range(start,start+60))==[i for i,v in saved] and all(len(v)==704 for i,v in saved)
                continue
            block=list(pool.map(work,((i,'all') for i in range(start,start+60)),chunksize=1))
            assert [i for i,v in block]==list(range(start,start+60)) and all(len(v)==704 for i,v in block)
            data=gzip.compress(json.dumps(block,separators=(',',':')).encode(),compresslevel=4,mtime=0)
            temp=filename.with_suffix('.part');temp.write_bytes(data);temp.replace(filename)
            print(start,len(block),len(data),hashlib.sha256(data).hexdigest(),flush=True)
