"""Exact finite-rank character dimensions in the archived omega grading."""
from pathlib import Path
import argparse,json,time
from rank_runtime import m
from trace_descent import dimensions
from mixed_incoming import sector

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes',type=int,nargs='+',default=[6,8,10],help='Defining dimensions 2N of Sp(N)')
    parser.add_argument('--lengths',type=int,nargs='+',default=[8,9,10])
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args()
    assert all(size>=2 and size%2==0 for size in args.sizes)
    start=time.monotonic();weight=(3,3,5,0,0);records=[]
    for size in args.sizes:
        for length in args.lengths:
            before=time.monotonic()
            d=dimensions(weight,length,size,'Sp')
            row={'group':'Sp','defining_dimension':size,'weight':weight,'length':length,
                 'fermion_degrees':d,'dimension':sum(d.values()),'sector':sector(weight,length),
                 'seconds':time.monotonic()-before}
            records.append(row)
            print(json.dumps(row),flush=True)
    args.output.write_text(json.dumps({'records':records,'seconds':time.monotonic()-start},indent=2)+'\n')

if __name__=='__main__':main()
