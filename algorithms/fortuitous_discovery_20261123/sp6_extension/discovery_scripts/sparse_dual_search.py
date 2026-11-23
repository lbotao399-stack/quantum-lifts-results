"""Search exact one- and two-coefficient dual witnesses."""
from fast_paths import *
from fractions import Fraction
import argparse,json,math,pickle,time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,default=5)
    args=ap.parse_args();rank=args.rank;start=time.monotonic()
    folder=HERE/f'frames_sp{rank}'
    defect=pickle.loads((HERE/f'small_sp{rank}/candidate_defect.pkl').read_bytes())
    source=[pickle.loads((folder/f'Q{2*rank+3}'/f'row_{i}.pkl').read_bytes())[1] for i in range(188)]
    universe=sorted(set(defect).union(*(r.keys() for r in source)),key=lambda w:(len(w),w))
    seen={};found=None;records=[]
    for word in universe:
        column=tuple(row.get(word,0) for row in source)
        divisor=math.gcd(*column)
        if not divisor:
            if defect.get(word,0):found=[(word,Fraction(1,defect[word]))];break
            continue
        if next(c for c in column if c)<0:divisor=-divisor
        signature=tuple(c//divisor for c in column)
        value=Fraction(defect.get(word,0),divisor)
        if signature in seen:
            old,old_divisor,old_value=seen[signature]
            if value!=old_value:
                found=[(word,1/(divisor*(value-old_value))),
                       (old,-1/(old_divisor*(value-old_value)))];break
        else:seen[signature]=(word,divisor,value)
    if found:
        assert all(sum(c*r.get(w,0) for w,c in found)==0 for r in source)
        assert sum(c*defect.get(w,0) for w,c in found)==1
    result=dict(rank=rank,coefficient_probes_examined=len(universe),
                proportional_column_signatures=len(seen),found=found is not None,
                dual=[(str(c),pf.m.encode(w)) for w,c in found] if found else None,
                exact_rational=True,seconds=time.monotonic()-start)
    (HERE/f'sparse_dual_sp{rank}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
