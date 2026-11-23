"""Two-stage highest weights for orthogonal BMN trace/Pfaffian graphs.

Impose E1 inside its invariant shape classes, then E2 inside each full
flavour shape. This is the exact intersection of the two kernels over
the stated field, without the large simultaneous raising matrix.
"""
import hw_scan as base
from hw_scan import *
from hw_modular import BoundedFrame


def e1_shape(inv):
    labels=(0,0,2,3,3,5,6)
    out=[]
    for typ,seq in enumerate(inv):
        words=[]
        for w in seq:
            t=tuple(labels[x] for x in w)
            words.append(min(t,t[::-1]) if typ==0 else
                min(v[i:]+v[:i] for v in (t,t[::-1]) for i in range(len(v))))
        out.append(tuple(sorted(words)))
    return tuple(out)


def raise_expr(e,i,p):
    out=defaultdict(int)
    for w,c in e.items():
        for v,d in raising(w,i).items():out[v]+=int(c)*d
    return {w:c%p for w,c in out.items() if c%p}


def kernel(exprs,i,p):
    H,_=formal_closed_kernel([raise_expr(e,i,p) for e in exprs],p)
    result=[]
    for row in H:
        out=defaultdict(int)
        for e,c in zip(exprs,row):
            if c:
                for w,d in e.items():out[w]+=int(c)*int(d)
        result.append({w:(c%p if c%p<=p//2 else c%p-p)
                       for w,c in out.items() if c%p})
    return result


@lru_cache(None)
def staged_basis(sec,parity,p):
    if sec[0]<0:return []
    shapes=defaultdict(list)
    for w in sector_basis(sec,parity):shapes[shape(w)].append(w)
    result=[]
    for words in shapes.values():
        groups=defaultdict(list)
        for w in words:groups[e1_shape(w)].append({w:1})
        e1=[]
        for es in groups.values():e1.extend(kernel(es,1,p))
        if e1:result.extend(kernel(e1,2,p))
    for e in result:
        for i in (1,2):assert not raise_expr(e,i,p)
    return result


def run(sec,parity,p):
    base.basis=lambda s,r:staged_basis(s,r,p)
    base.Frame=BoundedFrame
    result=base.sector_hw(sec,parity,p)
    result['highest_weight_coefficient_field']='F_'+str(p)
    result['highest_weight_algorithm']='E1 then E2, exact finite-field kernels'
    result['source_sha256']['hw_staged.py']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sector',required=True)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--seconds',type=int,default=300);ap.add_argument('--output',required=True)
    a=ap.parse_args();signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('staged HW wall limit')))
    signal.alarm(a.seconds)
    result=run(tuple(map(int,a.sector.split(','))),a.parity,a.prime)
    pathlib.Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
