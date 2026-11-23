"""Highest-weight screening over one finite field, without rational basis blowup.

Completeness is still controlled by the same independent exact characters.
No rational cohomology certificate is inferred from this screen.
"""
import hw_scan as base
from hw_scan import *
from collections import OrderedDict


class BoundedFrame(Frame):
    """Retain at most 64 MiB of matrix-valued word prefixes per sample frame."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.matrices=OrderedDict();self.matrix_bytes=0

    def matrix(self,w):
        if w in self.matrices:
            self.matrices.move_to_end(w);return self.matrices[w]
        if len(w)==1:d,v=self.L[w[0]]
        else:
            v,d=self.matrix(w[:-1]);v,d=self.mat_mul_letter(v,d,w[-1])
        size=v.nbytes
        while self.matrices and self.matrix_bytes+size>64*1024**2:
            _,(old,_) = self.matrices.popitem(last=False);self.matrix_bytes-=old.nbytes
        if size<=64*1024**2:
            self.matrices[w]=(v,d);self.matrix_bytes+=size
        return v,d


def modular_hw_basis(sec,parity,p):
    groups=defaultdict(list)
    for w in sector_basis(sec,parity):groups[shape(w)].append(w)
    result=[]
    for words in groups.values():
        images=[{(i,v):c for i in (1,2) for v,c in raising(w,i).items()} for w in words]
        cols={v:k for k,v in enumerate(sorted({v for im in images for v in im}))}
        if not cols:
            result.extend({w:1} for w in words);continue
        M=nmod_mat(len(cols),len(words),p)
        for j,im in enumerate(images):
            for v,c in im.items():M[cols[v],j]=int(c)%p
        R,r=M.rref();piv=[next(i for i in range(len(words)) if R[k,i]) for k in range(r)]
        pivset=set(piv)
        for free in range(len(words)):
            if free in pivset:continue
            e={words[free]:1}
            for k,pivot in enumerate(piv):
                c=-int(R[k,free])%p
                if c:e[words[pivot]]=c if c<=p//2 else c-p
            for i in (1,2):
                check=defaultdict(int)
                for w,c in e.items():
                    for v,d in raising(w,i).items():check[v]+=c*d
                assert all(c%p==0 for c in check.values())
            result.append(e)
    return result


def run(sec,parity,p):
    @lru_cache(None)
    def modular_basis(s,r):
        return modular_hw_basis(s,r,p) if s[0]>=0 else []
    base.basis=modular_basis
    base.Frame=BoundedFrame
    r=base.sector_hw(sec,parity,p)
    r['highest_weight_coefficient_field']='F_'+str(p)
    r['source_sha256']['hw_modular.py']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    return r


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sector',required=True)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--seconds',type=int,default=120);ap.add_argument('--output',required=True)
    a=ap.parse_args();signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('wall limit')))
    signal.alarm(a.seconds)
    r=run(tuple(map(int,a.sector.split(','))),a.parity,a.prime)
    pathlib.Path(a.output).write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r),flush=True)
