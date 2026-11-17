"""Staged HW with component-separated exact formal kernels.

Two input columns can interact only if their images share an output word.
Connected components of this incidence relation give independent kernel
blocks. No random compression or rank truncation is used.
"""
import hw_scan as base
import hw_staged as staged
from hw_staged import *
import gc

DENSE_KERNEL=base.formal_closed_kernel


def split_kernel(images,p):
    n=len(images);parent=list(range(n));sizes=[1]*n
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):
        i,j=root(i),root(j)
        if i==j:return
        if sizes[i]<sizes[j]:i,j=j,i
        parent[j]=i;sizes[i]+=sizes[j]
    owners={}
    for i,e in enumerate(images):
        for w,c in e.items():
            if not int(c)%p:continue
            if w in owners:union(i,owners[w])
            else:owners[w]=i
    groups=defaultdict(list)
    for i in range(n):groups[root(i)].append(i)
    del owners;parts=[];rank=0
    for ids in groups.values():
        H,r=DENSE_KERNEL([images[i] for i in ids],p);rank+=r
        if len(H):parts.append((ids,H))
    result=np.zeros((n-rank,n),dtype=np.int64);offset=0
    for ids,H in parts:
        result[np.ix_(range(offset,offset+len(H)),ids)]=H;offset+=len(H)
    assert offset==n-rank
    return result,rank


def run_blocked(sec,parity,p):
    staged.formal_closed_kernel=split_kernel
    original_image_rank=base.faithful_image_rank
    def cleaned_image_rank(*args,**kwargs):
        clear_caches();gc.collect()
        return original_image_rank(*args,**kwargs)
    base.faithful_image_rank=cleaned_image_rank
    def cleaned_kernel(images,prime):
        clear_caches();gc.collect()
        return split_kernel(images,prime)
    base.formal_closed_kernel=cleaned_kernel
    result=staged.run(sec,parity,p)
    result['highest_weight_algorithm']='E1 then E2; formal kernels split into disjoint incidence components'
    result['source_sha256']['hw_blocked.py']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sector',required=True)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--seconds',type=int,default=300);ap.add_argument('--output',required=True)
    a=ap.parse_args();signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('blocked HW wall limit')))
    signal.alarm(a.seconds)
    r=run_blocked(tuple(map(int,a.sector.split(','))),a.parity,a.prime)
    pathlib.Path(a.output).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
