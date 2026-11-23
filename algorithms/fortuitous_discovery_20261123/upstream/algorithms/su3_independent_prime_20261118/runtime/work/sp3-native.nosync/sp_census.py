"""Isolated symplectic BMN census using the existing formal trace kernels."""
import sp_environment
import sys
sys.path.insert(0,str(sp_environment.WORK/'cluster-next-v2.nosync/work/so8-protected'))
from sp_dimensions import multiplicity
import hw_scan as base
import hw_staged as staged
import hw_blocked as blocked
import hw_streamed_kernel as streamed
from fortuity.hwcoh import evalg
from fortuity.bmn.words import contents_of_sector
import argparse,hashlib,json,pathlib,time,gc
from collections import defaultdict


def configure(rank,prime):
    def expected(sector,parity):
        assert parity==0
        result=defaultdict(int)
        for content in contents_of_sector(sector):
            result[sum(content[3:6])]+=multiplicity(content,rank)
        return dict(result)

    class Frame(evalg.Frame):
        def __init__(self,degree,npts,seed,prime,k=None):
            super().__init__(rank,npts,degree if k is None else k,[degree],seed=seed,p=prime,group='Sp')
        def values(self,words):
            assert all(not paths for paths,traces in words)
            return self.eval_words([traces for paths,traces in words],dtype=base.np.int64)

    base.expected=expected
    base.exact_hw.cache_clear()
    def basis(sector,parity):
        result=staged.staged_basis(sector,parity,prime)
        print(json.dumps(dict(stage='basis',sector=sector,rows=len(result))),flush=True)
        return result
    base.basis=basis
    base.Frame=Frame
    blocked.DENSE_KERNEL=lambda images,prime:streamed.bounded_kernel(images,prime,chunk_rows=128,dense_entries=1000000,verbose=True)
    staged.formal_closed_kernel=blocked.split_kernel
    base.formal_closed_kernel=blocked.split_kernel
    original=base.faithful_image_rank
    def faithful(*args,**kwargs):
        base.clear_caches();gc.collect()
        return original(*args,**kwargs)
    base.faithful_image_rank=faithful


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rank',type=int,required=True)
    parser.add_argument('--sector',required=True)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    configure(args.rank,args.prime)
    sector=tuple(map(int,args.sector.split(',')))
    result=base.sector_hw(sector,0,args.prime)
    result.pop('reflection_parity',None)
    result.update(group='Sp(%d)'%args.rank,matrix_size=2*args.rank,
                  highest_weight_coefficient_field='F_%d'%args.prime,
                  classification='Finite-field cohomology and quotient by universally closed formal traces. Positive residuals require known-family subtraction and quantum tests.',
                  source_sha256={**result['source_sha256'],**{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (pathlib.Path(__file__),pathlib.Path(__file__).parent/'sp_dimensions.py')}})
    pathlib.Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
