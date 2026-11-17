"""Highest-weight BMN cohomology with exact character dimension controls.

The coefficient field for differential ranks is the reported prime field.
This program does not classify positive cohomology as fortuitous.
"""
from highest import *
from functools import lru_cache
from itertools import permutations
import hashlib


@lru_cache(None)
def exact_hw(sec,parity):
    if sec[0]<0:return {}
    out=defaultdict(int);rho=(2,1,0)
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[k] for i in range(3) for k in range(i+1,3))
        q=tuple(sec[i+1]+rho[i]-rho[perm[i]] for i in range(3))
        if min(q)<0:continue
        for d,m in expected((sec[0],)+q,parity).items():out[d]+=sign*m
    assert all(m>=0 for m in out.values()),(sec,parity,out)
    return {d:m for d,m in out.items() if m}


@lru_cache(None)
def basis(sec,parity):
    if sec[0]<0:return []
    return hw_basis(sec,parity)


def faithful_image_rank(sec,parity,images,prime,seed,include_basis=False):
    exact=exact_hw(sec,parity)
    base=basis(sec,parity);exprs=base+images
    all_words=sorted({w for e in exprs for w in e})
    degrees={w:sum(count(w)[3:6]) for w in all_words}
    all_blocks=[];basis_blocks=[];controls={}
    for d,dim in sorted(exact.items()):
        words=[w for w in all_words if degrees[w]==d];idx={w:i for i,w in enumerate(words)}
        rr=[];cc=[];vv=[]
        for i,e in enumerate(exprs):
            for w,c in e.items():
                if w in idx:rr.append(i);cc.append(idx[w]);vv.append(int(c)%prime)
        A=csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(exprs),len(words)))
        blocks=[];used=0;k=0 if d==0 else d+2;perpoint=comb(k,d)
        for attempt in range(4):
            target=dim+8 if attempt==0 else 2*used
            while used<target:
                npts=min(max(1,32//perpoint),8,ceil((target-used)/perpoint))
                fr=Frame(d,npts,seed+100003*d+7919*used,prime,k=k)
                X=fr.values(words)
                blocks.append((np.asarray(A@X)%prime).astype(np.int32))
                used+=fr.ncols
                del fr,X
            E=np.concatenate(blocks,axis=1).astype(np.int64)
            r=rank_mod(E[:len(base)],prime)
            assert r<=dim,(sec,parity,d,r,dim)
            if r==dim:break
        assert r==dim,('highest-weight evaluation incomplete',sec,parity,d,r,dim)
        # Every image is in the highest-weight target. Check it does not
        # enlarge the dimension certified from the character.
        assert rank_mod(E,prime)==dim
        controls[d]={'dimension':dim,'evaluation_rank':r,'columns':used}
        all_blocks.append(E[len(base):])
        if include_basis:basis_blocks.append(E[:len(base)])
    rank=rank_mod(np.concatenate(all_blocks,axis=1),prime) if all_blocks and images else 0
    if include_basis:
        return rank,controls,np.concatenate(basis_blocks,axis=1) if basis_blocks else np.zeros((len(base),0),dtype=np.int64)
    return rank,controls


def formal_closed_kernel(images,prime):
    """Kernel before finite-rank identities, in the canonical trace algebra."""
    cols={w:k for k,w in enumerate(sorted({w for e in images for w in e}))}
    M=nmod_mat(len(cols),len(images),prime)
    for i,e in enumerate(images):
        for w,c in e.items():M[cols[w],i]=int(c)%prime
    R,r=M.rref()
    piv=[next(i for i in range(len(images)) if R[k,i]) for k in range(r)]
    free=[i for i in range(len(images)) if i not in piv]
    H=np.zeros((len(free),len(images)),dtype=np.int64)
    for k,i in enumerate(free):
        H[k,i]=1
        for j,pivot in enumerate(piv):H[k,pivot]=-int(R[j,i])%prime
    return H,int(r)


def sector_hw(sec,parity,prime=1000003,seed=2201):
    started=time.monotonic();exact=exact_hw(sec,parity);dim=sum(exact.values())
    result={'sector':sec,'reflection_parity':(-1)**parity,'prime':prime,
            'highest_weight_dimension':dim,'exact_dimensions_by_fermions':exact}
    if dim:
        B=basis(sec,parity)
        incoming=basis((sec[0]+1,)+sec[1:],parity)
        rank_in,source_controls,E=faithful_image_rank(sec,parity,[action_expr(e) for e in incoming],prime,seed,True)
        outgoing=[action_expr(e) for e in B]
        rank_out,target_controls=faithful_image_rank((sec[0]-1,)+sec[1:],parity,outgoing,prime,seed+200003)
        result.update(formal_source_highest_weights=len(B),formal_incoming_highest_weights=len(incoming),
                      rank_in=rank_in,rank_out=rank_out,cohomology_multiplicity=dim-rank_in-rank_out,
                      source_controls=source_controls,target_controls=target_controls)
        if parity==0:
            H,formal_rank=formal_closed_kernel(outgoing,prime)
            closed_rank=rank_mod(H@E%prime,prime)
            result.update(universal_Q_rank=formal_rank,universal_closed_dimension=len(H),
                          universal_closed_evaluated_rank=closed_rank,
                          rank_dependent_trace_multiplicity_mod_prime=dim-rank_out-closed_rank)
            assert 0<=result['rank_dependent_trace_multiplicity_mod_prime']<=result['cohomology_multiplicity']
    else:result.update(rank_in=0,rank_out=0,cohomology_multiplicity=0)
    assert result['cohomology_multiplicity']>=0
    result.update(seconds=time.monotonic()-started,maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  classification='Ranks over the stated finite field. The rank-dependent trace quotient can contain Pfaffian products; it is not a count of new fortuitous cores.',
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'highest.py',ROOT/'hw_scan.py')})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--sector',required=True);parser.add_argument('--parity',type=int,default=0)
    parser.add_argument('--prime',type=int,default=1000003);parser.add_argument('--seconds',type=int,default=600)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('highest-weight sector limit')))
    signal.alarm(args.seconds)
    result=sector_hw(tuple(map(int,args.sector.split(','))),args.parity,args.prime)
    pathlib.Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
