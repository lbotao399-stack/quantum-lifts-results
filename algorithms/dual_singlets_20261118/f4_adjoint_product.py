"""Contract two adjoint words through symmetric products in the Albert representation."""
import argparse
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from f4_adjoint_traces import adjoint_frame
from f4_colour import basis
from f4_first_source import multiply, singlet
from f4_formal_cycles import coefficient_matrix, probe
from f4_joined_boundaries import evaluate, native, parent, vary
from f4_scalar_correction import evaluate_rows
from mixed33_exact import read
from mixed_middle import pin


def projection(prime):
    generators,_=basis()
    trace=np.einsum("iab,jba->ij",generators,generators)
    assert np.array_equal(trace,np.diag(trace.diagonal()))
    dual=np.array([pow(int(g),-1,prime) for g in trace.diagonal()],np.int64)
    products=generators[:,None] @ generators[None,:]
    symmetric=products+products.swapaxes(0,1)
    metric=np.array([1,1,1]+[2]*24,np.int64)
    assert np.array_equal(symmetric*metric[None,None,:,None],
                          symmetric.swapaxes(-1,-2)*metric[None,None,None,:])
    i,j=np.triu_indices(27)
    # A^k_l g^{ll} {T_k,T_l} is a self-adjoint 27x27 matrix. The native
    # trace reads each coefficient against A in reverse matrix-index order.
    covectors=(symmetric[:,:,i,j].transpose(2,1,0)*dual[None,:,None])%prime
    weight=np.array([1 if x==y else 2*int(metric[x])*pow(int(metric[y]),-1,prime)%prime
                     for x,y in zip(i,j)],np.int64)
    rng=np.random.default_rng(886720)
    a,b=rng.integers(0,prime,(2,52,52),dtype=np.int64)
    def direct(x):
        return np.einsum("kl,klij->ij",x*dual%prime,symmetric)%prime
    left,right=direct(a),direct(b)
    av=np.einsum("hkl,lk->h",covectors,a)%prime
    bv=np.einsum("hkl,lk->h",covectors,b)%prime
    assert np.array_equal(av,left[i,j]) and np.array_equal(bv,right[i,j])
    assert int((av*weight%prime)@bv)%prime==int(np.trace(left@right))%prime
    return covectors,weight


def run(plan,previous,work,workers,count,mode,seed_start,batch_size=90):
    work.mkdir(exist_ok=True,parents=True)
    z,old=read(plan),read(previous)
    assert pin(plan)==old["formal_plan_sha256"]
    p=old["prime"]
    lib=native(work,p,"jordan")
    transform=projection(p)
    joined=native(work,p,"adjoint_product",52,len(transform[1]))
    closed=m.library(work,27,p)
    source=mode=="source"
    which="source" if source else "image"
    saved=old["stages"][1]["ranks"][which]
    idx=old["ordinary_rows"]
    ordinary=[idx[i] for i in saved["rows"] if i<len(idx)]
    seeds=[old["stages"][1]["seeds"][i-len(idx)] for i in saved["rows"] if i>=len(idx)]
    degree,dim,salt=(7,666,507260) if source else (8,144,507100)
    dictionary=list(map(m.decode,z["source_dictionary" if source else "target_dictionary"]))
    used=[i for i,w in enumerate(dictionary) if m.nferm(w)==degree]
    words=[dictionary[i] for i in used]
    remap={v:i for i,v in enumerate(used)}
    original=z["sources" if source else "images"]
    rows=[[(remap[j],c) for j,c in original[i]] for i in ordinary]
    if not source:
        lookup={w:i for i,w in enumerate(words)}
        rhs=m.q(multiply(singlet(),4))
        rows.append([(lookup[w],c) for w,c in rhs.items()])
    base,frames=probe(rows,words,{degree:dim},p,closed,workers,salt)
    if not source:
        rhs=base[-1]
        base=base[:-1]
    baseline=[parent(j,1,"jordan",source_only=True) for j in seeds]
    alphabet={l for w in words for t in w for l in t}
    extra,_=evaluate_rows(baseline if source else [vary(x,0,"jordan") for x in baseline],
                          alphabet,frames,p,work,workers,lib,False)
    base=np.vstack((base,extra))
    assert int(nmod_mat(base[:,saved["columns"]].tolist(),p).det())==saved["determinant"]
    frames_with_colours=[adjoint_frame(alphabet,s["points"],s["k"],s["degree"],p,s["seed"]) for s in frames]
    all_values=[base]
    chosen=[]
    graph_counts=[]
    independent=[]
    for start in range(seed_start,seed_start+count,batch_size):
        recipes=[]
        for seed in range(start,min(start+batch_size,seed_start+count)):
            g=parent(seed,1,"adjoint_product",source_only=True)
            if g and m.nferm(next(iter(g)))==7:
                chosen.append(seed)
                recipes.append(g if source else vary(g,0,"adjoint_product"))
        if not recipes:
            continue
        graphs=sorted({g for row in recipes for g in row})
        graph_counts.append(len(graphs))
        lookup={g:i for i,g in enumerate(graphs)}
        co=coefficient_matrix([[(lookup[g],c) for g,c in row.items()] for row in recipes],len(graphs),p)
        blocks=[]
        for f in frames_with_colours:
            values=evaluate(graphs,f,work,workers,joined,"adjoint_product",transform)
            blocks.append(np.asarray(co@values)%p)
        all_values.append(np.column_stack(blocks))
        mat=np.vstack(all_values)
        independent=m.minor(mat,p)
        out=dict(prime=p,plan_sha256=pin(plan),mode=mode,attempted=start+min(batch_size,seed_start+count-start)-seed_start,
                 seeds=chosen,graphs_per_batch=graph_counts,frames=frames,baseline=saved["rank"],minor=independent)
        if not source:
            out["with_required_image"]=m.minor(np.vstack((mat,rhs)),p)
        m.write(work/(mode+".json"),out)
        print(mode,"parents",len(chosen),"rank",independent["rank"],"augmented",
              out.get("with_required_image",{}).get("rank"),flush=True)
        if independent["rank"]==dim:
            break


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    for key in ("plan","previous","work"):
        a.add_argument("--"+key,type=Path,required=True)
    a.add_argument("--workers",type=int,default=12)
    a.add_argument("--count",type=int,default=600)
    a.add_argument("--seed-start",type=int,default=784550)
    a.add_argument("--mode",choices=("image","source"),default="image")
    a.add_argument("--batch-size",type=int,default=90)
    args=a.parse_args()
    assert __debug__
    run(**vars(args))
