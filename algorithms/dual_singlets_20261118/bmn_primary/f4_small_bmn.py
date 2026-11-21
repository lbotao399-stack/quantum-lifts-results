"""Classical image minors in the exceptional level-18 singlet and its lower BMN sectors."""
import argparse
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import mixed_incoming as m
from bmn_primary_index import cochains
from f4_formal_cycles import probe


def compute(weight, work, prime, workers):
    weight=tuple(weight)
    mu=(0,0,*weight)
    levels=cochains(sum(weight),"F4")[1]
    sizes={y:d[weight][0] for y,d in levels.items() if weight in d}
    library=m.library(work,27,prime)
    previous=[]
    data=[]
    for y in range(min(sizes),max(sizes)+1):
        blocks=m.words_by_shape(y,mu)
        sources=[e for _,block in sorted(blocks.items()) for e in m.singlets(block)]
        images=[m.q(e) for e in sources]
        if any(m.q(e) for e in images):
            raise ArithmeticError("Classical action did not square to zero")
        closed=[src for src,img in zip(sources,images) if not img]
        families=previous+sources
        words,coefficients=m.flatten(families)
        dim=sizes.get(y,0)
        if words:
            # The overall dimension is a valid sampling ceiling for each fermion degree.
            dims={m.nferm(w):dim for w in words}
            salt=298771+sum((i+1)*w for i,w in enumerate(weight))*757+y*99
            evaluations,frames=probe(coefficients,words,dims,prime,library,workers,salt)
            zero={"rank":0,"rows":[],"columns":[],"determinant":1}
            incoming=m.minor(evaluations[:len(previous)],prime) if previous else zero
            source=m.minor(evaluations[len(previous):],prime) if sources else zero
            closed_ids=[i for i,e in enumerate(images) if not e]
            combined=list(range(len(previous)))+[len(previous)+i for i in closed_ids]
            closed_minor=m.minor(evaluations[combined],prime) if combined else zero
            if max(incoming["rank"],source["rank"],closed_minor["rank"])>dim:
                raise ArithmeticError("An image exceeded the independent character dimension")
        else:
            incoming=source=closed_minor={"rank":0}
            frames=[]
            closed_ids=[]
        data.append(dict(length=y,dimension=dim,formal_sources=len(sources),formal_closed=len(closed),
                         formal_closed_indices=closed_ids,incoming=incoming,source=source,
                         incoming_and_formally_closed=closed_minor,frames=frames))
        print(weight,y,dim,incoming["rank"],source["rank"],closed_minor["rank"],flush=True)
        previous=images
    for i,row in enumerate(data):
        outgoing=data[i+1]["incoming"]["rank"] if i+1<len(data) else 0
        row["cohomology_upper"]=row["dimension"]-row["incoming"]["rank"]-outgoing
        if row["cohomology_upper"]<0:
            raise ArithmeticError("The two adjacent rank bounds are incompatible")
    return dict(group="F4",weight=weight,prime=prime,rows=data)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weight",nargs=3,type=int,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prime",type=int,default=1000003)
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    a.work.mkdir(exist_ok=True,parents=True)
    result=compute(a.weight,a.work,a.prime,a.workers)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
