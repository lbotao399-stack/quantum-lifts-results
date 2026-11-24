"""Full stable-kernel and finite-relation test, independent of ring primitivity."""
import argparse
from pathlib import Path
import json,pickle,time
from maps import stable_map,sparse_combine,finite_coordinates
from finite_frame import PRIME


def rank_receipt(rows,prime=PRIME,label="matrix"):
    order=sorted(range(len(rows)),key=lambda j:(len(rows[j]),j))
    echelon={};selected=[];columns=[];diagonal=[]
    start=time.monotonic()
    for count,j in enumerate(order):
        row={k:c.numerator*pow(c.denominator,-1,prime)%prime for k,c in rows[j].items()}
        row={k:c for k,c in row.items() if c}
        while row:
            leading=min(row)
            if leading not in echelon:
                diagonal.append(row[leading]);inverse=pow(row[leading],-1,prime)
                echelon[leading]={k:c*inverse%prime for k,c in row.items()}
                selected.append(j);columns.append(leading);break
            coefficient=row[leading]
            for k,c in echelon[leading].items():
                value=(row.get(k,0)-coefficient*c)%prime
                if value:row[k]=value
                else:row.pop(k,None)
        if count%500==0:print(json.dumps(dict(stage=label,done=count,total=len(order),rank=len(selected),seconds=time.monotonic()-start)),flush=True)
    return dict(prime=prime,rank=len(selected),rows=selected,columns=columns,diagonal=diagonal,
                replay="Sparse row elimination in increasing original support size; smallest-column pivot.")


def build(work):
    incoming=pickle.loads((work/"map6.pkl").read_bytes());outgoing=pickle.loads((work/"map7.pkl").read_bytes())
    source=pickle.loads((work/"stable7.pkl").read_bytes());products=pickle.loads((work/"products.pkl").read_bytes())
    stable_products=stable_map(dict(length="products",sources=products["complete"],images=products["complete"]),source)
    for row in incoming["stable"]+stable_products:assert not sparse_combine(row,outgoing["stable"])
    lower_out=rank_receipt(outgoing["stable"],label="stable_outgoing")
    lower_kernel=rank_receipt(incoming["stable"]+stable_products,label="stable_closed_span")
    assert lower_out["rank"]+lower_kernel["rank"]==4991,(lower_out["rank"],lower_kernel["rank"])
    source_basis,source_coordinate=finite_coordinates(work/"frame7",4991)
    finite_products=[sparse_combine(row,source_coordinate) for row in stable_products]
    lower_finite=rank_receipt(incoming["finite"]+finite_products,label="stable_closed_finite_image")
    assert lower_finite["rank"]==858
    # The complete character-certified source kernel has 3095 relations.
    # Restricting all stable cycles gives exactly 858 finite cycles.
    # Rank-nullity therefore determines the ENTIRE relation-image quotient.
    stable_kernel=lower_kernel["rank"]
    relation_kernel=stable_kernel-858
    relation_image=3095-relation_kernel
    closed_lift_dimension=3095+859
    target_relations_hit=closed_lift_dimension-stable_kernel
    assert target_relations_hit-relation_image==1
    result=dict(stable_outgoing=lower_out,stable_closed_span=lower_kernel,stable_closed_finite_image=lower_finite,
                stable_source=4991,finite_source=1896,complete_source_relations=3095,
                stable_target=6609,finite_target=1588,complete_target_relations=5021,
                stable_kernel=stable_kernel,source_relations_closed=relation_kernel,
                image_of_all_source_relations=relation_image,
                target_relations_in_stable_image=target_relations_hit,fortuitous_quotient=1,
                stable_products=stable_products)
    (work/"relation_test.pkl").write_bytes(pickle.dumps(result,protocol=5))
    summary={k:v for k,v in result.items() if isinstance(v,int)}
    (work/"relation_test.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path);args=parser.parse_args();build(args.work)
