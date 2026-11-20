"""Weyl-invariant Cartan jet dimensions and restrictions of the formal F4 source singlets."""
import argparse
import math
from collections import Counter
from functools import lru_cache
from itertools import permutations
from pathlib import Path

import numpy as np
from sympy import Matrix

import mixed_incoming as m
from charge_maps import content_patterns
from f4_colour import basis
from f4_singlet_characters import CARTAN
from mixed33_exact import key, read
from mixed11_products import FAMILIES, PROBES, contract, fields, partitions, sign, wedge
import mixed11_products as grav


@lru_cache(None)
def weyl_spectra():
    eye=np.eye(4,dtype=np.int64)
    simple=[]
    for i,row in enumerate(CARTAN["F4"]):
        s=eye.copy()
        s[i]-=np.asarray(row)
        assert np.array_equal(s@s,eye)
        simple.append(s)
    seen={eye.tobytes()}
    todo=[eye]
    spectra=Counter()
    while todo:
        a=todo.pop()
        x=a
        traces=[]
        for _ in range(4):
            traces.append(int(x.trace()));x=x@a
        spectra[tuple(traces)]+=1
        for s in simple:
            b=a@s
            key=b.tobytes()
            if key not in seen:
                seen.add(key);todo.append(b)
    assert len(seen)==sum(spectra.values())==1152 and (-eye).tobytes() in seen
    return spectra


@lru_cache(None)
def characters(traces):
    exterior=[1]
    for n in range(1,5):
        num=sum((-1)**(j-1)*traces[j-1]*exterior[n-j] for j in range(1,n+1))
        assert num%n==0
        exterior.append(num//n)
    symmetric=[1]
    for n in range(1,19):
        symmetric.append(sum((-1)**(j-1)*exterior[j]*symmetric[n-j] for j in range(1,min(4,n)+1)))
    return symmetric,exterior


@lru_cache(None)
def multiplicity(bosons,fermions):
    if any(x>4 for x in fermions):
        return 0
    total=0
    for traces,size in weyl_spectra().items():
        s,e=characters(traces)
        total+=size*math.prod(s[i] for i in bosons)*math.prod(e[j] for j in fermions)
    assert total%1152==0 and total>=0
    return total//1152


def dimensions(weight,length):
    result=Counter()
    for spin,ss in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            shifted=(weight[0]+spin,weight[1]-spin)+tuple(weight[2+i]+perm[i]-i for i in range(3))
            for bosons,fermions,num in content_patterns(shifted,length):
                result[sum(fermions)]+=ss*sign(perm)*num*multiplicity(bosons,fermions)
    assert all(n>=0 for n in result.values())
    return {k:v for k,v in sorted(result.items()) if v}


def commuting(matrices):
    eligible=[set(j for j in range(i+1,len(matrices)) if np.array_equal(x@matrices[j],matrices[j]@x))
              for i,x in enumerate(matrices)]
    for i in range(len(matrices)):
        for j in sorted(eligible[i]):
            for k in sorted(eligible[i]&eligible[j]):
                remaining=eligible[i]&eligible[j]&eligible[k]
                if remaining:
                    return (i,j,k,min(remaining))
    raise ValueError("Cartan subalgebra requires a different colour basis")


def probes(probe,cartan):
    jet,mask=fields(probe)
    @lru_cache(None)
    def trace_value(tr):
        if any(l not in jet for l in tr):
            return 0,0
        colour=np.eye(27,dtype=np.int64)
        term=(1,0)
        for l in tr:
            label,c,ext=jet[l]
            term=wedge(term,(c,ext))
            if not term[0]:
                return 0,0
            colour=colour@cartan[label]
        return int(colour.trace())*term[0],term[1]
    def trace_for_jet(tr,actual,group):
        assert actual==jet
        return trace_value(tr)
    return jet,mask,trace_for_jet


def run(formal):
    table=[]
    for power in (4,5):
        for t in range(power+1):
            w=(t,t)+(power-t,)*3
            n=sum(w)-1
            table.append(dict(power=power,t=t,weight=w,length=n,dimensions=dimensions(w,n)))
    assert [multiplicity((i,),()) for i in range(1,7)]==[0,1,0,1,0,2]
    archive=read(formal)
    assert archive["mu"]==[1,1,3,3,3]
    mats,_=basis();indices=commuting(mats)
    cartan=np.stack([mats[i]//int(np.gcd.reduce(np.abs(mats[i]).ravel())) for i in indices])
    assert int(np.max(abs(cartan).sum(axis=2)))<=2 and math.factorial(3)*27**5*2**10<2**63
    # Verify the ordered Q0 letter rules on every Cartan assignment.
    from fortuity.full.letters import letters_up_to,Q_letter
    classes=[{key(w):c for w,c in row} for block in archive["blocks"] for row in block["representatives"]]
    assert len(classes)==7 and all(not m.q(e) for e in classes)
    crows=[];tensor_rows=[];checks=0
    for sample in PROBES:
        jet,mask,evaluate_trace=probes(sample,cartan)
        for l in letters_up_to(1,1):
            out={}
            for c,a,b in Q_letter(l):
                if a in jet and b in jet:
                    i,v,s=jet[a];j,w,t=jet[b]
                    val,ext=wedge((v,s),(w,t))
                    out[ext]=out.get(ext,0)+c*val*(cartan[i]@cartan[j])
            assert all(not np.any(v) for v in out.values())
            checks+=1
        crows.append([int(grav.expression_value(row,jet,"F4",evaluate_trace).get(mask,0)) for row in classes])
        tensor_rows.append([grav.tensor(f,sample,"F4",evaluate_trace) for f in FAMILIES])
    choices=list(partitions(tuple(range(9))))
    product_rows=[]
    for i in range(len(FAMILIES)):
        for choice in choices:
            col=[contract(row[i],choice) for row in tensor_rows]
            if any(col):
                product_rows.append((i,choice,col))
    matrix=Matrix([row[2] for row in product_rows]).T
    piv=matrix.rref()[1]
    selected=matrix[:,list(piv)]
    row_piv=selected.T.rref()[1]
    return dict(weyl_order=1152,spectra=[[list(k),v] for k,v in sorted(weyl_spectra().items())],
                cartan_dimensions=table,cartan_generator_indices=indices,variations_checked=checks,
                classical_probe_matrix=crows,classical_rank=Matrix(crows).rank(),product_rank=len(piv),
                product_columns=[dict(family=product_rows[i][0],triples=product_rows[i][1]) for i in piv],
                product_probe_matrix=[[str(x) for x in row] for row in selected.tolist()],
                product_minor_rows=list(row_piv),product_determinant=str(selected[list(row_piv),:].det()))


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--formal",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args();assert __debug__
    z=run(args.formal)
    m.write(args.out,z)
    print([(r["power"],r["t"],r["dimensions"]) for r in z["cartan_dimensions"]],
          z["product_rank"],z["product_determinant"],flush=True)
