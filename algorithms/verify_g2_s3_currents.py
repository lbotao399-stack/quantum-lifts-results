"""Component evaluations of two candidate six-letter singlets in the complete incoming G2 charge sector of S^3."""
import argparse
import contextlib
import io
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from tempfile import TemporaryDirectory

from verify_g2_current_product import FIRST, SECOND
from fortuity.full.letters import Q_letter, letters_up_to


def add(*terms):
    out={}
    for term in terms:
        for monomial,coefficient in term.items():
            out[monomial]=out.get(monomial,0)+coefficient
    return {m:c for m,c in out.items() if c}


def scale(value, scalar):
    return {k:c*scalar for k,c in value.items() if c*scalar}


def exterior(left,right):
    out={}
    for m,c in left.items():
        for n,h in right.items():
            if not m&n:
                # Move every right-hand generator past the larger left-hand ones.
                bits=m; swaps=0
                while bits:
                    bit=bits & -bits; swaps+=(n & (bit-1)).bit_count(); bits-=bit
                out[m|n]=out.get(m|n,0)+(-1)**swaps*c*h
    return {k:v for k,v in out.items() if v}


def diagonal(first,second):
    return [add(scale(first,u),scale(second,v)) for u,v in zip(FIRST,SECOND)]


def trace(left,right):
    # Diagonal entries use a seven-element sequence. Sparse matrices use
    # (row, column) keys so that an off-diagonal product is a true trace.
    if isinstance(left,dict) or isinstance(right,dict):
        lhs=left if isinstance(left,dict) else {(i,i):v for i,v in enumerate(left) if v}
        rhs=right if isinstance(right,dict) else {(i,i):v for i,v in enumerate(right) if v}
        return add(*(exterior(v,rhs.get((j,i),{})) for (i,j),v in lhs.items()))
    assert len(left)==len(right)==7
    return add(*(exterior(a,b) for a,b in zip(left,right)))


AXES=((1,0),(0,1))


def differentiated_current(scalar,d_scalar,gaugino,d_gaugino):
    # The last argument is the differentiated field (not its Taylor coefficient).
    return add(trace(d_scalar,gaugino),trace(scalar,d_gaugino))


def differentiated_jet(jets,flavour,k,h):
    zero=diagonal({},{})
    ak,ah=AXES[k],AXES[h]; second=jets.get((ak[0]+ah[0],ak[1]+ah[1],()),zero)
    # jets store coefficients of ordinary powers of z in the odd superfield.
    factor=1+int(k==h)
    scaled=({ij:scale(v,factor) for ij,v in second.items()} if isinstance(second,dict)
            else [scale(v,factor) for v in second])
    return differentiated_current(jets.get((0,0,(flavour,)),zero),jets.get((*ak,(flavour,)),zero),
                                  jets.get((*ah,()),zero),scaled)


def diagonal_probe(matter):
    zero=diagonal({},{})
    u,v=diagonal({0:1},{}),diagonal({},{0:1})
    jets={(1,0,()):diagonal({1:1},{2:1}),(0,1,()):diagonal({4:1},{8:1}), (0,0,(0,)):u}
    if matter:jets[0,0,(1,2)]=diagonal({16:1},{})
    else:jets.update({(0,0,(1,)):v,(1,0,(2,)):u})
    return jets


def components(jets,matter_sign):
    # matter_sign=-1: literal ordered component comparison; +1: archive
    # Taylor convention, with the increasing-theta orientation restored below.
    assert matter_sign in (-1,1)
    zero=diagonal({},{})
    phi=[jets.get((0,0,(i,)),zero) for i in range(3)]
    gauge=[jets.get((*axis,()),zero) for axis in AXES]
    l=[[trace(phi[i],gauge[k]) for k in range(2)] for i in range(3)]
    g=[]
    for i in range(3):
        m,n=(i+1)%3,(i+2)%3
        fermion=jets.get((0,0,tuple(sorted((m,n)))),zero)
        row=[]
        for k,axis in enumerate(AXES):
            row.append(add(scale(trace(phi[m],jets.get((*axis,(n,)),zero)),Fraction(1,2)),
                           scale(trace(phi[n],jets.get((*axis,(m,)),zero)),Fraction(-1,2)),
                           scale(trace(fermion,gauge[k]),matter_sign*(-1)**i)))
        g.append(row)
    old=scale(add(*(add(exterior(l[i][0],g[i][1]),scale(exterior(l[i][1],g[i][0]),-1))
                    for i in range(3))),Fraction(1,6))
    first=exterior(trace(gauge[0],gauge[1]),old)
    terms=[]
    for order in permutations(range(3)):
        sign=(-1)**sum(order[i]>order[j] for i in range(3) for j in range(i+1,3))
        a,b,c=order
        for i,k,eps1 in ((0,1,1),(1,0,-1)):
            for j,h,eps2 in ((0,1,1),(1,0,-1)):
                dl=differentiated_jet(jets,c,k,h)
                terms.append(scale(exterior(exterior(l[a][i],l[b][j]),dl),sign*eps1*eps2))
    return first,add(*terms)


def probe(matter):
    jets=diagonal_probe(matter); zero=diagonal({},{})
    checks=0
    for letter in letters_up_to(2,2):
        terms=Q_letter(letter)
        for colour in range(7):
            value=add(*(scale(exterior(jets.get(x,zero)[colour],jets.get(y,zero)[colour]),c) for c,x,y in terms))
            assert not value
            checks+=1
    return (*components(jets,-1),checks)


def run():
    rows=[]
    checked=0
    for matter, mask in ((True,31),(False,13)):
        a,b,n=probe(matter)
        rows.append([int(Fraction(v.get(mask,0))) for v in (a,b)])
        assert all(Fraction(v.get(mask,0)).denominator==1 for v in (a,b))
        checked+=n
    det=rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]
    assert rows==[[-64,0],[16,576]] and det==-36864
    return dict(source_doubled_charges=[1,1,6,6,6],length=6,
                Q0_letter_evaluations=checked,cartan_minor=rows,determinant=det)


def classical_ranks():
    from g2_full_classical_screen_v1 import map_ranks, sector
    records=[]
    expected=[(15,6,9),(48,28,20),(30,17,11),(3,0,3)]
    with TemporaryDirectory() as scratch:
        for t in range(4):
            weight=(t,t,3-t,3-t,3-t)
            length=8-t
            assert sector(weight,length)==(1,1,6,6,6)
            for prime in (1000003,1000033):
                folder=Path(scratch)/str(prime)/str(t)
                folder.mkdir(parents=True)
                with contextlib.redirect_stdout(io.StringIO()):
                    data=map_ranks(weight,length,"G2",prime,folder)
                dim,incoming,outgoing=expected[t]
                assert (data["dimensions"]["0"],data["incoming_rank_lower_bound"],
                        data["outgoing_rank_lower_bound"]) == (dim,incoming,outgoing)
                records.append(dict(prime=prime,weight=weight,length=length,dimension=dim,
                                    incoming_rank_lower=incoming,outgoing_rank_lower=outgoing,
                                    cohomology_upper=dim-incoming-outgoing))
    return records


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path)
    p.add_argument("--recompute-ranks",action="store_true")
    args=p.parse_args()
    result=run()
    assert result["determinant"]
    if args.recompute_ranks:
        result["classical_maps"]=classical_ranks()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
