#!/usr/bin/env python3
"""Rebuild both SO13 fortuitous relation quotients over Q, using stdlib only.

The optional cochain query returns an actual coefficient of the explicitly
specified odd-trace primitive. It does not merely name an unspecified lift.
"""
import argparse,hashlib,json,time
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import pfaffian_scalar_candidates as p
from so13_character_certificate import flavour,predecessor_shapes,transpose,certificate as character_certificate

HERE=Path(__file__).resolve().parent

def raising(poly,i,j):
    """E_ij Phi_j=Phi_i; E_ij psi_i=-psi_j for i != j."""
    out=Counter()
    for m,c in poly.items():
        for a in range(3):
            for old,new,s in [(3*j+a,3*i+a,1),(9+3*i+a,9+3*j+a,-1)]:
                if not m[old]:continue
                n=list(m);n[old]-=1;n[new]+=1
                out[tuple(n)]+=c*s*m[old]
    return {m:c for m,c in out.items() if c}

def physical_weight(m):
    return tuple(sum(m[3*i:3*i+3])+sum(sum(m[9+3*j:12+3*j]) for j in range(3) if j!=i) for i in range(3))

def flavour_action(expr,i,j):
    out=Counter()
    for w,c in expr.items():
        for k,t in enumerate(w):
            for a,l in enumerate(t):
                if l==j:new,s=i,1
                elif l==3+i:new,s=3+j,-1
                else:continue
                orientation,u=p.canonical(w[:k]+(t[:a]+(new,)+t[a+1:],)+w[k+1:])
                if orientation:out[u]+=c*s*orientation
    return {u:c for u,c in out.items() if c}

def raw_primitive(w,k):
    """The repository's adjacent-swap primitive on the specified odd trace."""
    assert all(l<3 for t in w for l in t)
    t=list(w[k]);out=Counter()
    for end in range(len(t)-1,0,-1):
        for i in range(end):
            a,b=t[i:i+2]
            if a!=b:
                j=next(j for j,pair in enumerate(p.PAIRS) if {a,b}==set(pair))
                orientation=1 if (a,b)==p.PAIRS[j] else -1
                s,u=p.canonical(w[:k]+(tuple(t[:i]+[3+j]+t[i+2:]),)+w[k+1:])
                if s:out[u]+=Fraction(s*orientation,2)
            t[i],t[i+1]=t[i+1],t[i]
    return {u:c for u,c in out.items() if c}

@lru_cache(None)
def primitive(w):
    """Flavour-equivariant cyclic/reversal average of the ordered h_B.

    Average over odd trace factors, all their cyclic rotations, and the
    antisymmetric reversal representative. This is independent of trace
    presentation and labels, so it intertwines the full flavour action.
    The underlying adjacent-swap primitive and physical Q0 are unchanged.
    """
    positions=[i for i,t in enumerate(w) if len(t)%2];assert positions
    out=Counter()
    for k in positions:
        t=w[k]
        factor=Fraction(1,2*len(t)*len(positions))
        for v,s in [(t,1),(t[::-1],-1)]:
            for j in range(len(t)):
                raw=w[:k]+(v[j:]+v[:j],)+w[k+1:]
                for u,c in raw_primitive(raw,k).items():out[u]+=factor*s*c
    return {u:c for u,c in out.items() if c}

def scalar_preimages(u):
    """Inverse fixed reversal steps; O(L^3) candidates rather than all words."""
    positions=[(k,j,l) for k,t in enumerate(u) for j,l in enumerate(t) if l>=3]
    assert len(positions)==1
    k,_,l=positions[0];t=u[k]
    if len(t)%2:return ()
    orbit={v[i:]+v[:i] for v in (t,t[::-1]) for i in range(len(t))}
    size=len(t)+1;permutation=list(range(size));candidates=set()
    for end in range(size-1,0,-1):
        for i in range(end):
            for raw in orbit:
                if raw[i]!=l:continue
                a,b=p.PAIRS[l-3]
                for pair in [(a,b),(b,a)]:
                    expanded=raw[:i]+pair+raw[i+1:]
                    original=[None]*size
                    for j,v in enumerate(expanded):original[permutation[j]]=v
                    s,w=p.canonical(u[:k]+(tuple(original),)+u[k+1:])
                    if s:candidates.add(w)
            permutation[i],permutation[i+1]=permutation[i+1],permutation[i]
    return tuple(sorted(candidates))

class PrimitiveCoefficient:
    def __init__(self,terminal):self.terminal=terminal
    @lru_cache(None)
    def __call__(self,u):
        s,u=p.canonical(u)
        if not s:return Fraction()
        return s*sum(primitive(w).get(u,0)*self.terminal(w) for w in scalar_preimages(u))

def q_primitive(w,terminal):
    lifted=PrimitiveCoefficient(terminal)
    return sum(c*lifted(u) for u,c in p.predecessors(w).items())

def record(shape):
    start=time.monotonic();sources,terminal=p.marker_polynomials(shape)
    source_polys=sources;sources=list(map(p.Coefficient,sources));terminal=p.Coefficient(terminal)
    for poly in source_polys+[terminal.raw]:
        assert {physical_weight(m) for m in poly}=={shape}
        assert not raising(poly,0,1) and not raising(poly,1,2)
    assert all({sum(m[:9]) for m in q}=={13} and {sum(m[9:]) for m in q}=={1} for q in source_polys)
    assert {sum(m) for m in terminal.raw}=={15}
    path=HERE/('so13_'+''.join(map(str,shape))+'_probe.json')
    old=json.loads(path.read_text());words=list(map(p.decode,old['basis_probe_words']))
    matrix=[[p.image(q,w) for q in sources]+[terminal(w)] for w in words]
    assert [[str(c) for c in r] for r in matrix]==old['exact_probe_matrix']
    assert p.determinant(matrix)==Fraction(old['exact_determinant'])
    srcdim=sum(flavour(s,1) for s in predecessor_shapes(shape))
    tgtdim=flavour(shape)
    assert len(sources)==srcdim and len(words)==tgtdim==srcdim+1
    diagram=transpose(shape)
    source_rank=p.matrix_rank([r[:-1] for r in matrix]);assert source_rank==srcdim
    dual=p.solve(list(zip(*matrix)),[0]*srcdim+[1])
    assert list(map(str,dual))==old['dual_coefficients']
    scores=[sum(dual[i]*matrix[i][j] for i in range(len(words))) for j in range(tgtdim)]
    assert scores==[0]*srcdim+[1]
    primitive_checks=[]
    for w in words:
        actual=Counter()
        for u,c in primitive(w).items():
            for v,d in p.differential(u).items():actual[v]+=c*d
        assert {v:d for v,d in actual.items() if d}=={w:Fraction(1)}
        for i,j in [(0,1),(1,2),(1,0),(2,1)]:
            left=flavour_action(primitive(w),i,j);right=Counter()
            for v,c in flavour_action({w:1},i,j).items():
                for u,d in primitive(v).items():right[u]+=c*d
            assert left=={u:c for u,c in right.items() if c}
        closure=q_primitive(w,terminal);assert closure==terminal(w)
        primitive_checks.append(dict(probe=p.encode(w),Q0_h_Pfaffian_coefficient=str(closure),Pfaffian_coefficient=str(terminal(w))))
    return dict(passed=True,group='SO(13)',scalar_highest_weight=shape,j=30,source_length=14,
                source_field_content={'scalars':13,'matter_fermions':1},Taylor_diagram=diagram,
                exact_character_source_relation_dimension=srcdim,exact_character_terminal_relation_dimension=tgtdim,
                predecessor_scalar_Pieri_shapes=predecessor_shapes(shape),character_previous_relation_dimension=0,
                source_highest_weight_checks=True,terminal_highest_weight_checks=True,
                primitive_flavour_equivariance_checks=True,
                exact_Q0_relation_rank=source_rank,exact_fortuitous_quotient_dimension=tgtdim-source_rank,
                ordered_probe_words=old['basis_probe_words'],
                matrix_columns=['Q0_source_'+str(i+1) for i in range(srcdim)]+['terminal_Pfaffian'],
                exact_probe_matrix=[[str(c) for c in r] for r in matrix],exact_determinant=str(p.determinant(matrix)),
                fortuity_dual_coefficients=list(map(str,dual)),dual_pairings=list(map(str,scores)),
                actual_representative='GL3-equivariant cyclic/reversal average primitive() applied coefficientwise to the explicitly specified terminal Pfaffian',
                cochain_coefficient_oracle='PrimitiveCoefficient(Coefficient(marker_polynomials(shape)[1]))',
                primitive_Q0_checks=primitive_checks,
                physical_scope='Full zero-dotted-weight classical Q0 cohomology modulo stable-image cohomology. No quantum survival claim.',
                seconds=time.monotonic()-start)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=HERE/'SO13_EXACT_CERTIFICATES.json')
    ap.add_argument('--query-shape',type=int,nargs=3);ap.add_argument('--query-word')
    ap.add_argument('--native-sorted-theta',action='store_true',help='Use a=Psi_BC,b=Psi_AC,c=Psi_AB; cyclic psi_B=-Psi_AC.')
    a=ap.parse_args()
    if a.query_word:
        assert a.query_shape
        _,terminal=p.marker_polynomials(tuple(a.query_shape))
        value=PrimitiveCoefficient(p.Coefficient(terminal))(p.decode(a.query_word))
        if a.native_sorted_theta:value*=(-1)**a.query_word.count('b')
        print(json.dumps({'shape':a.query_shape,'word':a.query_word,'native_sorted_theta':a.native_sorted_theta,'cochain_coefficient':str(value)}));return
    rows=[record(shape) for shape in [(5,5,5),(7,5,3)]]
    data={'passed':True,'field':'Q','characters':character_certificate(),'records':rows,'code_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['pfaffian_scalar_candidates.py','so13_exact_certificate.py','so13_character_certificate.py']}}
    a.out.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
