#!/usr/bin/env python3
"""Exact minimal-degree scalar Pfaffian classes for arbitrary odd SO colours.

Run the specified sectors only. All coefficients retain the native scalar
Pfaffian/ordinary-differentiation conventions already certified for SO13.
"""
import argparse,hashlib,json,random,time
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import pfaffian_scalar_candidates as p
import so13_exact_certificate as h
from so13_character_certificate import character,partitions,centralizer,predecessor_shapes,transpose

@lru_cache(None)
def relation_trace(cycles,colours):
    degree=sum(cycles)
    if degree<=colours:return 0
    if degree==colours+1:shapes=((degree,degree),)
    elif degree==colours+2:shapes=((degree,degree),(degree-1,degree-1,1,1))
    else:raise ValueError((degree,colours))
    types=sorted(Counter(cycles).items());value=0
    for crossed in product(*(range(m+1) for _,m in types)):
        vc=[];coefficient=1
        for (k,m),a in zip(types,crossed):
            vc += [k]*(2*(m-a))+[2*k]*a;coefficient*=comb(m,a)
        value+=coefficient*sum(character(shape,tuple(sorted(vc))) for shape in shapes)
    value,remainder=divmod(value,2**len(cycles));assert not remainder
    return value

def flavour(shape,distinguished,colours):
    value=sum(Fraction(character(shape,c)*relation_trace(tuple(sorted(c+(1,)*distinguished)),colours),centralizer(c)) for c in partitions(sum(shape)))
    assert value.denominator==1 and value>=0
    return int(value)

def operators(shape):
    assert len(shape)==3 and all(n>=3 and n%2 for n in shape) and list(shape)==sorted(shape,reverse=True)
    k=shape[2];b=(shape[1]-k)//2;a=(shape[0]-shape[1])//2
    terminal=p.mul(p.power(p.D,k),p.mul(p.power(p.AA,a),p.power(p.ZZ,b)))
    sources=[p.mul(p.power(p.D,k-1),p.mul(p.PAIR,p.mul(p.power(p.AA,a),p.power(p.ZZ,b))))]
    recipes=['D^'+str(k-1)+' T AA^'+str(a)+' ZZ^'+str(b)]
    if a:
        contraction=p.add(*(p.mul(p.dot(p.X[0],p.X[i]),p.dot(p.X[0],p.Y[i])) for i in range(3)))
        sources.append(p.mul(p.power(p.D,k-1),p.mul(contraction,p.mul(p.power(p.AA,a-1),p.power(p.ZZ,b)))))
        recipes.append('D^'+str(k-1)+' sum_i dot(A,X_i)dot(A,Y_i) AA^'+str(a-1)+' ZZ^'+str(b))
    if b:
        contraction=p.add(*(p.mul(p.dot(p.Z,p.X[i]),p.dot(p.Z,p.Y[i])) for i in range(3)))
        sources.append(p.mul(p.power(p.D,k-1),p.mul(contraction,p.mul(p.power(p.AA,a),p.power(p.ZZ,b-1)))))
        recipes.append('D^'+str(k-1)+' sum_i dot(Z,X_i)dot(Z,Y_i) AA^'+str(a)+' ZZ^'+str(b-1))
    return sources,terminal,recipes,'D^'+str(k)+' AA^'+str(a)+' ZZ^'+str(b)

def run(colours,shape):
    start=time.monotonic();degree=sum(shape);assert degree==colours+2 and colours%2
    raw_sources,raw_terminal,recipes,terminal_recipe=operators(shape)
    for poly in raw_sources+[raw_terminal]:
        assert {h.physical_weight(m) for m in poly}=={shape}
        assert not h.raising(poly,0,1) and not h.raising(poly,1,2)
    source_dims={','.join(map(str,s)):flavour(s,1,colours) for s in predecessor_shapes(shape)}
    source_dimension=sum(source_dims.values());terminal_dimension=flavour(shape,0,colours)
    assert source_dimension==len(raw_sources) and terminal_dimension==source_dimension+1
    sources=list(map(p.Coefficient,raw_sources));terminal=p.Coefficient(raw_terminal)
    probes=[]
    for triples in [1,3,5]:
        if triples>min(shape):continue
        w=[(0,1,2)]*triples
        w += [(i,i) for i,n in enumerate(shape) for _ in range((n-triples)//2)]
        s,w=p.canonical(w);assert s;probes.append(w)
    layouts=[q for q in partitions(degree,2) if max(q)<=7 and len(q)>=5]
    rng=random.Random(20261123+101*colours+sum((i+1)*n for i,n in enumerate(shape)))
    samples=[i for i,n in enumerate(shape) for _ in range(n)]
    for _ in range(300):
        rng.shuffle(samples);layout=rng.choice(layouts);w=[];j=0
        for n in layout:w.append(tuple(samples[j:j+n]));j+=n
        s,w=p.canonical(w)
        if s:probes.append(w)
    matrix=[];words=[];examined=0;seen=set()
    for w in probes:
        if w in seen:continue
        seen.add(w);examined+=1
        values=[p.image(q,w) for q in sources]+[terminal(w)]
        if p.matrix_rank(matrix+[values])>len(matrix):
            matrix.append(values);words.append(w)
            print('PROBE',colours,shape,p.encode(w),list(map(str,values)),flush=True)
        if len(matrix)==terminal_dimension:break
    assert len(matrix)==terminal_dimension,(colours,shape,len(matrix))
    rank=p.matrix_rank([r[:-1] for r in matrix]);assert rank==source_dimension
    dual=p.solve(list(zip(*matrix)),[0]*source_dimension+[1])
    assert [sum(dual[i]*matrix[i][j] for i in range(len(words))) for j in range(len(words))]==[0]*source_dimension+[1]
    primitive_checks=[]
    for w in words:
        image=Counter()
        for u,c in h.primitive(w).items():
            for v,d in p.differential(u).items():image[v]+=c*d
        assert {v:c for v,c in image.items() if c}=={w:Fraction(1)}
        for i,j in [(0,1),(1,2)]:
            left=h.flavour_action(h.primitive(w),i,j);right=Counter()
            for v,c in h.flavour_action({w:1},i,j).items():
                for u,d in h.primitive(v).items():right[u]+=c*d
            assert left=={u:c for u,c in right.items() if c}
        value=h.q_primitive(w,terminal);assert value==terminal(w)
        primitive_checks.append({'word':p.encode(w),'Q0_cochain_coefficient':str(value)})
    return dict(passed=True,group='SO('+str(colours)+')',colours=colours,shape=shape,j=2*degree,
                cochain_length=degree-1,Taylor_diagram=transpose(shape),source_scalar_count=degree-2,matter_count=1,
                source_dimension=source_dimension,terminal_dimension=terminal_dimension,source_Pieri_multiplicities=source_dims,
                exact_Q0_rank=rank,exact_fortuitous_quotient_dimension=terminal_dimension-rank,
                source_marker_recipes=recipes,terminal_marker_recipe=terminal_recipe,
                source_marker_term_counts=[len(q) for q in raw_sources],terminal_marker_term_count=len(raw_terminal),
                ordered_probes=list(map(p.encode,words)),probe_matrix=[[str(c) for c in r] for r in matrix],
                exact_determinant=str(p.determinant(matrix)),dual_coefficients=list(map(str,dual)),
                source_and_terminal_highest_weight_checks=True,equivariant_primitive_checks=primitive_checks,
                representative='h.PrimitiveCoefficient(p.Coefficient(operators(shape)[1]))',
                scope='Exact classical Q0 fortuitous quotient. Historical novelty and quantum fate are separate questions.',
                probes_examined=examined,seconds=time.monotonic()-start)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--colours',type=int,required=True)
    ap.add_argument('--shape',type=int,nargs=3,action='append',required=True)
    ap.add_argument('--out',type=Path,default=Path(__file__).with_name('MINIMAL_SCALAR_CERTIFICATES.json'))
    ap.add_argument('--query-word');ap.add_argument('--native-sorted-theta',action='store_true')
    a=ap.parse_args();rows=[]
    if a.query_word:
        assert len(a.shape)==1 and sum(a.shape[0])==a.colours+2
        terminal=p.Coefficient(operators(tuple(a.shape[0]))[1])
        value=h.PrimitiveCoefficient(terminal)(p.decode(a.query_word))
        if a.native_sorted_theta:value*=(-1)**a.query_word.count('b')
        print(json.dumps({'colours':a.colours,'shape':a.shape[0],'word':a.query_word,'native_sorted_theta':a.native_sorted_theta,'cochain_coefficient':str(value)}));return
    for shape in a.shape:
        row=run(a.colours,tuple(shape));rows.append(row)
        a.out.write_text(json.dumps({'passed':True,'field':'Q','records':rows},indent=2)+'\n')
        print('CERTIFIED',a.colours,shape,row['exact_determinant'],row['dual_coefficients'],flush=True)

if __name__=='__main__':main()
