"""Compare left derivatives of the odd superfield, archive words and (4.84) over Q.
Run at the repository root: python3 algorithms/verify_g2_superfield_sign.py .
Taylor coordinates are always immediately to the left of the coefficient matrix.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations, product
import json, gzip
from pathlib import Path
from verify_g2_current_product import HERE, FIRST, SECOND, canonical, fw, polynomial, factors
from verify_g2_s3_currents import add, scale, exterior, diagonal, probe
from fortuity.full.letters import letters_up_to, is_ferm, Q_letter

AXES=((1,0),(0,1))
ZERO=diagonal({}, {})

def accum(items):
    d=defaultdict(F)
    for k,v in items:d[k]+=v
    return {k:v for k,v in d.items() if v}

def epsilon(i,j,k):
    return (-1)**((i>j)+(i>k)+(j>k))

def odd_derivative(terms, axis):
    # Each key contains the ordered theta monomial followed by the operator word.
    return accum((((s[:i]+s[i+1:],w),c*(-1)**i) for (s,w),c in terms.items()
                  for i in range(len(s)) if s[i]==axis))

def super_product(u,v):
    terms=[]
    for (s,w),a in u.items():
        for (t,y),b in v.items():
            if set(s).intersection(t):continue
            swaps=sum(is_ferm(l) for l in w)*len(t)+sum(i>j for i in s for j in t)
            terms.append(((tuple(sorted(s+t)),w+y),a*b*(-1)**swaps))
    return accum(terms)

def at_origin(u):return canonical({(w,):c for (s,w),c in u.items() if not s})

def graviton(flavour, direction):
    psi={(l[2],(l,)):1 for l in letters_up_to(0,0)}
    # At z=0 the even derivative removes exactly one power of that coordinate.
    dz={(l[2],(l,)):1 for l in letters_up_to(*AXES[direction]) if l[:2]==AXES[direction]}
    def bilinear(s):return super_product(odd_derivative(psi,s),dz)
    s,t=(flavour+1)%3,(flavour+2)%3
    a,b=at_origin(odd_derivative(bilinear(s),t)),at_origin(odd_derivative(bilinear(t),s))
    return at_origin(bilinear(flavour)),add(scale(a,F(1,2)),scale(b,F(-1,2)))

def mult(*factors_):
    out={():F(1)}
    for factor in factors_:
        out=canonical(accum((w+v,c*b) for w,c in out.items() for v,b in factor.items()))
    return out

def even_derivative(expr,direction):
    terms=[]
    for w,c in expr.items():
        for k,tr in enumerate(w):
            for i,l in enumerate(tr):
                nxt=(l[0]+AXES[direction][0],l[1]+AXES[direction][1],l[2]);rep=tr[:i]+(nxt,)+tr[i+1:]
                terms.append((w[:k]+(rep,)+w[k+1:],c*(l[direction]+1)))
    return canonical(accum(terms))

def apply(expr,fields):
    vals=[]
    for w,c in expr.items():
        term={0:c}
        for tr in w:
            components=[]
            for color in range(7):
                v={0:1}
                for l in tr:v=exterior(v,fields.get(l,ZERO)[color])
                components.append(v)
            term=exterior(term,add(*components))
        vals.append(term)
    return add(*vals)

def qzero(fields,n):
    count=0
    for l in letters_up_to(n,n):
        for color in range(7):
            val=add(*(scale(exterior(fields.get(x,ZERO)[color],fields.get(y,ZERO)[color]),c) for c,x,y in Q_letter(l)))
            assert not val;count+=1
    return count

def run():
    data=json.loads(gzip.decompress((HERE.parent/'data/N4_G2_np_partner_replay_20261118.json.gz').read_bytes()))
    archived=polynomial(data['representative']);fields={(0,0,(0,)):diagonal({0:1},{}),(1,0,()):diagonal({1:1},{}),
                    (0,1,()):diagonal({2:1},{}),(0,0,(1,2)):diagonal({4:1},{})}
    checks=qzero(fields,1); assert sum(x*x for x in FIRST)==4
    # Evaluate the printed minus sign in matter-gaugino order; both flavour orientations +1.
    printed=add(*(scale(exterior({1<<i:4},exterior({4:1},{1<<(1-i):4})),F(-(-1)**i,6)) for i in (0,1)))
    components={(s,i):graviton(s,i) for s,i in product(range(3),range(2))}
    for k,(l,g) in components.items():
        assert (l,g)==factors(*k) and not canonical(fw.Q_expr(l)) and not canonical(fw.Q_expr(g))
    lg=add(*(scale(mult(components[s,i][0],components[s,1-i][1]),F((-1)**i,6)) for s,i in components))
    assert lg==archived and apply(lg,fields)=={7:F(-16,3)} and printed=={7:F(16,3)}
    bilinear=canonical({(((1,0,()),(0,1,())),):1})
    six1=mult(bilinear,lg);six2=add(*(scale(mult(components[s,i][0],components[t,j][0],
                even_derivative(components[h,1-j][0],1-i)),epsilon(s,t,h)*(-1)**(i+j))
                for s,t,h in permutations(range(3)) for i,j in product(range(2),repeat=2)))
    for expr in (lg,six1,six2):assert not canonical(fw.Q_expr(expr))
    arch_rows=[]; print_rows=[]; checks6=0
    for matter,mask in ((True,31),(False,13)):
        jet={(1,0,()):diagonal({1:1},{2:1}),(0,1,()):diagonal({4:1},{8:1}),(0,0,(0,)):diagonal({0:1},{})}
        if matter:jet[0,0,(1,2)]=diagonal({16:1},{})
        else:jet.update({(0,0,(1,)):diagonal({},{0:1}),(1,0,(2,)):diagonal({0:1},{})})
        checks6+=qzero(jet,2)
        arch_rows.append([apply(x,jet).get(mask,0) for x in (six1,six2)])
        a,b,_=probe(matter);print_rows.append([a.get(mask,0),b.get(mask,0)])
    assert arch_rows==[[64,0],[16,576]] and print_rows==[[-64,0],[16,576]]
    det=lambda rows:int(rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]); assert det(arch_rows)==-det(print_rows)==36864
    return dict(archive_terms=len(lg),exterior_generator_order=['gaugino_plus','gaugino_minus','first_matter_fermion'],
                three_fermion_Q0_checks=checks,three_fermion_archive_coefficient=str(apply(archived,fields)[7]),
                three_fermion_printed_both_flavour_orientations_positive=str(printed[7]),
                superfield_bilinear_pairs_equal_to_archived=6,closed_six_letter_expressions=[len(six1),len(six2)],
                six_letter_Q0_checks=checks6,superfield_six_letter_minor=[[int(x) for x in r] for r in arch_rows],
                superfield_determinant=det(arch_rows),printed_minor=[[int(x) for x in r] for r in print_rows],
                printed_determinant=det(print_rows),instanton_matrix_element_computed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);args=p.parse_args();assert __debug__
    text=json.dumps(run(),indent=2)+'\n';print(text,end='')
    if args.output:args.output.write_text(text)
