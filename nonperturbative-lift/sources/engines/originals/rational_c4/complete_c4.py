"""Complete the C4 cocycle and the hidden-product exact boundary over Q."""
from pathlib import Path
import pickle,json,time,resource,hashlib
from fractions import Fraction as F
import certify_c4 as C
import exact_product as E
import exact_star as S

ROOT=Path(__file__).resolve().parent
C.CHECK_PATH=ROOT/'completion_checks.json'
START=time.time()


def save(name,data):
    with (ROOT/name).open('wb') as f:pickle.dump(data,f)


def q(n,p):
    return E.Q0(p) if n==0 else E.T.qn(n,p)


def product_coefficient(left,right,r):
    out={}
    for i,A in enumerate(left):
        for j,B in enumerate(right):
            n=r-i-j
            if n<0:continue
            C.say('  product (%d,%d,%d), input terms (%d,%d)'%(i,j,n,len(A),len(B)))
            op=S.star_n(n,A,B)
            out=E.add(out,op)
    return out


def main():
    certificate=C.load('c4_hidden_product_exact.pkl')
    bases=C.load('completion_bases_mod.pkl')
    C0=certificate['C0']
    C2=certificate['C2']
    C4=certificate['C4_partial'][:]
    for index,key in [(3,(1,2,5)),(4,(0,2,5))]:
        target=E.sc(E.add(*(q(index-i,p) for i,p in enumerate(C4))),-1)
        source=[C.poly(p) for p in bases[key]]
        coeffs,_=C.rational_solve([E.Q0(p) for p in source],target)
        C4.append(E.add(*(E.sc(p,c) for p,c in zip(source,coeffs))))
        C.check('C4 full closure at loop %d'%index,E.add(E.Q0(C4[-1]),E.sc(target,-1)))
    save('C4_full_exact.pkl',dict(C0=C0,C2=C2,C4_source=C4))
    C.say('FULL C4 COCYCLE COMPLETE: terms '+repr(list(map(len,C4))))
    P=certificate['C2_squared'][:]
    L=certificate['C0_C4'][:]
    for r in (1,2):
        C.check('direct forest matches C2 product loop %d'%r,
                E.sub(product_coefficient(C2,C2,r),P[r]))
        C.check('direct forest matches C0C4 product loop %d'%r,
                E.sub(product_coefficient([C0],C4,r),L[r]))
    for r in (3,4):
        P.append(product_coefficient(C2,C2,r))
        L.append(product_coefficient([C0],C4,r))
        C.check('C2 product full closure loop %d'%r,
                E.add(*(q(r-i,p) for i,p in enumerate(P))))
        C.check('C0C4 product full closure loop %d'%r,
                E.add(*(q(r-i,p) for i,p in enumerate(L))))
        save('full_products_partial.pkl',dict(C0=C0,C2=C2,C4_source=C4,P=P,L=L))
    C.say('FULL PRODUCTS COMPLETE')
    prim=[E.sub(certificate['C2_squared_primitives'][0],E.sc(certificate['C0_C4_primitives'][0],3)),
          E.add(E.sub(certificate['C2_squared_primitives'][1],E.sc(certificate['C0_C4_primitives'][1],3)),certificate['boundary_Q1']),
          certificate['boundary_Q0']]
    for r in range(3):
        C.check('hidden product full boundary loop %d'%r,
                E.sub(E.sub(P[r],E.sc(L[r],3)),E.add(*(q(r-i,p) for i,p in enumerate(prim[:r+1])))))
    for r,key in [(3,(1,3,6)),(4,(0,3,6))]:
        target=E.sub(E.sub(P[r],E.sc(L[r],3)),E.add(*(q(r-i,p) for i,p in enumerate(prim))))
        source=[C.poly(p) for p in bases[key]]
        coeffs,_=C.rational_solve([E.Q0(p) for p in source],target)
        prim.append(E.add(*(E.sc(p,c) for p,c in zip(source,coeffs))))
        C.check('hidden product full boundary loop %d'%r,E.sub(target,E.Q0(prim[-1])))
    save('full_hidden_product_exact.pkl',dict(C0=C0,C2=C2,C4_source=C4,P=P,L=L,boundary=prim))
    summary=dict(passed=True,seconds=time.time()-START,checks=len(C.REPORT),
        maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        relation='C2 star C2 - 3 C0 star C4_source = Qpert boundary',
        article_normalization='C4_article = -C4_source, so C2 star C2 + 3 C0 star C4_article = Qpert boundary')
    (ROOT/'full_completion_summary.json').write_text(json.dumps(summary,indent=2))
    C.say(json.dumps(summary))


if __name__=='__main__':main()
