#!/usr/bin/env python3
"""Verify the frozen Sp(6) cochain and finite-dimensional proof over Q."""
from fast_paths import *
from fractions import Fraction
from flint import fmpz_mat
import gzip,hashlib,json,math,pickle,time
import abelian_obstruction as ab
from assemble_operator import operator_bytes


def check_combination(row,coefficients,basis):
    coefficients=list(map(Fraction,coefficients))
    denominator=math.lcm(*(c.denominator for c in coefficients))
    result=Counter()
    for c,source in zip(coefficients,basis):
        if c:
            multiplier=int(c*denominator)
            for w,d in source.items():result[w]+=multiplier*d
    assert clean(result)=={w:c*denominator for w,c in row.items()}


def main():
    cache_native_canonicalization();start=time.monotonic()
    foreign=[]
    for module in sys.modules.values():
        path=getattr(module,'__file__',None)
        if path and ('quantum_lifts_frontier_source' in path or 'so11_pin' in path):foreign.append(path)
    assert not foreign,foreign
    frozen=CODE_DIR/'certificates';data=CODE_DIR/'data'
    certificate=json.loads((frozen/'sp6_fortuitous_certificate.json').read_text())
    small_bytes=gzip.decompress((frozen/'sp6_small_complex.json.gz').read_bytes())
    assert hashlib.sha256(small_bytes).hexdigest()==certificate['small_complex_certificate_sha256']
    small=json.loads(small_bytes)
    coordinates=json.loads(gzip.decompress((frozen/'sp6_rational_coordinates.json.gz').read_bytes()))
    packed_operator=operator_bytes(data);dp=data/'sp6_defect.json.gz'
    assert hashlib.sha256(packed_operator).hexdigest()==certificate['operator_sha256']
    assert hashlib.sha256(dp.read_bytes()).hexdigest()==certificate['defect_sha256']
    operator=json.loads(gzip.decompress(packed_operator))
    rational={pf.m.decode(w):Fraction(c) for c,w in operator['terms']}
    denominator=math.lcm(*(c.denominator for c in rational.values()))
    assert denominator==certificate['least_trace_coefficient_denominator']
    integral={w:int(c*denominator) for w,c in rational.items()}
    defect={pf.m.decode(w):int(c) for c,w in json.loads(gzip.decompress(dp.read_bytes()))['terms']}
    assert len(integral)==certificate['operator_terms'] and len(defect)==certificate['defect_terms']
    for word in integral:
        letters=[l for tr in word for l in tr]
        assert len(letters)==15
        assert (sum(l[0] for l in letters),sum(l[1] for l in letters),
                *(sum(i in l[2] for l in letters) for i in range(3)))==(3,3,11,0,0)
    for raising in (1,2,3):
        result=Counter()
        for w,c in integral.items():
            for v,d in pf.m.raise_word(w,raising).items():result[v]+=c*d
        assert not any(result.values())
    assert pf.m.q(integral)=={w:c*denominator for w,c in defect.items()}
    assert not pf.m.q(defect) and not ab.restrict(defect)
    frames={length:pickle.loads((HERE/f'frames_sp6/relations_L{length}.pkl').read_bytes())[0]
            for length in (14,15,16)}
    for length,expected in zip((14,15,16),(67,188,183)):
        assert len(frames[length])==expected
        e=pf.Echelon()
        for row in frames[length]:assert e.add(row)
    check_combination(defect,small['candidate_relation_coefficients'],frames[16])
    print('frozen cochain, highest weights and exact finite closure passed',flush=True)
    from fortuity.full.letters import letters_up_to,Q_letter
    for letter in letters_up_to(3,3):
        twice=Counter()
        for c,left,right in Q_letter(letter):
            for d,a,b in Q_letter(left):twice[a,b,right]+=c*d
            for d,a,b in Q_letter(right):twice[left,a,b]+=c*d*(-1)**pf.is_ferm(left)
        assert not any(twice.values())
    e=pf.Echelon()
    for i,row in enumerate(frames[14]):
        image=pf.m.q(row);assert e.add(image)
        check_combination(image,coordinates['preceding_image_source_coordinates'][i],frames[15])
    assert len(e.rows)==67
    dual=[(pf.m.decode(w),Fraction(c)) for c,w in small['dual_support']]
    assert all(len(w)==1 for w,c in dual)
    e=pf.Echelon()
    for row in frames[15]:
        image=pf.m.q(row);e.add(image)
        assert sum(c*image.get(w,0) for w,c in dual)==0
    assert len(e.rows)==121
    assert sum(c*defect.get(w,0) for w,c in dual)==1
    print('67 independent first images; all 188 dual boundary pairings passed',flush=True)
    columns=[(kind,pf.m.decode(w)) for kind,w in small['combined_minor_columns']]
    matrix=[]
    for i in small['combined_minor_rows']:
        q=pf.m.q(frames[16][i]);restriction=ab.restrict(frames[16][i])
        matrix.append([(q if kind=='Q0' else restriction).get(w,0) for kind,w in columns])
    determinant=int(fmpz_mat(matrix).det())
    assert len(matrix)==len(columns)==61
    assert str(determinant)==small['combined_exact_minor_determinant']!=str(0)
    from trace_descent import dimensions
    characters=[]
    for length,expected in zip((14,15,16),(67,188,183)):
        finite=sum(dimensions((3,3,11,0,0),length,12,'Sp').values())
        stable=sum(dimensions((3,3,11,0,0),length,2*(length+1),'Sp').values())
        assert stable-finite==expected
        characters.append(dict(length=length,stable=stable,finite=finite,identities=expected))
    assert characters==certificate['exact_characters']
    result=dict(passed=True,group='Sp(6)',native_weight=[3,3,11,0,0],length=15,j=40,
        operator_terms=len(integral),defect_terms=len(defect),
        exact_operator_Q0_defect_identity=True,exact_Pfaffian_closure=True,
        exact_dual_source_images_checked=188,dual_defect_pairing='1',
        source_relation_Q0_rank_over_Q=121,target_relation_dimension=183,
        combined_integer_minor_dimension=61,combined_integer_minor_determinant=str(determinant),
        exact_fortuitous_dimension=(183-61)-121,field='Q',classical_Q0_only=True,
        full_finite_sector_dimension_NOT_claimed=True,exact_characters=characters,
        native_source_relative_path=str(Path(pf.m.__file__).resolve().relative_to(BASE_PACKAGE)),
        elapsed_seconds=time.monotonic()-start)
    (HERE/'portable_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
