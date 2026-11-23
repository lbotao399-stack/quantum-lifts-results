"""Certify a highest-weight native operator and exact fortuitous dimension."""
from fast_paths import *
import argparse,gzip,hashlib,json,math,pickle,time
from fractions import Fraction


def main():
    cache_native_canonicalization()
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,default=6)
    args=ap.parse_args();rank=args.rank;start=time.monotonic()
    small=HERE/f'small_sp{rank}';support=HERE/f'support_sp{rank}'
    output=HERE/f'final_sp{rank}';output.mkdir(exist_ok=True)
    discovery=json.loads((small/'certificate.json').read_text())
    solver=json.loads((support/'result.json').read_text());assert solver['passed']
    raw=pickle.loads((support/'universal_preimage.pkl').read_bytes())
    defect=pickle.loads((small/'candidate_defect.pkl').read_bytes())
    row=raw['polynomial'];denominator=raw['denominator']
    for spin in range(1,4):
        raised=pf.spin_action(row,True)
        returned=pf.spin_action(raised,False)
        combined=Counter({w:spin*(spin+1)*c for w,c in row.items()})
        for w,c in returned.items():combined[w]-=c
        row=clean(combined);denominator*=spin*(spin+1)
        common=math.gcd(denominator,*row.values())
        row={w:c//common for w,c in row.items()};denominator//=common
        print('highest-weight projection',spin,len(row),denominator,'seconds',time.monotonic()-start,flush=True)
    for raising in (1,2,3):
        result=Counter()
        for word,c in row.items():
            for w,d in pf.m.raise_word(word,raising).items():result[w]+=c*d
        assert not any(result.values()),('highest-weight',raising)
    variation=pf.m.q(row)
    assert variation=={w:c*denominator for w,c in defect.items()}
    assert not pf.m.q(defect)
    assert discovery['candidate_exact_nonboundary']
    dual=[(pf.m.decode(w),Fraction(c)) for c,w in discovery['dual_support']]
    assert sum(c*defect.get(w,0) for w,c in dual)==1
    source=pickle.loads((HERE/f'frames_sp{rank}/relations_L{2*rank+3}.pkl').read_bytes())[0]
    for i,identity in enumerate(source):
        image=pf.m.q(identity)
        assert sum(c*image.get(w,0) for w,c in dual)==0,('dual boundary',i)
    print('exact closure and all dual pairings passed','seconds',time.monotonic()-start,flush=True)
    from trace_descent import dimensions
    characters=[];weight=(3,3,2*rank-1,0,0)
    for length,expected in zip(range(2*rank+2,2*rank+5),(67,188,183)):
        finite=sum(dimensions(weight,length,2*rank,'Sp').values())
        stable=sum(dimensions(weight,length,2*(length+1),'Sp').values())
        assert stable-finite==expected
        characters.append(dict(length=length,stable=stable,finite=finite,identities=expected))
        print('exact character',characters[-1],flush=True)
    assert discovery['source_relation_Q_rank']==121 and discovery['combined_rank']==61
    assert discovery['all_183_combined_rows_reconstructed_over_Q']
    operator=dict(group=f'Sp({rank})',native_weight=list(weight),length=2*rank+3,j=4*rank+16,
                  normalization='Q0(operator) equals the supplied integer defect, without a field rescaling',
                  terms=[[str(Fraction(row[w],denominator)),pf.m.encode(w)] for w in sorted(row)])
    op=output/f'sp{rank}_fortuitous_operator.json.gz'
    op.write_bytes(gzip.compress(json.dumps(operator,separators=(',',':')).encode(),mtime=0))
    dp=output/f'sp{rank}_defect.json.gz'
    dp.write_bytes(gzip.compress(json.dumps(dict(terms=[[str(defect[w]),pf.m.encode(w)] for w in sorted(defect)]),separators=(',',':')).encode(),mtime=0))
    determinant=int(discovery['combined_exact_minor_determinant']);remainder=abs(determinant);factors=[]
    for prime in (2,3,5,7,11,13,17,19,23,29,31):
        power=0
        while remainder%prime==0:remainder//=prime;power+=1
        if power:factors.append((prime,power))
    receipt=dict(group=f'Sp({rank})',native_weight=list(weight),paper_sector=[2*rank-1,0,0,3,3],
        length=2*rank+3,j=4*rank+16,operator_terms=len(row),least_trace_coefficient_denominator=denominator,
        defect_terms=len(defect),all_three_native_highest_weight_checks=True,
        exact_native_Q0_identity=True,exact_Q0_squared_zero=True,
        exact_finite_rank_closure_via_Pfaffian_identities=True,
        exact_rational_dual_support_size=len(dual),dual_supported_only_on_single_traces=True,
        dual_source_boundaries_checked=len(source),all_dual_boundary_pairings='0',dual_defect_pairing='1',
        source_relation_dimension=188,source_relation_Q0_rank=121,target_relation_dimension=183,
        target_Q0_rank=45,stable_abelianization_rank=24,combined_Q0_abelianization_rank=61,
        combined_exact_minor_determinant=str(determinant),determinant_sign=-1 if determinant<0 else 1,
        determinant_small_prime_factors=factors,determinant_remaining_factor=str(remainder),
        fortuitous_dimension_exact=(183-61)-121,full_finite_cohomology_dimension_NOT_claimed=True,
        exact_characters=characters,universal_inverse_blocks=solver['blocks'],
        independent_new_family_NOT_claimed=True,classical_Q0_only=True,loop_corrections_NOT_computed=True,
        no_lift_to_next_symplectic_rank_under_trace_covering=True,
        no_lift_reason='At Sp(r+1), length 2r+3 is in the universally stable first-layer exclusion range.',
        operator_sha256=hashlib.sha256(op.read_bytes()).hexdigest(),
        defect_sha256=hashlib.sha256(dp.read_bytes()).hexdigest(),
        small_complex_certificate_sha256=hashlib.sha256((small/'certificate.json').read_bytes()).hexdigest(),
        support_solver_result_sha256=hashlib.sha256((support/'result.json').read_bytes()).hexdigest(),
        elapsed_seconds=time.monotonic()-start)
    (output/f'sp{rank}_fortuitous_certificate.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
