#!/usr/bin/env python3
"""Direct flavour-character count, compared with the Taylor-character count."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json,sys,time

HERE=Path(__file__).resolve().parent
DEFAULT=HERE.parents[1]/'quantum_lifts_frontier_work/root_audit/so11_pin/algorithms/dual_singlets_20261118'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--character-source',type=Path,default=DEFAULT);args=ap.parse_args()
    sys.path.insert(0,str(args.character_source));t=time.monotonic()
    from b2_schur_complex import centralizer,character,partitions
    from orthogonal_glueball_sources import missing_adjoint_trace,relations
    def direct(shape,distinguished):
        value=sum(Fraction(character(shape,c)*missing_adjoint_trace(tuple(sorted((1,)*distinguished+c)),11),centralizer(c)) for c in partitions(sum(shape)))
        assert value.denominator==1;return int(value)
    terminal=direct((5,5,5),0);source=direct((5,4,4),1)
    assert (source,terminal)==(21,17)
    taylor={str(y):relations((3,3,3,3,3),y,11) for y in range(12,16)}
    assert taylor['14']==source and taylor['15']==terminal
    stable={str(y):relations((3,3,3,3,3),y,0) for y in [14,15]}
    r={'passed':True,'colours':11,'scalar_weight':[5,5,5],
       'terminal_scalar_Schur_shape':[5,5,5],'predecessor_scalar_Schur_shape':[5,4,4],
       'predecessor_distinguished_adjoint_letters':1,'direct_flavour_relation_dimensions':{'source':source,'terminal':terminal},
       'Taylor_relation_dimensions':taylor,'stable_cochain_dimensions':stable,
       'finite_cochain_dimensions':{'14':stable['14']-source,'15':stable['15']-terminal},
       'input_sha256':{name:hashlib.sha256((args.character_source/name).read_bytes()).hexdigest() for name in ['b2_schur_complex.py','orthogonal_glueball_sources.py']},
       'seconds':time.monotonic()-t,
       'scope':'Direct GL3flavour branching and full Taylor-character routes agree; both use the same pinned orthogonal invariant-tensor character. This is not an independent derivation of the orthogonal second fundamental theorem.'}
    (HERE/'character_completeness.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
