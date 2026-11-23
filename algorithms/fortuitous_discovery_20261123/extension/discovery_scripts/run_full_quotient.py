"""Repeat the established stable-image test on all 183 target identities.

All shared discovery checkpoints are read-only. The links below avoid
duplicating expensive inputs while isolating every new output.
"""
from pathlib import Path
import argparse
import importlib
import sys

HERE=Path(__file__).resolve().parent
ORIGINAL=HERE.parent/'fortuity_discovery_20261123'
sys.path.insert(0,str(ORIGINAL/'jets'))


def link(destination,source):
    if destination.is_symlink():
        assert destination.resolve()==source.resolve()
    elif not destination.exists():
        destination.symlink_to(source)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,choices=[4,5],required=True)
    args=ap.parse_args();rank=args.rank
    work=HERE/'jets';work.mkdir(exist_ok=True)
    link(HERE/'sp5_universal',ORIGINAL/'sp5_universal')
    stem='337' if rank==4 else '339'
    names=[f'exact_relations_sp{rank}_L{2*rank+2}.pkl',f'exact_relations_sp{rank}_L{2*rank+3}.pkl',
           f'exact_q_relations_sp{rank}_L{2*rank+2}.pkl']
    if rank==4:
        names += ['formal_337_L11.pkl','q_337_L11.pkl','exact_q_relations_sp4_L11.pkl']
    else:
        names += ['sp5_q13_identity_rows']
    for name in names:link(work/name,ORIGINAL/'jets'/name)
    target=f'complete_target_sp{rank}.pkl'
    link(work/target,HERE/target)
    module=importlib.import_module('relation_quotient' if rank==4 else 'sp5_relation_quotient')
    module.HERE=work
    module.main(2147483647,300,target)


if __name__=='__main__':main()
