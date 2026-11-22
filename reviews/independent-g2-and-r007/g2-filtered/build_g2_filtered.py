"""Use the published jet-pair kernel with the two exceptional colour contractions."""
from argparse import ArgumentParser
from pathlib import Path
import subprocess


def build(repo, native, root):
    header = (repo/'algorithms/dual_singlets_20261118/taylor_hooks/native/q1.h').read_text()
    begin = header.index('static void psi3(')
    end = header.index('static u64 q1_eval_fast(', begin)
    header = header[:begin] + (root/'g2_filter_colour.h').read_text()+'\n'+header[end:]
    header = header.replace('#include "d1.h"', '')
    (native/'g2_filter_kernel.h').write_text(header)
    subprocess.run(['c++','-O3','-std=c++17','-pthread','-I'+str(native),
                    str(root/'g2_filter_q1.cpp'),'-o',str(native/'g2_filter_q1')],check=True)


if __name__ == '__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--native',type=Path,required=True)
    a=p.parse_args()
    build(a.repo,a.native,Path(__file__).resolve().parent)
