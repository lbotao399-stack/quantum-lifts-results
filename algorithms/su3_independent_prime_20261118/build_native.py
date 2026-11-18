"""Build the supplied SU(3) production libraries with a C++17 compiler."""
from pathlib import Path
import argparse,subprocess

ROOT=Path(__file__).resolve().parent/'runtime/work/sp3-native.nosync'
parser=argparse.ArgumentParser();parser.add_argument('--compiler',default='c++');args=parser.parse_args()
targets=[('su_full_projection_native_v1.cpp','su_full_projection_native_v1.dylib',[]),
         ('su_quantum_integer_v1.cpp','su_quantum_integer_v1.dylib',[])]
for prime in (1000003,1000033):
    for suffix,source in (('v1','bmn_derivation_rank_v1.cpp'),('threadsafe_v1','bmn_derivation_rank_threadsafe_v1.cpp')):
        targets.append((source,'su_derivation_d3_p%d_%s.dylib'%(prime,suffix),
            ['-DTRACE_FRAME_MATRIX_SIZE=3','-DTRACE_FRAME_PRIME='+str(prime)]))
for source,output,flags in targets:
    command=[args.compiler,'-O3','-std=c++17','-shared','-fPIC',*flags,source,'-o',output]
    subprocess.run(command,cwd=ROOT,check=True);print(output,flush=True)
