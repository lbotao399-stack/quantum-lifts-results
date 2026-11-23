"""Reproduce the support-based inverse-Q validation on the published Sp(5) class."""
from fast_paths import *
from fractions import Fraction
import gzip,json,math,pickle,subprocess,time

cache_native_canonicalization();assert BASE_PACKAGE is not None
start=time.monotonic()
data=json.loads(gzip.decompress((BASE_PACKAGE/'symplectic/sp5/sp5_new_fortuitous_operator.json.gz').read_bytes()))
polynomial={pf.m.decode(w):Fraction(c) for c,w in data['terms']}
denominator=math.lcm(*(c.denominator for c in polynomial.values()))
defect=pf.m.q({w:int(c*denominator) for w,c in polynomial.items()})
assert len(defect)==42480
path=HERE/'known_defect_sp5.pkl';path.write_bytes(pickle.dumps(defect,protocol=4))
output=HERE/'support_sp5'
while not (output/'result.json').exists():
    subprocess.run([sys.executable,str(CODE_DIR/'support_solver.py'),'--defect',str(path),
                    '--output',str(output),'--max-sources','12000','--seconds','600'],check=True)
result=json.loads((output/'result.json').read_text());assert result['passed']
print(json.dumps(dict(passed=True,blocks=result['blocks'],exact_native_Q0_identity=True,
                      existing_class_only=True,seconds=time.monotonic()-start)),flush=True)
