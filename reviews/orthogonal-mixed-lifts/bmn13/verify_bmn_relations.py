"""Rank certificates for the two dressed SO(13) relations in flavour (7,5,5)."""
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import argparse,json,subprocess,tempfile
from pathlib import Path
import sys
from flint import nmod_mat

def pivots(matrix):
    reduced,rank=matrix.rref()
    return [next(j for j in range(matrix.ncols()) if reduced[i,j]) for i in range(rank)]

def witness(rows,p):
    matrix=nmod_mat(rows,p)
    ri=pivots(matrix.transpose());ci=pivots(matrix)
    det=int(nmod_mat([[rows[i][j] for j in ci] for i in ri],p).det())
    assert det
    return dict(rank=len(ri),rows=ri,columns=ci,determinant=det)

def run(a):
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from orthogonal_glueball_sources import relations
    from b2_schur_complex import transpose
    shape=transpose((7,5,5));dims=[relations(shape,y,13) for y in range(13,18)]
    assert dims==[0,9,35,69,42]
    counts={1:100,2:8,3:8,4:2,6:40,7:90,10:20,11:5}
    specs=sorted(counts)
    proof=dict(colours=13,flavour=[7,5,5],lengths=list(range(13,18)),relation_dimensions=dims,
               recipe_counts=counts,primes={})
    with tempfile.TemporaryDirectory() as tmp:
        binary=Path(tmp)/'probe';subprocess.run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('bmn_relation_probe.cpp')),'-o',str(binary)],check=True)
        for p in (1000003,1000033):
            def collect(mode):
                result=subprocess.run([str(binary),'13','5','2','100','0',str(counts[mode]),str(mode),str(p)],
                                      check=True,capture_output=True,text=True)
                fields=[list(map(int,line.split())) for line in result.stdout.splitlines()]
                assert [r[0] for r in fields]==list(range(counts[mode])) and all(len(r)==101 for r in fields)
                return [r[1:] for r in fields]
            with ThreadPoolExecutor(a.workers) as pool:bank=dict(zip(specs,pool.map(collect,specs)))
            kernel=bank[6]+bank[7]+bank[10]+bank[11]
            image=bank[1]+bank[2]+bank[3]
            cert={name:witness(rows,p) for name,rows in [('kernel',kernel),('image',image),('with_products',image+bank[4])]}
            ranks=[d['rank'] for d in cert.values()]
            assert ranks==[29,40,42] and ranks[0]+ranks[1]==dims[-2] and ranks[-1]==dims[-1]
            assert cert['with_products']['rows'][-2:]==[len(image),len(image)+1]
            proof['primes'][p]=cert
            print(p,[(name,r['rank'],r['determinant']) for name,r in cert.items()],flush=True)
    a.out.write_text(json.dumps(proof,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();assert 1<=args.workers<=4;run(args)
