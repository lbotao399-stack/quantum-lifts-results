"""Exact non-exactness detector for C0 star C4 and absence of incoming higher pages."""
from pathlib import Path
from fractions import Fraction as F
import pickle,json,time,subprocess
from flint import fmpq,fmpq_mat
import certify_c4 as C
import exact_product as E
ROOT=Path(__file__).resolve().parent
C.CHECK_PATH=ROOT/'nonzero_checks.json'
START=time.time()


def indices(polys):
    rows=[];selected=[]
    for i,poly in enumerate(polys):
        row={w:C.modq(c) for w,c in poly.items() if c}
        for pivot,old in rows:
            value=row.get(pivot)
            if value:
                for w,c in old.items():
                    new=(row.get(w,0)-value*c)%C.P
                    if new:row[w]=new
                    else:row.pop(w,None)
        if row:
            pivot=min(row);scale=pow(row[pivot],C.P-2,C.P)
            rows.append((pivot,{w:c*scale%C.P for w,c in row.items()}))
            selected.append(i)
    return selected


def dimension(key):
    Nb,Nc,n=key;values=[]
    for n1,n2,weight in [(n,n,(0,0)),(n,n,(1,-1)),(n+1,n-1,(0,0)),(n+1,n-1,(1,-1))]:
        run=subprocess.run(['/Users/libotao/np_lifts_engine/bf','2',str(n1),str(n2),'DIM',str(Nb),str(Nc),*map(str,weight)],capture_output=True,text=True,check=True)
        values.append(json.loads(run.stdout)['dim'])
    return values,values[0]-values[1]-values[2]+values[3]


def main():
    full=C.load('full_spaces_mod.pkl')
    incoming=C.load('incoming_spaces_mod.pkl')
    allspaces={**full,**incoming}
    verified={};dimensions={}
    for key,space in allspaces.items():
        basis=[C.poly(p) for p in space['basis']]
        closed=[C.poly(p) for p in space['closed']]
        counts,dim=dimension(key)
        assert len(basis)==dim
        assert len(indices(basis))==dim
        assert len(indices(closed))==len(closed)
        images=[E.Q0(p) for p in basis]
        assert len(indices(images))==space['rank']
        for i,p in enumerate(closed):C.check('%s exact kernel vector %d'%(key,i),E.Q0(p))
        assert len(closed)+space['rank']==dim
        verified[key]=dict(basis=basis,closed=closed,images=images,rank=space['rank'])
        dimensions[str(key)]=dict(weight_counts=counts,dimension=dim,Q0_rank=space['rank'],kernel=len(closed))
        C.say('FULL EXACT SPACE %s: dim %d, Q0 rank %d, kernel %d'%(key,dim,space['rank'],len(closed)))
    assert verified[(4,2,6)]['rank']+verified[(4,3,6)]['rank']==len(verified[(4,3,6)]['basis'])
    for p in verified[(4,2,6)]['images']:C.check('incoming primitive differential squares to zero',E.Q0(p))
    for key in [(5,3,6),(6,3,6)]:
        counts,dim=dimension(key);assert dim==0
        dimensions[str(key)]=dict(weight_counts=counts,dimension=dim)
    C.say('E1 incoming Nb>=4,Nc=3 at (6,6) VANISHES over Q')

    closed=verified[(3,3,6)]['closed']
    images=[E.cartan(E.Q1(p)) for p in closed]
    select=indices(images);independent=[images[i] for i in select]
    assert len(independent)==4
    certificate=C.load('c4_hidden_product_exact.pkl')
    target=E.cartan(certificate['r2'])
    coordinates=C.independent_coordinate_rows(independent+[target])
    fq=lambda x:fmpq(F(x).numerator,F(x).denominator)
    matrix=fmpq_mat([[fq(p.get(w,0)) for w in coordinates] for p in independent+[target]])
    value=matrix.solve(fmpq_mat([[0],[0],[0],[0],[1]]))
    functional={w:F(int(value[i,0].numerator),int(value[i,0].denominator)) for i,w in enumerate(coordinates) if value[i,0]}
    evaluate=lambda p:sum(c*p.get(w,0) for w,c in functional.items())
    assert all(evaluate(p)==0 for p in images)
    assert evaluate(target)==1
    C.say('Exact Cartan functional: annihilates all 28 Q1(kernel Q0) images, evaluates r2 to1')
    for w,c in functional.items():C.say(str(c)+' '+repr(w))
    result=dict(spaces=verified,dimensions=dimensions,cartan_images=images,functional=functional,
        target=target,value=F(1),coordinates=coordinates,
        statement='C0 star C4 is nonzero on E2; incoming higher differentials have zero E1 sources, so the full product cocycle is nonzero in H(Qpert).')
    with (ROOT/'C0C4_nonzero_exact.pkl').open('wb') as f:pickle.dump(result,f)
    (ROOT/'nonzero_summary.json').write_text(json.dumps(dict(passed=True,seconds=time.time()-START,
        dimensions=dimensions,cartan_image_rank=4,functional_terms=len(functional),value=1,
        scope=result['statement']),indent=2))


if __name__=='__main__':main()
