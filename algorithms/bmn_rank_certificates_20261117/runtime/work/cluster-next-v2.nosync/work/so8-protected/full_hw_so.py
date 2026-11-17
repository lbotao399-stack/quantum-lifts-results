"""SO trace highest weights for full letters, including dotted SU(2).

Flavour rules and the E1-then-E2 reduction follow pinned fortuity/full/hw.py
and fullpar.py. Orthogonal transpose identities and the derivative raising
operator are imposed here. All kernel coefficients use the reported prime.
"""
from quantum_j18 import *
from hw_scan import formal_closed_kernel
from flint import nmod_mat
import pickle,os


def canonical_sequence(seq):
    return min(v[i:]+v[:i] for v in (seq,seq[::-1]) for i in range(len(v)))


def shape(w,mode='flavour'):
    def lab(L):
        a,b,S=L
        if mode=='e1':return a,b,int(0 in S)+int(1 in S),int(2 in S)
        if mode=='spin':return a+b,len(S)
        return a,b,len(S)
    return tuple(sorted(canonical_sequence(tuple(lab(L) for L in t)) for t in w))


def raise_word(w,i):
    out=defaultdict(int)
    for k,t in enumerate(w):
        for pos,(a,b,S) in enumerate(t):
            if i==3:
                if not b:continue
                factor=a+1;L=(a+1,b-1,S)
            else:
                src,dst=(1,0) if i==1 else (2,1)
                if src not in S or dst in S:continue
                factor=1;L=(a,b,tuple(sorted(dst if x==src else x for x in S)))
            sg,v=so_canon(w[:k]+(t[:pos]+(L,)+t[pos+1:],)+w[k+1:])
            if sg:out[v]+=sg*factor
    return {w:c for w,c in out.items() if c}


def raise_expr(e,i,p):
    out=defaultdict(int)
    for w,c in e.items():
        for v,d in raise_word(w,i).items():out[v]+=int(c)*d
    return {w:c%p for w,c in out.items() if c%p}


def combine(exprs,row,p):
    out=defaultdict(int)
    for e,c in zip(exprs,row):
        if c:
            for w,d in e.items():out[w]+=int(c)*int(d)
    return {w:c%p for w,c in out.items() if c%p}


def kernel(exprs,i,p):
    images=[raise_expr(e,i,p) for e in exprs]
    H,_=formal_closed_kernel(images,p)
    return [combine(exprs,row,p) for row in H]


def hw_flavour(words,p):
    classes=defaultdict(list)
    for w in words:classes[shape(w,'e1')].append({w:1})
    E1=[]
    for es in classes.values():E1.extend(kernel(es,1,p))
    return kernel(E1,2,p) if E1 else []


def enumerate_shapes(sec,Y,D):
    shapes=defaultdict(set);contents=[c for c in fw.contents_of_sector(sec,Y)
        if len(c)==Y and sum(L[0]+L[1] for L in c)==D]
    for j,cs in enumerate(contents):
        for w in multitraces_with_content(cs,Y):
            sg,v=so_canon(w)
            if sg:shapes[shape(v)].add(v)
        print(json.dumps({'content':j+1,'of':len(contents),'words':sum(map(len,shapes.values()))}),flush=True)
    return shapes


def build(sec,Y,D,p,out):
    start=time.monotonic();out=pathlib.Path(out);out.mkdir(exist_ok=True)
    shapes=enumerate_shapes(sec,Y,D)
    nwords=sum(map(len,shapes.values()));spin_groups=defaultdict(list);nflavour=0
    for index,sh in enumerate(sorted(shapes,key=lambda s:len(shapes[s]))):
        block=out/('flavour_%05d.pkl'%index)
        if block.exists():es=pickle.loads(block.read_bytes())
        else:
            es=hw_flavour(sorted(shapes[sh]),p);block.write_bytes(pickle.dumps(es,protocol=5))
        for e in es:spin_groups[shape(next(iter(e)),'spin')].append(e)
        nflavour+=len(es);del shapes[sh]
        if index%50==0:print(json.dumps({'flavour_block':index,'highest_weights':nflavour,'seconds':time.monotonic()-start}),flush=True)
    final=[]
    for index,(sh,es) in enumerate(sorted(spin_groups.items())):
        block=out/('spin_%05d.pkl'%index)
        if block.exists():rs=pickle.loads(block.read_bytes())
        else:rs=kernel(es,3,p);block.write_bytes(pickle.dumps(rs,protocol=5))
        final.extend(rs)
    for e in final:
        for i in (1,2,3):assert not raise_expr(e,i,p)
    (out/'basis.pkl').write_bytes(pickle.dumps(final,protocol=5))
    result=dict(sector=sec,Y=Y,derivatives=D,prime=p,formal_words=nwords,
        flavour_highest_weights=nflavour,flavour_and_spin_highest_weights=len(final),
        seconds=time.monotonic()-start,maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Formal highest-weight kernel modulo the stated prime. Finite-rank completeness still requires exact-dimension evaluation.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sector',required=True);ap.add_argument('--Y',type=int,required=True)
    ap.add_argument('--derivatives',type=int,required=True);ap.add_argument('--prime',type=int,default=1000003)
    ap.add_argument('--output',required=True);ap.add_argument('--seconds',type=int,default=600)
    a=ap.parse_args();signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('full highest-weight wall limit')))
    signal.alarm(a.seconds);build(tuple(map(int,a.sector.split(','))),a.Y,a.derivatives,a.prime,a.output)
