"""Native colour frames and the selected polarized trace recipes."""
from collections import defaultdict
import gzip
import json
import numpy as np
from scipy.sparse import csr_matrix
import exceptional_rows
from flag_rows import ADJOINT, MARK, polynomial


def matrices(folder):
    fundamental=json.loads(gzip.decompress((folder/'representation26.json.gz').read_bytes()))
    rep=np.zeros((52,26,26),dtype=np.int64)
    for i,row in enumerate(fundamental['matrices']):
        for j,k,v in row:rep[i,j,k]=v
    algebra=json.loads(gzip.decompress((folder/'inputs.json.gz').read_bytes()))
    adj=np.zeros((52,52,52),dtype=np.int64)
    for i,j,row in algebra['brackets']:
        for k,v in row:
            adj[i,k,j]=v
            adj[j,k,i]=-v
    for i in range(52):
        for j in range(i):
            assert np.array_equal(rep[i]@rep[j]-rep[j]@rep[i],(adj[i,:,j]@rep.reshape(52,-1)).reshape(26,26))
    return rep,adj


def inputs(mode,seeds,two_reps):
    monomials={}
    row_indices=[];data=[];ends=[0]
    for s in seeds:
        p=polynomial(mode,s,two_reps)
        for word,coeff in p.items():
            ix=monomials.setdefault(word,len(monomials))
            row_indices.append(ix);data.append(coeff)
        ends.append(len(data))
    matrix=csr_matrix((np.array(data,dtype=np.int64),row_indices,ends),shape=(len(seeds),len(monomials)))
    atoms=sorted(set(v for m in monomials for v in m))
    ix={k:i for i,k in enumerate(atoms)}
    factors=defaultdict(list)
    for m,j in monomials.items():factors[len(m)].append((j,[ix[k] for k in m]))
    factors={k:(np.array([j for j,v in rows]),np.array([v for j,v in rows])) for k,rows in factors.items()}
    groups=defaultdict(list)
    for pos,x in enumerate(atoms):
        rep,tail=divmod(x,ADJOINT)
        pair,tail=divmod(tail,MARK)
        n=(tail.bit_length()-1)//2
        code=tail-4**n
        assert 0<=code<4**n and rep<2 and pair<=6
        groups[rep,pair,n].append((pos,code))
    groups={k:np.asarray(v,dtype=np.int64).T for k,v in groups.items()}
    return matrix,groups,factors,len(atoms)


def pivots(rows,prime):
    basis={}; selected=[];columns=[];det=1
    for i,r in enumerate(rows):
        v=r.copy()
        for p in sorted(basis):
            if v[p]:v=(v-v[p]*basis[p])%prime
        nz=np.flatnonzero(v)
        if not len(nz):continue
        p=int(nz[0]);pivot=int(v[p])
        det=det*pivot%prime
        basis[p]=v*pow(pivot,-1,prime)%prime
        selected.append(i);columns.append(p)
    return selected,columns,det


def prepare(mode,recipes):
    terms={};idx=[];co=[];ptr=[0]
    for k,s in recipes:
        polynomial=exceptional_rows.polynomial(mode,s,k)
        for word,v in polynomial.items():
            ix=terms.setdefault((k,*word),len(terms))
            idx.append(ix);co.append(v)
        ptr.append(len(idx))
    matrix=csr_matrix((np.asarray(co,dtype=np.int64),idx,ptr),shape=(len(recipes),len(terms)))
    closed=sorted({c for k,a,b,c in terms if c})
    allopen=sorted({a for k,x,y,c in terms for a in (x,y)})
    op={j:i for i,j in enumerate(allopen)}
    cf={j:i for i,j in enumerate([0]+closed)}
    by_mode=defaultdict(list)
    for (k,a,b,c),pos in terms.items():by_mode[k].append([pos,op[a],op[b],cf[c]])
    by_mode={k:np.asarray(v,dtype=np.int64).T for k,v in by_mode.items()}
    levels=defaultdict(set)
    for a in allopen:
        v=a
        while v>1:
            n=(v.bit_length()-1)//4
            levels[n].add(v);v//=16
        assert v==1
    previous={1:0};layers=[]
    for n in range(1,max(levels,default=0)+1):
        layer=sorted(levels[n])
        annotated=[(i,op[a]) for i,a in enumerate(layer) if a in op]
        use,dest=zip(*annotated) if annotated else ((),())
        layers.append((np.array([previous[a//16] for a in layer]),np.array([a%16 for a in layer]),
                       np.array(use,dtype=np.int64),np.array(dest,dtype=np.int64)))
        previous={a:i for i,a in enumerate(layer)}
    groups=defaultdict(list)
    for id in closed:
        pair,z=divmod(id,MARK)
        n=(z.bit_length()-1)//2
        groups[0,pair,n].append((cf[id],z-4**n))
    groups={k:np.asarray(v,dtype=np.int64).T for k,v in groups.items()}
    return matrix,by_mode,layers,groups,len(allopen),len(closed)+1
