"""Exact integer flavour contractions and Q0; modular minors for the zero-jet singlet.

Source slot tags 0..8 label positions in alternating triples or tied phi/psi pairs.
Words store the seven BMN letters in the order of fortuity.bmn.words.QLET.
Only the final target (zero or one fermion) is evaluated by commuting arrays.
"""
import argparse, gzip, hashlib, json, math, random
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from pathlib import Path
import numpy as np
from flint import nmod_mat

HERE=Path(__file__).resolve().parent
P=1000003
SEED=202611187
PERMS=[(p,(-1)**sum(a>b for i,a in enumerate(p) for b in p[i+1:])) for p in permutations(range(3))]
Q={3:((1,1,2),(-1,2,1)),4:((1,2,0),(-1,0,2)),5:((1,0,1),(-1,1,0)),
   6:((1,0,3),(-1,3,0),(1,1,4),(-1,4,1),(1,2,5),(-1,5,2))}

def shapes(n,small=2):
    if not n:yield ()
    for s in range(small,n+1):
        for tail in shapes(n-s,s):yield (s,)+tail

@lru_cache(None)
def assignments(kind):
    # Triples carry epsilon_{ijk}; a pair is phi^i and psi_i in this order.
    triples,pairs,fs={'f':(3,0,1),'pp':(2,2,0),'p':(3,1,0),'scalar':(4,0,0)}[kind]
    out=[]
    for eps in product(PERMS,repeat=triples):
        coefficient=math.prod(s for _,s in eps)
        prefix=tuple(v for perm,_ in eps for v in perm)
        for matched in product(range(3),repeat=pairs):
            out.append((prefix+tuple(v for i in matched for v in (i,3+i))+(6,)*fs,coefficient))
    return out

@lru_cache(maxsize=140000)
def necklace(word):
    # At most one odd letter remains in these target words: no odd-odd exchange.
    assert sum(3<=i<=5 for i in word)<=1 and 6 not in word
    if len(word)<2:return 0,()
    rev=word[::-1]; sign=(-1)**len(word)
    a=min(word[i:]+word[:i] for i in range(len(word)))
    b=min(rev[i:]+rev[:i] for i in range(len(word)))
    if a==b and sign==-1:return 0,()
    return (sign,b) if b<a else (1,a)

@lru_cache(maxsize=60000)
def dtrace(tr):
    out=Counter();s=1
    for j,l in enumerate(tr):
        for u,b,c in Q.get(l,()):
            h,t=necklace(tr[:j]+(b,c)+tr[j+1:]);out[t]+=s*u*h
        if 3<=l<=5:s=-s
    return tuple((t,c) for t,c in out.items() if c)

def expression(recipe,apply_q=True):
    kind,layout=recipe;out=Counter()
    for species,co in assignments(kind):
        traces=tuple(tuple(species[pos] for pos in tr) for tr in layout)
        if not apply_q:
            original=[necklace(t) for t in traces];sg=math.prod(s for s,_ in original)
            if sg:out[tuple(sorted(w for _,w in original))]+=co*sg
            continue
        phase=1
        for j,t in enumerate(traces):
            if any(l>=3 for l in t):
                other=[necklace(tr) for k,tr in enumerate(traces) if k!=j];sg=math.prod(s for s,_ in other)
                if sg:
                    rest=tuple(w for _,w in other)
                    for w,c in dtrace(t):out[tuple(sorted(rest+(w,)))]+=co*phase*sg*c
            if sum(3<=l<=5 for l in t)%2:phase=-phase
    return {w:c for w,c in out.items() if c}

class Frame:
    def __init__(self,group,degree,points):
        self.group,self.degree,self.points=group,degree,points
        n=7 if group=='SO_odd' else 6; self.n=n
        basis=[]
        for i,j in combinations_with_replacement(range(n),2):
            e=np.zeros((n,n),dtype=np.int64)
            if group=='SO_odd':
                if i==j:continue
                e[i,j]=1;e[j,i]=-1
            else:
                e[i,j]=e[j,i]=1
                e=np.concatenate((-e[3:],e[:3]),axis=0)
            basis.append(e)
        self.basis=np.stack(basis); assert len(basis)==21
        self.rng=np.random.default_rng(SEED+degree)
        self.letters=(self.rng.integers(P,size=(3,points,21),dtype=np.int64)@self.basis.reshape(21,n*n)%P).reshape(3,points,n,n)
        self.width=points*(63 if degree else 1)
        # A trace with a single fermion has one flavour and exactly 21 free matrix components.
        self.word_cache=lru_cache(maxsize=7200)(self._word);self.trace_cache=lru_cache(maxsize=15000)(self._trace)
        self.multi_cache=lru_cache(maxsize=12000)(self._multi)
    def _word(self,w):
        if len(w)==1:return self.letters[w[0]]
        return self.word_cache(w[:-1])@self.letters[w[-1]]%P
    def _trace(self,w):
        odd=[i for i,l in enumerate(w) if l>=3]
        if not odd:return -1,self.word_cache(w).trace(axis1=-2,axis2=-1)%P
        assert len(odd)==1; i=odd[0];rot=w[i+1:]+w[:i]
        A=self.word_cache(rot).swapaxes(1,2).reshape(self.points,-1)
        return w[i]-3,(A@self.basis.reshape(21,-1).T)%P
    def _multi(self,w):
        scalar=1; odd=[]
        for t in w:
            l,v=self.trace_cache(t)
            if l<0:scalar=scalar*v%P
            else:odd.append((l,v))
        assert len(odd)==self.degree
        if not odd:return -1,scalar
        l,v=odd[0];return l,(scalar[:,None]*v%P if not isinstance(scalar,int) else v)
    def value(self,expr):
        n=sum(abs(c) for c in expr.values());assert n*(P-1)<2**63
        result=np.zeros((3,self.points,21) if self.degree else (self.points,),dtype=np.int64)
        for w,c in expr.items():
            l,v=self.multi_cache(w)
            if l<0:result+=c*v
            else:result[l]+=c*v
        return result.ravel()%P

class Echelon:
    def __init__(self,width):self.width=width;self.pivots=[];self.rows=[];self.source=[];self.raw=[]
    def add(self,value,recipe):
        row=value.copy()
        for j,r in zip(self.pivots,self.rows):
            if row[j]:row=(row-row[j]*r)%P
        nz=np.flatnonzero(row)
        if not len(nz):return
        pivot=int(nz[0]);row=row*pow(int(row[pivot]),-1,P)%P
        self.rows.append(row);self.pivots.append(pivot);self.source.append(recipe);self.raw.append(value)
    def determinant(self):
        if not self.rows:return 1
        return int(nmod_mat(np.stack(self.raw)[:,self.pivots].tolist(),P).det())

def decode(r):return r[0],tuple(tuple(t) for t in r[1])

def dimensions(group):
    entries=json.loads((HERE/'charge_maps.json').read_text()); ds=next(e['by_length_and_fermion_degree'] for e in entries if e['group']==group and e['mu']==[0,0,4,4,4]);return sum(ds['11'].values()),sum(ds['12'].values())

def main(args):
    groups=['SO_odd','Sp'];maps=[(g,l) for g in groups for l in (10,11)]
    frames={(g,l):Frame(g,11-l,8 if l==10 else 100) for g,l in maps}
    caps={(g,l):(dimensions(g)[1]-1 if l==11 else dimensions(g)[0]-dimensions(g)[1]+1) for g,l in maps}
    spaces={k:Echelon(frames[k].width) for k in maps}
    def save(attempts):
        out={}
        for (g,l),s in spaces.items():
            det=s.determinant();assert det
            out[g+':'+str(l)]={'group':g,'source_length':l,'target_length':l+1,'points':frames[g,l].points,'rank_lower_bound':len(s.rows),'rank_ceiling':caps[g,l], 'minor_determinant':det,'columns':s.pivots,'source_contractions':s.source}
        data={'prime':P,'frame_seed':SEED,'character_input_sha256':hashlib.sha256((HERE/'charge_maps.json').read_bytes()).hexdigest(),'attempts':attempts,'maps':out}
        path=args.output;temp=path.with_suffix('.tmp');temp.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0));temp.replace(path)
    if args.replay:
        saved=json.loads(gzip.decompress(args.replay.read_bytes()));assert saved['prime']==P and saved['frame_seed']==SEED and saved['character_input_sha256']==hashlib.sha256((HERE/'charge_maps.json').read_bytes()).hexdigest()
        for g,l in maps:
            row=saved['maps'][g+':'+str(l)]; assert row['points']==frames[g,l].points and row['rank_ceiling']==caps[g,l]
            for r in row['source_contractions']:
                r=decode(r);assert sorted(i for t in r[1] for i in t)==list(range(l));assert r[0] in (('f','pp') if l==10 else ('p',));v=frames[g,l].value(expression(r));spaces[g,l].add(v,r)
            assert row['columns']==spaces[g,l].pivots and row['rank_lower_bound']==len(spaces[g,l].rows) and row['minor_determinant']==spaces[g,l].determinant();print(g,l,row['rank_lower_bound'],row['minor_determinant'],flush=True)
        return
    rng=random.Random(SEED);choices={l:list(shapes(l)) for l in (10,11)}
    for i in range(args.draws):
        if all(len(spaces[k].rows)==caps[k] for k in maps):break
        l=10 if i%3!=2 else 11
        if all(len(spaces[g,l].rows)==caps[g,l] for g in groups):continue
        kind=('f' if i%6==0 else 'pp') if l==10 else 'p'; positions=list(range(l));rng.shuffle(positions);shape=rng.choice(choices[l]);layout=[];j=0
        for size in shape:layout.append(tuple(positions[j:j+size]));j+=size
        recipe=(kind,tuple(layout));e=expression(recipe)
        if e:
            for g in groups:
                if len(spaces[g,l].rows)<caps[g,l]:spaces[g,l].add(frames[g,l].value(e),recipe)
        if i%150==0:
            print(i,{str(k):len(v.rows) for k,v in spaces.items()},flush=True);save(i+1)
    save(i+1);print('final',{str(k):(len(v.rows),v.determinant()) for k,v in spaces.items()},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--draws',type=int,default=4500);p.add_argument('--output',type=Path,default=HERE/'scalar_minors.json.gz');p.add_argument('--replay',type=Path);main(p.parse_args())
