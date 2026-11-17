"""Q1 obstruction in the complete reflection-even full-letter target space."""
from highest import *
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm as full_odd
from fortuity.full.wordenum import multitraces_with_content
from fortuity.full.q1 import Q1_expr
from fortuity.full.evalg import Frame as FullFrame,nferm
from exact_dims import _multiplicity
CONV={0:(1,(0,0,(0,))),1:(1,(0,0,(1,))),2:(1,(0,0,(2,))),
      3:(1,(0,0,(1,2))),4:(-1,(0,0,(0,2))),5:(1,(0,0,(0,1))),6:(-1,(0,0,(0,1,2)))}


def to_full(expr):
    out=defaultdict(int)
    for (paths,traces),c in expr.items():
        assert not paths
        ts=[]
        for w in traces:
            t=[]
            for l in w:
                s,L=CONV[l];c*=s;t.append(L)
            ts.append(tuple(t))
        fw.add_term(out,tuple(ts),c)
    return {w:c for w,c in out.items() if c}


def so_canon(m):
    ts=[];total=1
    for w in m:
        s,c=fw.canon_word(w);sr,cr=fw.canon_word(w[::-1])
        d=sum(full_odd(L) for L in w);sr*=(-1)**(len(w)+d*(d-1)//2)
        if not s or (c==cr and s!=sr):return 0,None
        if cr<c:c,s=cr,sr
        total*=s;ts.append(c)
    s,v=fw.canon_multi(ts)
    return total*s,v


def canonical(expr):
    out=defaultdict(int)
    for w,c in expr.items():
        s,v=so_canon(w)
        if s:out[v]+=s*c
    return {w:c for w,c in out.items() if c}


def all_words(sec,Y):
    words=set();contents=[c for c in fw.contents_of_sector(sec,Y)
                          if len(c)==Y and sum(L[0]+L[1] for L in c)==2]
    for content in contents:
        for w in multitraces_with_content(content,Y):
            s,v=so_canon(w)
            if s:words.add(v)
    return sorted(words),contents


def exact_dimension(contents):
    total=0
    from collections import Counter
    for content in contents:
        cs=Counter(content)
        bos=tuple(sorted([v for L,v in cs.items() if not full_odd(L)],reverse=True))
        fer=tuple(sorted([v for L,v in cs.items() if full_odd(L)],reverse=True))
        total+=_multiplicity(bos,fer,0,8)
    return total


def main():
    started=time.monotonic()
    O=deserialize(json.loads((ROOT/'so8_j18_representative.json').read_text())['terms'])
    full=to_full(O);q1=canonical(Q1_expr(full,8,'SO_odd'))
    assert {fw.sector_of(w) for w in q1}=={(0,0,6,6,6)}
    assert {sum(len(w) for w in m) for m in q1}=={8}
    assert not canonical(fw.Q_expr(q1))
    print('Q1 terms',len(q1),'seconds',time.monotonic()-started,flush=True)
    S,sc=all_words((1,1,6,6,6),7);T,tc=all_words((0,0,6,6,6),8)
    dim=exact_dimension(tc);assert dim==848
    idx={w:i for i,w in enumerate(T)}
    Q=[canonical(fw.Q_expr({w:1})) for w in S]
    assert all(w in idx for e in Q+[q1] for w in e)
    rows=[];cols=[];vals=[]
    for i,e in enumerate(Q):
        for w,c in e.items():rows.append(i);cols.append(idx[w]);vals.append(c)
    mat=csr_matrix((np.array(vals,dtype=np.int64),(rows,cols)),shape=(len(S),len(T)))
    c=np.array([q1.get(w,0) for w in T],dtype=np.int64)
    letters=sorted({L for w in T for t in w for L in t});degrees=sorted({nferm(w) for w in T})
    assert degrees==[2]
    print('source words',len(S),'target words',len(T),'dimension',dim,flush=True)
    result=[]
    for p in (1000003,1000033):
        for seed in (831,1831):
            # Six Grassmann components/point, at least exact_dimension+16 columns.
            fr=FullFrame(8,ceil((dim+16)/6),4,degrees,letters,seed=seed,p=p,group='SO_odd')
            X=fr.eval_words(T,dtype=np.int64);r=rank_mod(X,p);assert r==dim
            B=np.asarray(mat@X)%p
            value=np.zeros(X.shape[1],dtype=np.int64)
            for i,coef in enumerate(c):
                if coef:value=(value+int(coef%p)*X[i])%p
            rin=rank_mod(B,p);added=rank_mod(np.vstack([B,value]),p)
            row={'prime':p,'seed':seed,'exact_target_dimension':dim,'target_evaluation_rank':r,
                 'classical_image_rank':rin,'image_plus_Q1O_rank':added,'obstruction_rank':added-rin,
                 'evaluation_columns':fr.ncols}
            result.append(row);print(json.dumps(row),flush=True)
    out={'source_sector':[1,3,3,3],'full_target_sector':[0,0,6,6,6],
         'target_letters':8,'target_derivatives':2,'source_words':len(S),'target_words':len(T),
         'checks':result,'Q1_normalization':'Pinned upstream q1.py; K=1 including its absorbed common factors. Only the nonvanishing obstruction is tested; no coefficient is identified with a physical anomalous dimension.',
         'seconds':time.monotonic()-started}
    (ROOT/'so8_j18_quantum.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
