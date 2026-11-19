"""Couple two failed evaluations to five complete, previously solved equations.

A restricted source includes all three parents that do not touch another success.
An expanded source includes every incident parent of cost <=3. Omitted coordinates
of the two failures are checked separately if a projected primitive exists.
The solver uses modular pivots only to propose an exact, rectangular certificate.
"""
import argparse, hashlib, json
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import lcm, gcd
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, vstack, hstack, load_npz
import joint_content_correction as j
import centralizer_scan as scan
import replay_centralizer_sweep as rc
s, E, f, P = j.s, j.E, j.l.f, j.P
FAIL = ('0000000000010405060712','0000000000010405060d10')

def graph(census):
    bad, back, cost, doc = scan.graph(census)
    passes = set(doc['exact_zero']); fw=defaultdict(set)
    for t, gs in back.items():
        for g in gs: fw[g].add(t)
    new = {t:[g for g in back[t] if cost[g]==3] for t in FAIL}
    assert all(len(v)==1 for v in new.values())
    side = sorted(set.union(*(fw[x[0]] & passes for x in new.values())))
    assert set(FAIL)<=bad and len(side)==5 and len(set.intersection(*(fw[x[0]] & passes for x in new.values())))==1
    hit = set(FAIL)|set(side); inc=sorted(g for g,ts in fw.items() if ts & hit)
    scopes={'preserve_303':sorted(g for g in inc if (fw[g]&passes)<=set(side)),
            'cost_at_most_three':sorted(g for g in inc if cost[g]<=3), 'all_36':inc}
    dims={g:E.parity_frame_dims(bytes.fromhex(g))[0] for g in inc+sorted(hit)}
    assert [(len(gs),sum(dims[g] for g in gs)) for gs in scopes.values()]==[(3,1983),(16,14218),(36,31903)]
    assert sum(dims[t] for t in side)==2244 and all(dims[t]==4012 for t in FAIL)
    return side,fw,back,cost,new,scopes,dims

def evaluator(t,witness):
    bits=[i for i in range(7) if witness['mask']>>i&1]; n=len(bits); ids=set(bytes.fromhex(t)); assert n==sum(f.ODD[i] for i in bytes.fromhex(t))
    M={(i,h):tuple(map(tuple,witness['matrices'][f'{i:02x}:{bits[h] if h>=0 else -1}'])) for i in ids for h in (range(n) if f.ODD[i] else [-1])}
    @lru_cache(150000)
    def trace(tr,order):
        it=iter(order); mats=[M[i,next(it) if f.ODD[i] else -1] for i in tr]; x=mats[0]
        for A in mats[1:]: x=scan.ft.mm(x,A)
        return sum(x[i][i] for i in range(3))
    @lru_cache(120000)
    def ev(w):
        v=0
        for perm in permutations(range(n)):
            a=scan.ft.sign(perm); k=0
            for tr in w:
                m=sum(f.ODD[i] for i in tr);a*=trace(tr,perm[k:k+m]);k+=m
            assert k==n;v+=a
        return v
    return ev,M,n

def parents_and_side(side,fw,scopes,dims,prior,receipt):
    meta=s.read(receipt);meta=meta.get('receipt',meta); assert meta['job']=='5664529'
    gs=[g for g in scopes['cost_at_most_three'] if fw[g]&set(side)]
    assert len(gs)==12
    frames={k:s.read(prior/'frames'/(k+'.json.gz')) for k in gs+side}
    for k,z in frames.items():
        a=meta['frames'][k];assert j.side.canonical_digest(z)==a['digest'] and z['content']==k and z['determinant_mod_prime']==a['det'] and len(z['words'])==len(z['polynomial_keys'])==dims[k]==a['dim']
    blocks=[];keys=[]
    for t in side:
        old=s.read(prior/(t+'.map.json')); matrix=load_npz(prior/(t+'.npz')).tocsr();offs={};idx=0
        for g in old['groups']:
            d=meta['frames'][g]['dim'];offs[g]=idx;idx+=d
        assert matrix.shape==(idx,dims[t]) and len(old['groups'])==len(offs) and old['frames'][t]==meta['frames'][t]['digest']
        for g in gs:
            if t in fw[g]:assert g in offs and old['frames'][g]==meta['frames'][g]['digest']
        blocks.append(vstack([matrix[offs[g]:offs[g]+dims[g]] if t in fw[g] else csr_matrix((dims[g],dims[t]),dtype='int64') for g in gs],format='csr')); keys.extend((t,int(k)) for k in frames[t]['polynomial_keys'])
    return gs,frames,hstack(blocks,format='csr'),keys

def solve_cert(A,y):
    # A: source words by observed coordinates. A rational primitive is checked
    # on every coordinate; an exact dual functional on every source word.
    oi=np.argsort(np.diff(A.indptr),kind='stable');oj=np.argsort(np.diff(A.T.tocsr().indptr),kind='stable'); C=[int(oj[k]) for k in j.piv(A[oi][:,oj])]; R=[int(oi[i]) for i in j.piv(A[oi][:,C].T.tocsr())]; r=len(C);assert r==len(R)>0
    B=A[R].T.tocsr();M=B[C].tocsr();b=[y[k] for k in C]; fld=j.prime_matrix(M);det=int(fld.det());assert det
    from flint import nmod_mat
    x=[int(a) for a in fld.solve(nmod_mat(r,1,[v%P for v in b],P)).entries()]; non=[k for k,(a,v) in enumerate(zip(j.multiply(B,x),y)) if (a-v)%P];rank=dict(prime=P,image_rank_mod_prime=r,augmented_rank_mod_prime=r+bool(non),minor_determinant_mod_prime=det)
    def integer(xs):
        den=lcm(*(v.denominator for v in xs));return [int(v*den) for v in xs],den
    def dot(sparse,x):return j.multiply(sparse,x)
    if non:
        k=min(non,key=lambda k:A[:,k].nnz); D=A[:,C].tocsr();v=A[:,k].toarray().ravel().tolist(); h=[v[i] for i in R];
        def accept(xs):
            nums,den=integer(xs);return dot(D,nums)==[den*a for a in v] and den*y[k]!=sum(p*n for p,n in zip(b,nums))
        xs=j.exact_solve(M.T.tocsr(),h,160,accept)
        if xs is None:return dict(**rank,kind='unreconstructed'),None
        eq={k:Fraction(1)};eq.update({c:-v for c,v in zip(C,xs) if v});pair=sum(y[i]*w for i,w in eq.items());assert pair
        return dict(**rank,kind='exact_separator',pairing_integer_rhs=str(pair)),eq
    def accept(xs):
        nums,den=integer(xs);return dot(B,nums)==[den*v for v in y]
    xs=j.exact_solve(M,b,160,accept)
    if xs is None:return dict(**rank,kind='unreconstructed'),None
    return dict(**rank,kind='exact_solution'),{i:v for i,v in zip(R,xs) if v}

def run(census,sweep,rat,shared,rhs,prior,receipt,out):
    out.mkdir(exist_ok=True,parents=True);side,fw,back,cost,new,scopes,dims=graph(census); gs,fr,side_matrix,sidekeys=parents_and_side(side,fw,scopes,dims,prior,receipt)
    z=s.read(sweep);r=s.read(rat); hv,rows=scan.old.residuals(rhs,shared,side+list(FAIL)); assert z['source_sha']==hv==r['source'] and z['rhs_sha']==r['rhs']==s.INPUT
    pool=[(g,w) for g in gs for w in fr[g]['words']];proj=[];true_values=[];gcds=[];records={};checked={}
    for t in FAIL:
        item=rc.test(z['targets'][t],2);ev,M,n=evaluator(t,item['witness']);value=sum((co*ev(w) for w,co in rows[t].items()),Fraction());assert value==Fraction(r['replayed'][t]['rational']);raw=[];base=new[t][0];trial=0
        for g,w in pool:
            if g!=base:raw.append(0);continue
            d={}
            for u,k in E.projected(s.unpack_hex(w)).items():
                for v,c in f.q(u).items():
                    if bytes(sorted(b''.join(v))).hex()==t:d[v]=d.get(v,0)+c*k
            if trial%83==0:
                for v in list(d)[:2]:assert ev(v)==scan.ft.exterior_value(v,M,n)
            trial+=1;raw.append(sum(k*ev(v) for v,k in d.items()))
        assert trial==dims[base] and any(raw)
        for g in scopes['cost_at_most_three']:
            if g in back[t] and g!=base:
                k,a,b,coef=back[t][g]; assert cost[g]<=2 and g in item['allowed']
                for (i,xi),X in M.items():
                    if i==a:
                        for (k,yi),Y in M.items():
                            if k==b:assert scan.ft.mm(X,Y)==scan.ft.mm(Y,X)
        div=gcd(*raw);assert div>0; vals=[x//div for x in raw]; true_values.append(-value/div);gcds.append(str(div)); proj.append(vals);checked[t]=sum(bool(x) for x in vals)
    # Integer normalization of the two *equations* only, retaining the source
    # coefficients in the archived unnormalized (1+C) convention. Source matrix
    # entries can exceed int64: rescale each integral witness by its exact gcd.
    mx=max(abs(v) for p in proj for v in p);assert mx<2**63
    A=hstack([side_matrix,csr_matrix(np.array(proj,dtype=np.int64).T)],format='csr'); keys=sidekeys+[(t,'cost_two_witness') for t in FAIL]; D=lcm(*(v.denominator for v in true_values));y=[0]*len(sidekeys)+[int(D*v) for v in true_values];assert D%P
    for name in ('preserve_303','cost_at_most_three'):
        scope=set(scopes[name]);indices=[i for i,(g,w) in enumerate(pool) if g in scope];hidden=sorted(scope-set(gs));assert all(not (fw[g]&set(side)) and cost[g]<=2 for g in hidden);result,coeff=solve_cert(A[indices],y);result.update(source_dimension=sum(dims[g] for g in scope),exactly_zero_columns=sum(dims[g] for g in hidden),stored_source_columns=len(indices),target_coordinates=len(keys),denominator=str(D));doc=dict(kind=result['kind'],source=hv,scope=name,targets=side+list(FAIL),side_frames={t:j.side.canonical_digest(fr[t]) for t in side},source_frames={g:j.side.canonical_digest(fr[g]) for g in sorted(scope&set(gs))},additional_zero_parent_contents=hidden,gcds=gcds,nonzero_witness_columns=checked)
        if result['kind']=='exact_separator':
            doc['functional']=[[keys[i][0],str(keys[i][1]),str(v)] for i,v in sorted(coeff.items())];doc['pairing_with_scaled_residual']=str(Fraction(result['pairing_integer_rhs'])/D)
        if result['kind']=='exact_solution':
            doc['frame_solution']=[[pool[indices[i]][0],pool[indices[i]][1],str(v/D)] for i,v in sorted(coeff.items())];delta={}
            for g,w,co in doc['frame_solution']:j.l.p.axpy(delta,E.projected(s.unpack_hex(w)),Fraction(co))
            v,hh=scan.old.vector(shared);assert hh==hv;base=dict(v);j.l.p.axpy(base,delta,1);j.l.save_vector(out/(name+'.vector.json.gz'),base);doc['vector_words']=len(base)
        s.write(out/(name+'.certificate.json.gz'),doc);records[name]=result;s.write(out/'results.json',records);print(name,result,flush=True)
    return records

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data',type=Path,default=s.HERE/'joint_70');ap.add_argument('--shared',type=Path);ap.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');ap.add_argument('--prior',type=Path,required=True);ap.add_argument('--receipt',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert __debug__;p=a.data
    run(p/'census_2724.json',p/'centralizer_sweep.json.gz',p/'centralizer_rational.json',a.shared or p/'vector_shared.json.gz',a.rhs,a.prior,a.receipt or p/'a8_5664529.json',a.out)
