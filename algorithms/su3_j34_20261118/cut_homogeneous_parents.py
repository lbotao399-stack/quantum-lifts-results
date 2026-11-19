"""Cut the entire source support of two homogeneous targets and solve the affected rows.

The charge-conjugation-even frame and its integer split tables are those of
A8 5664529. The residual of the new starting vector is computed from scratch
with the upstream Q_multi. A separator is accepted only after it annihilates
every allowed column over Q. No column is discarded on numerical tolerance.
"""
import argparse
import hashlib
from fractions import Fraction
from math import lcm
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, hstack, load_npz, vstack
import shared_success_coupling as shared
import complete_joint_minor as dixon
j=shared.j; s=j.s; P=j.P

def content(w):return bytes(sorted(b''.join(w))).hex()
def setup(census, original, graph):
    side,fw,back,cost,_,_,_=shared.graph(census);inc=set().union(*(back[t] for t in shared.FAIL)); old=s.read(census);passed=set(old['exact_zero']);base,oldhash=shared.scan.old.vector(original);assert oldhash==old['vector_sha']
    removed={g:sum(content(w)==g for w,v in base) for g in sorted(inc)};removed={g:n for g,n in removed.items() if n}; six=set.union(*(fw[g]&passed for g in removed));assert len(inc)==19 and sorted(removed.values())==[16,55,57] and len(six)==6 and set(side)<six<=set(graph['targets']);assert set.union(*(set(back[t]) for t in graph['targets']))==set(graph['groups']);return dict(base),fw,inc,passed,removed,sorted(six)

def prepare(rhs,original,census,prior,meta):
    graph=s.read(prior/'graph.json'); receipt=s.read(meta)['receipt'];assert receipt['job']=='5664529' and j.side.canonical_digest(graph)==receipt['graph'];orig,fw,inc,passed,removed,six=setup(census,original,graph);targets=graph['targets']; assert set(targets)<=passed
    # Start from the known solution. Apart from forbidden groups, zero only
    # the groups whose entire previous successful support lies in these 12 rows.
    # This makes the rational numbers far smaller (here all retained words are integral).
    stable={g for g in graph['groups'] if not ((fw[g]&passed)-set(targets))};base={w:v for w,v in orig.items() if content(w) not in (inc|stable)};assert all(v.denominator==1 for v in base.values());b={t:{} for t in targets}; counts={t:0 for t in shared.FAIL};h=hashlib.sha256();count=0
    for w,co in s.stream_rows(rhs,h):
        t=content(w);count+=1
        if t in b:b[t][w]=Fraction(co)
        if t in counts:counts[t]+=1
    assert count==517468 and h.hexdigest()==s.INPUT and not any(counts.values())
    for w,co in base.items():
        for u,m in j.side.replay.Q_multi(s.expanded(w)).items():
            v=s.packed(u);t=content(v);assert t not in counts
            if t in b:j.l.p.axpy(b[t],{v:m},-co)
    allgs=sorted(set(graph['groups'])-inc); conf=sorted(stable-inc); assert (len(conf),len(allgs))==(11,32)
    return graph,receipt,orig,base,b,allgs,conf,dict(forbidden=sorted(inc),removed_words_from_2724=removed,affected_original_successes=six,additional_reset_groups=sorted(stable-inc),targets=targets,homogeneous_rhs_word_counts=counts,starting_words=len(base),unaffected_previous_successes_if_confined=len(passed-set(targets)))

def algebra(prior,graph,receipt,gs,b,exe,out):
    ts=graph['targets'];dims=graph['dimensions'];fr={k:s.read(prior/'frames'/(k+'.json.gz')) for k in gs+ts}; assert all(j.side.canonical_digest(z)==receipt['frames'][k]['digest'] and len(z['words'])==z['C_even_dimension']==len(set(z['polynomial_keys']))==dims[k] and z['determinant_mod_prime']==receipt['frames'][k]['det']!=0 for k,z in fr.items())
    D=lcm(*(a.denominator for rr in b.values() for a in rr.values()));y=[];blocks=[];keys=[]
    for t in ts:
        k0=[int(v) for v in fr[t]['polynomial_keys']];jx={k:i for i,k in enumerate(k0)};rhs=[0]*len(k0);words=sorted(b[t]);
        for a in range(0,len(words),160):
            take=words[a:a+160]
            for w,rr in zip(take,shared.E.poly_iter(take,exe,out),strict=True):
                c=int(b[t][w]*D)
                for k,v in rr:
                    if k in jx:rhs[jx[k]]+=c*v
        m=s.read(prior/(t+'.map.json'));mat=load_npz(prior/(t+'.npz')).tocsr();offset=0;parts={}
        for g in m['groups']:
            assert m['frames'][g]==receipt['frames'][g]['digest'];parts[g]=mat[offset:offset+dims[g]];offset+=dims[g]
        assert mat.shape==(offset,len(k0)) and m['frames'][t]==receipt['frames'][t]['digest'] and mat.nnz==receipt['maps'][t]['entries'];blocks.append(vstack([parts[g] if g in parts else csr_matrix((dims[g],len(k0)),dtype='int64') for g in gs],format='csr'));y.extend(rhs);keys.extend((t,k) for k in k0)
    A=hstack(blocks,format='csr');pool=[(g,w) for g in gs for w in fr[g]['words']];assert A.shape==(len(pool),len(keys))==(9291,5475);return A,y,D,pool,keys

def verify_rect(B,nums,y,den):return j.multiply(B,nums)==[den*v for v in y]
def certificate(A,y,report):
    oi=np.argsort(np.diff(A.indptr),kind='stable');oj=np.argsort(np.diff(A.T.tocsr().indptr),kind='stable');C=[int(oj[k]) for k in j.piv(A[oi][:,oj])];R=[int(oi[k]) for k in j.piv(A[oi][:,C].T.tocsr())];assert len(C)==len(R)>0;B=A[R].T.tocsr();M=B[C];ref=j.prime_matrix(M);det=int(ref.det());assert det
    from flint import nmod_mat
    x=[int(a) for a in ref.solve(nmod_mat(len(C),1,[y[k]%P for k in C],P)).entries()];del ref;nz=[k for k,v in enumerate(j.multiply(B,x)) if (v-y[k])%P];r=dict(prime=P,image_rank_mod_prime=len(C),augmented_rank_mod_prime=len(C)+bool(nz),minor_determinant_mod_prime=det);report(dict(stage='modular ranks',**r))
    if nz:
        # Pick a failing row with a sparse candidate functional; prove it on the full rectangle.
        ns=np.diff(A.T.tocsr().indptr);k=min(nz,key=lambda i:(ns[i],i));d=A[:,C].tocsr();v=[int(x) for x in A[:,k].toarray().ravel()];info,ans=dixon.solve(d,v,R,report)
        if ans is None:return dict(**r,kind='uncertified',reason=info),None
        nums,den=ans;assert verify_rect(d,nums,v,den);ff={c:-Fraction(a,den) for c,a in zip(C,nums) if a};ff[k]=Fraction(1);pair=sum(y[k]*a for k,a in ff.items());assert pair!=0;return dict(**r,kind='exact_separator',integer_rhs_pairing=str(pair),digits=info['digits']),ff
    info,ans=dixon.solve(B,y,C,report)
    if ans is None:return dict(**r,kind='uncertified',reason=info),None
    nums,den=ans;assert verify_rect(B,nums,y,den);return dict(**r,kind='exact_solution',digits=info['digits']),{i:Fraction(a,den) for i,a in zip(R,nums) if a}

def run(rhs,original,census,prior,meta,exe,out):
    out.mkdir(exist_ok=True,parents=True);graph,receipt,orig,base,b,gs,conf,meta0=prepare(rhs,original,census,prior,meta);s.write(out/'setup.json',meta0);A,y,D,pool,keys=algebra(prior,graph,receipt,gs,b,exe,out);record={}
    for name,allowed in [('preserve_296',conf),('all_nonincident',gs)]:
        idx=[i for i,(g,w) in enumerate(pool) if g in allowed];history=[]
        def progress(v):
            history.append(v);s.write(out/(name+'.progress.json'),history);print(name,v,flush=True)
        res,values=certificate(A[idx],y,progress);res.update(source_words=len(idx),source_contents=allowed,denominator=str(D),target_coordinates=len(keys));doc=dict(kind=res['kind'],forbidden=meta0['forbidden'],allowed=allowed);record[name]=res
        if res['kind']=='exact_separator':doc.update(functional=[[keys[k][0],keys[k][1],str(v)] for k,v in sorted(values.items())],pairing=str(Fraction(res['integer_rhs_pairing'])/D))
        if res['kind']=='exact_solution':
            sol=dict(base);terms=[]
            for i,a in sorted(values.items()):
                g,w=pool[idx[i]];terms.append([g,w,str(a/D)]);j.l.p.axpy(sol,shared.E.projected(s.unpack_hex(w)),a/D)
            assert len(sol) and all(content(w) not in meta0['forbidden'] for w in sol);j.l.save_vector(out/(name+'.vector.json.gz'),sol);doc.update(frame_solution=terms,vector_words=len(sol));res['vector_words']=len(sol)
        s.write(out/(name+'.certificate.json.gz'),doc);s.write(out/'results.json',record)
    return record

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('rhs','original','census','prior','meta','exe','out'):ap.add_argument('--'+k,type=Path,required=True)
    run(**vars(ap.parse_args()))
