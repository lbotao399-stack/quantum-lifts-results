"""Rational interval ranks in adjacent full-letter cochain degrees of SO and USp."""
import argparse
import hashlib
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from charge_maps import _multiplicity, content_patterns, multiplicity_parts
from rank_four_s5 import Q1d_multi
from terminal_interval_ranks import sample


def dimensions(weight, length, size, group):
    if group in ("SU+", "SU-", "U+", "U-"):
        from unitary_characters import parity_parts
        assert size >= 2
    result = Counter()
    for ds, spin_sign in ((0, 1), (1, -1)):
        for perm in permutations(range(3)):
            shifted = (weight[0]+ds, weight[1]-ds) + tuple(
                weight[2+j]+perm[j]-j for j in range(3))
            sign = spin_sign * (-1)**sum(i>j for k,i in enumerate(perm) for j in perm[k+1:])
            for bosons, fermions, count in content_patterns(shifted, length):
                if group in ("SU+", "SU-", "U+", "U-"):
                    multiplicity = parity_parts(bosons, fermions, size, group.startswith("U"))[group.endswith("-")]
                elif group=="SO":
                    multiplicity = _multiplicity(bosons, fermions, 0, size)
                else:
                    assert group=="Sp"
                    multiplicity = multiplicity_parts(bosons, fermions, size//2)
                result[sum(fermions)] += sign*count*multiplicity
    assert all(n>=0 for n in result.values())
    return {k:n for k,n in sorted(result.items()) if n}


def highest_weights(weight, length, parity=None):
    if parity is not None:
        from unitary_trace_weights import highest_weights as unitary_weights
        return unitary_weights(weight,length,parity)
    rows = defaultdict(list)
    for block in m.words_by_shape(length, weight).values():
        for row in m.singlets(block):
            degrees = {m.nferm(w) for w in row}
            assert len(degrees)==1
            rows[degrees.pop()].append(row)
    return dict(rows)


def multiply_rows(rows, words, values):
    """Return integral evaluations, switching to Python integers before overflow."""
    lookup = {w:i for i,w in enumerate(words)}
    bound = max((sum(abs(n) for n in row.values()) for row in rows), default=0)
    maximum = max(abs(int(values.min(initial=0))), int(values.max(initial=0)))
    if bound*max(maximum,1) < 2**63:
        entries, columns, starts = [], [], [0]
        for row in rows:
            for w,n in row.items():
                entries.append(n)
                columns.append(lookup[w])
            starts.append(len(entries))
        sparse = csr_matrix((np.asarray(entries,np.int64),columns,starts),
                            shape=(len(rows),len(words)))
        return sparse @ values
    return np.asarray([[sum(n*int(values[lookup[w],j]) for w,n in row.items())
                        for j in range(values.shape[1])] for row in rows],dtype=object)


def frame_values(words, degree, size, points, seed, group, alphabet, dense_modular):
    if not dense_modular:
        return sample(words,degree,size,points,seed,group,alphabet)
    assert group.startswith("SU") and size < m.P
    assert size * math.comb(degree+1,(degree+1)//2) * (m.P-1)**2 < 2**63
    frame=m.Frame(size,points,degree+1,[degree],alphabet,seed=seed,p=m.P,group="SU")
    return frame.eval_words(words,dtype=np.int64)


def full_frame(rows, dim, degree, size, group, seed, alphabet, dense_modular=False):
    words = sorted({w for row in rows for w in row})
    chunks = []
    batches = []
    columns = 0
    # Take enough samples to cover the character dimension before counting
    # additional batches for linearly dependent evaluations.
    first_full_batch = math.ceil(dim / (6 * (degree + 1)))
    piv = []
    next_check = dim
    for k in range(first_full_batch + 128):
        points = min(6,max(2,math.ceil((next_check+4-columns)/(degree+1))))
        chosen_seed = seed+4099*k
        values = frame_values(words,degree,size,points,chosen_seed,group,alphabet,dense_modular)
        chunks.append(multiply_rows(rows,words,values) % m.P)
        batches.append([points,chosen_seed])
        columns += values.shape[1]
        if columns < next_check:
            continue
        mod = nmod_mat(np.column_stack(chunks).tolist(),m.P)
        piv = m.pivotcols(mod)
        assert len(piv)<=dim
        if len(piv)==dim:
            # Restrict the transpose to the columns already selected.
            selected = nmod_mat([[int(mod[i,j]) for i in range(mod.nrows())] for j in piv],m.P)
            ids = m.pivotcols(selected)
            determinant = int(nmod_mat([[int(mod[i,j]) for j in piv] for i in ids],m.P).det())
            assert len(ids)==dim and determinant
            record=dict(batches=batches,columns=piv,rows=ids,determinant_mod_prime=determinant,
                        dimension=dim,formal=len(rows),alphabet=alphabet)
            if dense_modular:
                record["dense_modular"]=True
            return [rows[i] for i in ids],record
        next_check = columns + max(degree+1, dim-len(piv))
    raise ValueError(("incomplete finite-colour frame",degree,dim,len(piv)))


def basis(rows, dims, group, size, seed, alphabet, dense_modular=False):
    chosen = []
    frames = {}
    for deg,dim in dims.items():
        values,frames[deg] = full_frame(rows[deg],dim,deg,size,group,seed+317*deg,alphabet,dense_modular)
        chosen.extend(values)
    return chosen,frames


@lru_cache(30000)
def quantum_word(word, order, size, group):
    if order==0:
        if group.startswith("SU"):
            return tuple((w,n) for w,n in m.fw.Q_multi(word).items() if n)
        return m.differential(word)
    name="SU" if group.startswith("SU") else ("SO_odd" if group=="SO" else "Sp")
    if order==1:
        if group!="Sp":
            op=Q1d_multi
        else:
            from q1d_sp_exact import Q1d_multi as op
            name = "Sp"
    elif order==2:
        if group!="Sp":
            from q2d_exact import Q2d_multi as op
        else:
            from q2d_sp_exact import Q2d_multi as op
            name="Sp"
    else:
        assert order==3
        from q3d_sp_exact import Q3d_multi as op
    # A fan of order n consumes all three flavour labels n times.
    flavours=sorted((len(l[2]) for tr in word for l in tr),reverse=True)
    if sum(flavours[:order+1]) < 3*order:
        return ()
    acted=op(word,size,name)
    if group.startswith("SU"):
        return tuple((w,n) for w,n in acted.items() if n)
    out=Counter()
    for term,n in acted.items():
        sig,w=m.so_canon(term)
        if sig:
            out[w]+=sig*n
    return tuple((w,n) for w,n in out.items() if n)


def action(row,order,size,group):
    result=Counter()
    for w,n in row.items():
        for v,c in quantum_word(w,order,size,group):
            result[v]+=n*c
    return {w:n for w,n in result.items() if n}


def integer_actions(rows, max_order, size, group):
    maps=[[] for _ in range(max_order+1)]
    for row in rows:
        images=[action(row,k,size,group) for k in range(max_order+1)]
        denominator=math.lcm(*(Fraction(n).denominator for img in images for n in img.values()))
        for block,img in zip(maps,images):
            block.append({w:int(n*denominator) for w,n in img.items()})
    return maps


def evaluate(rows, frames, size, group):
    pieces=[]
    for degree,record in frames.items():
        restricted=[{w:c for w,c in row.items() if m.nferm(w)==degree} for row in rows]
        words=sorted({w for row in restricted for w in row})
        if not words:
            pieces.append(np.zeros((len(rows),record["dimension"]),np.int64))
            continue
        evals=[]
        previous=0
        for points,seed in record["batches"]:
            width=points*(degree+1)
            cols=[k-previous for k in record["columns"] if previous<=k<previous+width]
            if cols:
                vals=frame_values(words,degree,size,points,seed,group,record["alphabet"],
                                  record.get("dense_modular",False))[:,cols]
                evals.append(multiply_rows(restricted,words,vals))
            previous+=width
        pieces.append(np.column_stack(evals))
    if not pieces:
        return np.zeros((len(rows),0),np.int64)
    return np.column_stack(pieces)


def interval_ranks(blocks, source_dims, target_dims):
    ranks={}
    pairings=[]
    for width in range(1,len(source_dims)+1):
        for lo in range(len(source_dims)-width+1):
            hi=lo+width-1
            row_offsets=np.cumsum([0]+source_dims[lo:hi+1]).tolist()
            col_offsets=np.cumsum([0]+target_dims[lo:hi+1]).tolist()
            matrix=fmpz_mat(row_offsets[-1],col_offsets[-1])
            for (a,b),block in blocks.items():
                if lo<=a<=b<=hi:
                    for i,row in enumerate(block.tolist(),row_offsets[a-lo]):
                        for j,value in enumerate(row,col_offsets[b-lo]):
                            if value:
                                matrix[i,j]=int(value)
            r=matrix.rank()
            ranks[lo,hi]=r
            count=(r-ranks.get((lo+1,hi),0)-ranks.get((lo,hi-1),0)
                   +ranks.get((lo+1,hi-1),0))
            assert count>=0
            pairings.append([lo,hi,r,count])
    return sorted(pairings)


def run(flavours,spin,deficit,specs,work,salt=0,screen=False,save_grades=False,dense_unitary=False):
    work.mkdir(parents=True,exist_ok=True)
    weights=[(t+spin,t)+tuple(f-t for f in flavours) for t in range(min(flavours)+1)]
    target_lengths=[sum(w)-deficit for w in weights]
    # A predecessor has total Taylor degree minus length equal to deficit+1.
    max_order=min((deficit+2)//2,len(weights)-1)
    assert max_order<=3
    if dense_unitary and (not screen or any(not g.startswith("SU") for g,_ in specs)):
        raise ValueError("Dense unitary frames require an SU finite-field screen.")
    modes={g if g.startswith("SU") else "orthosymplectic" for g,n in specs}
    if len(modes)>1:
        for mode in sorted(modes):
            same=[(g,n) for g,n in specs if (g if g.startswith("SU") else "orthosymplectic")==mode]
            run(flavours,spin,deficit,same,work,salt,screen,save_grades,dense_unitary)
        return
    mode=modes.pop()
    parity=1 if mode=="SU+" else (-1 if mode=="SU-" else None)
    formal=[]
    for w,y in zip(weights,target_lengths):
        formal.append((highest_weights(w,y-1,parity),highest_weights(w,y,parity)))
        print("formal",w,y,[sum(map(len,v.values())) for v in formal[-1]],flush=True)
    for group,size in specs:
        frames=[]; source_dims=[]; target_dims=[]; maps=[]; blocks={}
        stem="band"+str(deficit)+"_"+"_".join(map(str,flavours))+"_s"+str(spin)+"_"+group.lower()+str(size)
        if save_grades:
            checkpoints=work/(stem+"_grades")
            checkpoints.mkdir(exist_ok=True)
        for t,(w,y,(parents,targets)) in enumerate(zip(weights,target_lengths,formal)):
            seed=9121801+53*size+3907*t+salt
            ds=dimensions(w,y-1,size,group)
            dt=dimensions(w,y,size,group)
            src_alphabet=sorted({l for content in m.contents(w,y-1) for l in content})
            tgt_alphabet=sorted({l for content in m.contents(w,y) for l in content})
            src,src_frames=basis(parents,ds,group,size,seed,src_alphabet,dense_unitary)
            _,tgt_frames=basis(targets,dt,group,size,seed+100000,tgt_alphabet,dense_unitary)
            frames.append(dict(weight=w,length=y,source=src_frames,target=tgt_frames))
            source_dims.append(sum(ds.values()))
            target_dims.append(sum(dt.values()))
            maps.append(integer_actions(src,max_order,size,group))
            for k in range(1,max_order+1):
                if t+k>=len(weights):
                    assert all(not row for row in maps[-1][k])
            for k in range(min(t,max_order)+1):
                a=t-k
                block=evaluate(maps[a][k],tgt_frames,size,group)
                assert tuple(block.shape)==(source_dims[a],target_dims[t])
                blocks[a,t]=block
            if save_grades:
                # Integer coordinates may exceed int64 for longer words.
                active={f"{a}_{b}":v for (a,b),v in blocks.items() if b==t}
                if all(v.dtype!=object for v in active.values()):
                    temp=checkpoints/f"grade{t}.part"
                    with temp.open("wb") as file:
                        np.savez_compressed(file,**active)
                    temp.replace(checkpoints/f"grade{t}.npz")
                    m.write(checkpoints/f"frame{t}.json",frames[-1])
            print(group,size,"t",t,"dimensions",source_dims[-1],target_dims[-1],flush=True)
        if screen:
            from filtered_rank import filtered_ranks
            ranks=filtered_ranks(blocks,source_dims,target_dims,m.P)
        else:
            ranks=interval_ranks(blocks,source_dims,target_dims)
        path=work/("band"+str(deficit)+"_"+"_".join(map(str,flavours))+"_s"+str(spin)+"_"+group.lower()+str(size)+".json")
        cert=dict(group=group,colours=size,flavours=flavours,spin=spin,target_deficit=deficit,
                  prime_for_frame_selection=m.P,maximum_direct_order=max_order,
                  ranks_by_order=[sum(n for a,b,r,n in ranks if b-a==j) for j in range(len(weights))],
                  interval_ranks_and_pairings=ranks,frames=frames,source_dimensions=source_dims,
                  target_dimensions=target_dims,block_hashes={f"{a}-{b}":hashlib.sha256(
                      m.compact(v.tolist())).hexdigest() for (a,b),v in blocks.items()})
        if screen:
            cert["rank_field"] = f"F_{m.P}"
            cert["screen_only"] = True
        m.write(path,cert)
        print(path.name,cert["ranks_by_order"],flush=True)
        quantum_word.cache_clear()


def settings():
    here=Path(__file__).resolve().parent
    sys.path.insert(0,str(here.parent/"so7_cartan_target_20261118"/"native"))
    sys.path.insert(1,str(here.parent/"su3_j34_20261118"/"third_action_rules"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flavours",type=int,nargs=3,required=True)
    p.add_argument("--spin",type=int,default=0)
    p.add_argument("--deficit",type=int,choices=range(6),required=True)
    p.add_argument("--groups",nargs="+",required=True,help="SO7 Sp6 SU+4 SU-5, using defining dimensions")
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--salt",type=int,default=0)
    p.add_argument("--screen",action="store_true",help="Use finite-field cancellation for triage")
    p.add_argument("--dense-unitary",action="store_true",help="Use dense traceless unitary frames in a screen")
    p.add_argument("--save-grades",action="store_true",help="Save numeric blocks as each target grade finishes")
    a=p.parse_args()
    assert __debug__ and a.spin>=0 and a.flavours==sorted(a.flavours,reverse=True) and min(a.flavours)>=0
    settings()
    specs=[(s[:3],int(s[3:] or "3")) if s.startswith(("SU+","SU-")) else (s[:2],int(s[2:]))
           for s in a.groups]
    assert all(g in ("SO","Sp","SU+","SU-") and n >= (2 if g.startswith("SU") else 3)
               and (g!="Sp" or n%2==0) for g,n in specs)
    run(a.flavours,a.spin,a.deficit,specs,a.work,a.salt,a.screen,a.save_grades,a.dense_unitary)
