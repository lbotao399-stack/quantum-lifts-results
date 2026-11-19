"""Rank-three Pfaffian relations in the mixed terminal singlet, using exact integer contractions.

Each of ten homogeneous symplectic adjoints is a symmetric two-slot tensor
followed by the symplectic form. Antisymmetrize eight slots (the fourth
exterior power of a 2-form on six dimensions vanishes); pair the other twelve
with six oriented copies of the form. Record the relation in ordinary traces.
"""
import argparse, json, random
from collections import Counter
from itertools import combinations
from pathlib import Path
from flint import fmpz_mat
from orthogonal_relations import coinvariants
from mixed_cartan import compute as cartan_witness
from replay_orthogonal import exchanges, integer_cartan_test

HERE=Path(__file__).resolve().parent

def matchings(slots):
    if not slots:
        yield (),1
        return
    for j in range(1,len(slots)):
        for rest,s in matchings(slots[1:j]+slots[j+1:]):
            yield ((slots[0],slots[j]),)+rest,s*(-1)**(j-1)

def polynomial(species,slots):
    assert sorted(slots)==list(range(20)) and sorted(species)==[0,0,1,1,2,2,3,3,4,4]
    fixed=tuple(zip(slots[8::2],slots[9::2]));ans=Counter()
    for pf,pf_sign in matchings(tuple(slots[:8])):
        edge={}
        for i,j in pf+fixed:edge[i]=(j,1);edge[j]=(i,-1)
        unused=set(range(10));terms=[];coef=pf_sign
        while unused:
            start=2*min(unused);pos=start;cycle=[]
            while True:
                l=pos//2;unused.remove(l);cycle.append(l);pos,sg=edge[pos^1];coef*=sg
                if pos==start:break
            g=tuple(sum(species[i]==k for i in cycle) for k in range(5))
            if len(cycle)%2 or max(g[:2])>1:break
            terms.append((g,sorted((i for i in cycle if species[i]<2),key=species.__getitem__)))
        else:
            terms.sort(key=lambda t:(sum(t[0]),t[0]));word=tuple(g for g,ol in terms)
            if any(g==h and sum(g[:2])%2 for g,h in zip(word,word[1:])):continue
            # Reorder original matrix entries into the final ordered traces.
            coef*=(-1)**exchanges([i for _,odd in terms for i in odd]);ans[word]+=coef
    return ans

def make_case(rng):
    species=list(range(5))*2;rng.shuffle(species);order=list(range(10));rng.shuffle(order)
    # Choosing two slots in one symmetric tensor annihilates the Pfaffian identically.
    free=[2*i+rng.randrange(2) for i in order[:8]];rng.shuffle(free)
    if rng.randrange(2):
        # Leave a disconnected two-matrix trace in half of the trials to include products of relations.
        rem=[i^1 for i in free]+[2*order[8],2*order[9],2*order[8]+1,2*order[9]+1]
    else:
        rem=sorted(set(range(20))-set(free));rng.shuffle(rem)
    return species,free+rem

def compute(cases=None):
    basis,cov=coinvariants();rng=random.Random(202611190608);rs=[];save=[]
    for k in range(8000 if cases is None else len(cases)):
        s,l=make_case(rng) if cases is None else (cases[k]['species'],cases[k]['slots']);rel=polynomial(s,l)
        raw=[rel[w] for w in basis];row=list(map(int,(fmpz_mat([raw])*cov).tolist()[0]))
        if any(row) and fmpz_mat(rs+[row]).rank()>len(rs):
            save.append(dict(species=s,slots=l,formal_trace_coefficients=raw,coinvariant_coordinates=row));rs.append(row)
            if len(rs)==5:break
    assert len(rs)==5
    cs,dt=next((cs,d) for cs in combinations(range(6),5) for d in [int(fmpz_mat([[r[j] for j in cs] for r in rs]).det())] if d)
    raw=[row['formal_trace_coefficients'] for row in save];zero=integer_cartan_test(raw,basis,3);ctrl=integer_cartan_test(raw,basis,5)
    # The skew Cartan used by the common evaluator is complex symplectically conjugate to a diagonal one;
    # the total homogeneous degree ten multiplies every ordinary split-Cartan value by i^10=-1.
    assert zero==[0]*5 and all(ctrl);w=cartan_witness()['USp(6)']['witness_top_coefficient'];assert w==-64
    return dict(group='USp(6)',weight=[2]*5,free_slots=8,pfaffian_terms=105,formal_weight_space=len(basis),lowering_rank=218,formal_singlets=6,relations=save,rank=5,minor_columns=list(cs),minor_determinant=dt,rank3_cartan=zero,rank5_control=ctrl,rank3_witness=w,terminal_dimension=1)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--replay',type=Path);a=p.parse_args();assert __debug__
    old=json.loads(a.replay.read_bytes()) if a.replay else None;r=compute(old['relations'] if old else None)
    if old:assert old==r
    if a.output:a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(r['rank'],r['minor_columns'],r['minor_determinant'],r['rank3_witness'],r['rank5_control'])
