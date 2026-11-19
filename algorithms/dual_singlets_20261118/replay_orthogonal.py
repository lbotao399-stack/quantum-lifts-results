"""Independent directed-loop and exterior-sign replay of all five integer relations and the Cartan tests."""
import argparse, json
from collections import Counter
from itertools import permutations,combinations
from pathlib import Path
from flint import fmpz_mat
from terminal_commutators import words,parity
from mixed_cartan import compute as witness
HERE=Path(__file__).resolve().parent

def exchanges(labels):return sum(a>b for i,a in enumerate(labels) for b in labels[i+1:])%2

def pairing(species, slots, rank):
    """Evaluate the orthogonal determinant after quotienting by the supercommutators; rank = determinant order."""
    ans=Counter();edges=list(zip(slots[2*rank::2],slots[2*rank+1::2]));N=len(species)
    for perm in permutations(range(rank)):
        pairs=edges+[(slots[i],slots[rank+perm[i]]) for i in range(rank)]; links=dict(pairs);links.update((y,x) for x,y in pairs)
        unused=set(range(N));trace=[];odd_order=[];neg=exchanges(perm)
        while unused:
            # Start at the largest remaining matrix and at its second slot: the opposite traversal convention to the C++ program.
            pos=2*max(unused)+1;start=pos;loop=[]
            while True:
                mat=pos//2;unused.remove(mat);loop.append(mat);neg^=(pos%2);pos=links[pos^1];
                if pos==start:break
            tally=tuple(sum(species[mat]==s for mat in loop) for s in range(5))
            if len(loop)%2 or max(tally[:2])>1:break
            oo=[mat for mat in loop if species[mat]<2];odd_order+=oo;neg^=exchanges([species[m] for m in oo]);trace.append(tally)
        else:
            neg^=exchanges(odd_order);deg_key=lambda g:(sum(g),g)
            odd_traces=[deg_key(g) for g in trace if parity(g)]
            if len(odd_traces)!=len(set(odd_traces)):continue
            neg^=exchanges(odd_traces);ans[tuple(sorted(trace,key=deg_key))]+=(-1)**neg
    return ans

def cartan(b,colour_dim):
    # Three independent 2x2 skew blocks; two odd coefficients of each species along the first two blocks.
    z={0:1}
    for g in b:
        gen=Counter();
        for i in range(min(3,colour_dim//2)):
            if any(g[2+j] and j!=i for j in range(3)) or (i>=2 and any(g[:2])):continue
            mask=(g[0]<<i)|(g[1]<<(i+2));gen[mask]+=2*(-1)**(sum(g)//2)
        temp=Counter()
        for a,val in z.items():
            for m,c in gen.items():
                if a&m:continue
                flip=sum(((m%2**j).bit_count()) for j in range(4) if a&(1<<j));temp[a|m]+=val*c*(-1)**flip
        z=temp
    return z.get(15,0)

def integer_cartan_test(raw,basis,rank):
    # Distinct integer weights on each Cartan slot and 4 arbitrary odd exterior coefficients distinguish the stable traces.
    def value(g):
        z=Counter()
        for c in range(rank):
            sc=(c+2,(2*c+1)**2,c*c-3*c+1);x=2*(-1)**(sum(g)//2)
            for k,p in enumerate(g[2:]):x*=sc[k]**p
            term={0:x}
            for t in range(2):
                if not g[t]:continue
                odd=((t+1)*c*c+2, (c+1)**3+3*t,1+(c+t)**2,2*c*c*c-c*t+4);n=Counter()
                for a,v in term.items():
                    for i,vv in enumerate(odd):
                        if not a&(1<<i):n[a|1<<i]+=v*vv*(-1)**((a>>(i+1)).bit_count())
                term=n
            z.update(term)
        return z
    gens={g:value(g) for w in basis for g in w};f=[]
    for w in basis:
        z={0:1}
        for g in w:
            q=Counter()
            for a,v in z.items():
                for m,c in gens[g].items():
                    if not a&m:q[a|m]+=v*c*(-1)**sum((a>>(i+1)).bit_count() for i in range(4) if m&(1<<i))
            z=q
        f.append(z.get(15,0))
    return [sum(v*y for v,y in zip(row,f)) for row in raw]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,default=HERE/'orthogonal_relations.json');ap.add_argument('--output',type=Path);a=ap.parse_args();assert __debug__
    d=json.loads(a.input.read_bytes());basis=words(tuple(d['weight']));raw=[];reduced=[]
    for rec in d['independent_relations']:
        species,slots=rec['species'],rec['slots'];assert sorted(slots)==list(range(20)) and sorted(species)==[0,0,1,1,2,2,3,3,4,4]
        vals=pairing(species,slots,8);row=[vals[x] for x in basis];assert row==rec['formal_trace_coefficients'];raw.append(row);reduced.append(rec['coinvariant_coordinates']);
    mat=fmpz_mat(reduced);idx,det=next((cols,int(m.det())) for cols in combinations(range(6),5) for m in [fmpz_mat([[r[i] for i in cols] for r in reduced])] if m.det());assert mat.rank()==5
    er=integer_cartan_test(raw,basis,3);big=integer_cartan_test(raw,basis,5);assert er==[0]*5 and any(big)
    w=witness();out=dict(permutations_per_relation=40320,integer_coefficients_recomputed=sum(map(len,raw)),coinvariant_rows=reduced,minor_columns=idx,minor_determinant=det,rank3_cartan=er,rank5_cartan_control=big,cartan_witness={k:v['witness_top_coefficient'] for k,v in w.items()});print(json.dumps(out,indent=2));
    if a.output:a.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
