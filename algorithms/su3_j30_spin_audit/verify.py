"""Exact dotted-spin branching and equivariance for the SU(3) R-007 third target.
Only the four pinned external rank/metadata records are read; no large rank is replayed here.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import argparse, json, sys

here = Path(__file__).resolve().parent
rt = here.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
sys.path.insert(0,str(rt))
from su3_charge_dimensions_v1 import parity_parts
def sector(mu, length):
    n,m,*theta=mu
    return (2*n+sum(theta)-length,2*m+sum(theta)-length)+tuple(n+m+2*x for x in theta)
from fortuity.full.letters import Q_letter, letters_up_to

@lru_cache(None)
def patterns(weight, length):
    # Enumerate occupied species in order. Skip branches even before choosing any
    # species if the remaining theta degrees exceed the number of letters left.
    alphabet=sorted(letters_up_to(weight[0],weight[1]),key=lambda l:(len(l[2]),l),reverse=True)
    species=[]
    for a,b,s in alphabet:
        w=(a,b)+tuple(int(i in s) for i in range(3))
        if all(x<=y for x,y in zip(w,weight)) and sum(w)<=sum(weight)-length+1:
            species.append((w, len(s)%2==0))
    out=Counter()
    def visit(start,rem,left,bos,fer):
        if not left:
            if not any(rem):out[tuple(sorted(bos)),tuple(sorted(fer))]+=1
            return
        if left>sum(rem) or max(rem[2:])>left:return
        for i in range(start,len(species)):
            w,odd=species[i]
            if any(x>y for x,y in zip(w,rem)):continue
            rr=rem
            for n in range(1,left+1):
                rr=tuple(y-x for x,y in zip(w,rr))
                if min(rr)<0:break
                visit(i+1,rr,left-n,bos if odd else bos+(n,),fer+(n,) if odd else fer)
    visit(0,weight,length,(),())
    return out


def fixed_spin(weights,length):
    co=[Counter(),Counter()];rho=(2,1,0)
    for perm in permutations(range(3)):
        w=weights[:2]+tuple(weights[i+2]+rho[i]-rho[perm[i]] for i in range(3))
        if min(w)<0:continue
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for (bos,fer),multi in patterns(w,length).items():
            for n,x in enumerate(parity_parts(bos,fer)):co[n][sum(fer)]+=sign*multi*x
    assert all(v>=0 for row in co for v in row.values())
    return [{k:v for k,v in sorted(row.items()) if v} for row in co]


def spin(l,direction):
    a,b,s=l
    if direction==1:return ((a+1,b-1,s),a+1) if b else (None,0)
    return ((a-1,b+1,s),b+1) if a else (None,0)


def equivariance():
    # Compare on free noncommuting letters, before imposing any trace relation.
    # The spin operators are even derivations, hence there is no Koszul sign.
    count=0
    for l in letters_up_to(6,6):
        if sum(l[:2])>6:continue
        for direction in (1,-1):
            v,f=spin(l,direction); delta=Counter()
            if f:
                for c,u,w in Q_letter(v):delta[u,w]+=f*c
            for c,u,w in Q_letter(l):
                uu,f=spin(u,direction)
                if f:delta[uu,w]-=c*f
                ww,f=spin(w,direction)
                if f:delta[u,ww]-=c*f
            assert all(v==0 for v in delta.values()),(l,direction,delta)
            count+=1
    return count


def run():
    weight=(3,3,4,1,1); length=8; assert sector(weight,length)==(4,4,14,8,8)
    fixed=fixed_spin(weight,length);shift=fixed_spin((4,2,4,1,1),length)
    sing=[{d:a.get(d,0)-b.get(d,0) for d in sorted(a.keys()|b.keys())} for a,b in zip(fixed,shift)]
    rec={};dig={}
    for f in sorted((here/'upstream_receipts').glob('*.json')):
        data=f.read_bytes();rec[f.stem]=json.loads(data);dig[f.name]=sha256(data).hexdigest()
    assert len(rec)==4
    inside=rec['r007_o3_in_codd_meta'];outside=rec['r007_o3_out_codd_meta']
    for m,y,s in [(inside,8,[5,5,14,8,8]),(outside,9,[4,4,14,8,8])]:
        assert m['SEC']==s and m['N']==3 and m['Y']==y and m['p']==67108859 and m['csec']==-1
    assert fixed[1]=={int(k):v for k,v in inside['target_m_d'].items()}
    assert all(v>=0 for row in sing for v in row.values()) and sum(map(sum,[row.values() for row in sing]))==52490
    dim=sum(fixed[1].values());assert dim==outside['tp_m']==131048
    ranks=[]
    for tag,m in [('in',inside),('out',outside)]:
        r=rec[f'r007_o3_{tag}_codd_rank_gpu']; assert r['case']==m['case'] and (r['nB'],r['ncols'])==(m['nB'],m['ncols'])
        ranks.append(r['rank_B'])
    assert ranks==[67869,63179] and dim==sum(ranks)
    return dict(upstream_pin='4256f207986029540060fb44dd494746828f65ce',source='R-007',order=3,length=length,
                target_doubled_charges=sector(weight,length),derivative_and_flavour_degrees=weight,
                fermion_degrees_are_cochain_statistics=True,
                components=[dict(C=p,spin_weight_zero=f,spin_weight_one=h,spin_singlet=s,singlet_total=sum(s.values())) for p,f,h,s in zip((1,-1),fixed,shift,sing)],
                checked_free_letter_spin_commutators=equivariance(),upstream_metadata_and_rank_receipt_sha256=dig,
                upstream_rank_prime=inside['p'],upstream_incoming_outgoing_lower_bounds=ranks,
                independently_replayed_large_rank_matrices=False,odd_full_weight_cohomology_upper_bound=dim-sum(ranks),
                odd_singlet_cohomology_dimension_using_upstream_ranks=0)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args();assert __debug__
    blob=json.dumps(run(),indent=2)+'\n';print(blob,end='')
    if args.output:args.output.write_text(blob)
