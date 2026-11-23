"""Exact comparison after stripping Q0-closed pure-scalar trace factors."""
from inverse_support import *
import hashlib,math

A=(0,0,(0,))

def reduced_blocks(defect):
    groups=defaultdict(dict)
    for w,c in defect.items():groups[block(w)][w]=c
    records=[]
    for key,part in sorted(groups.items()):
        reduced={tuple(tr for tr in w if not all(l==A for l in tr)):c for w,c in part.items()}
        divisor=math.gcd(*reduced.values())
        if reduced[min(reduced)]<0:divisor=-divisor
        normalized={w:c//divisor for w,c in reduced.items()}
        serialized=json.dumps([[c,pf.m.encode(w)] for w,c in sorted(normalized.items())],separators=(',',':'))
        records.append((key,normalized,divisor,hashlib.sha256(serialized.encode()).hexdigest()))
    return records

def main():
    start=time.monotonic()
    old=pickle.loads((HERE/'known_defect_sp5.pkl').read_bytes())
    new=pickle.loads((HERE/'small_sp6/candidate_defect.pkl').read_bytes())
    previous=reduced_blocks(old);current=reduced_blocks(new)
    reference={digest:(key,row,divisor) for key,row,divisor,digest in previous}
    unique={};shared=[]
    for key,row,divisor,digest in current:
        if digest in unique:assert row==unique[digest][1]
        else:unique[digest]=(key,row,divisor)
        if digest in reference:
            old_key,old_row,old_divisor=reference[digest];assert row==old_row
            shared.append(dict(current_trace_weights=key,previous_trace_weights=old_key,
                               current_integer_factor=divisor,previous_integer_factor=old_divisor,
                               terms=len(row),sha256=digest))
    result=dict(Sp5_blocks=len(previous),Sp6_blocks=len(current),
        Sp5_distinct_stripped_polynomials=len(reference),Sp6_distinct_stripped_polynomials=len(unique),
        Sp6_blocks_reusing_a_Sp5_stripped_polynomial=len(shared),
        exact_coefficient_comparisons=True,closed_spectator='Tr((0,0,(0,))^(even length))',
        scopes='Linear inverse-Q tasks only; not an equivalence or independence test for physical families',
        seconds=time.monotonic()-start,shared_blocks=shared)
    (HERE/'spectator_reuse.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='shared_blocks'}),flush=True)

if __name__=='__main__':main()
