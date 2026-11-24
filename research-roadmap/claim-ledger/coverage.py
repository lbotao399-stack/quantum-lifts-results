#!/usr/bin/env python3
"""Enumerate the full finite charge envelope, without a signed-index screen.

Every row denotes an actual native cochain slot, which may have zero gauge
invariant dimension. Character dimensions are deliberately not guessed.
The letter alphabet is the seven nonempty subsets of three odd coordinates,
with arbitrary nonnegative dotted jet degrees, and the empty subset with
strictly positive total dotted degree. Its parity is 1-|subset| modulo two.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def twice_physical(native):
    a,b,c,u,v,length = native
    total = a+b+c
    return (total+length+2*(u+v), 2*a-total+length,
            2*b-total+length, 2*c-total+length,
            total-length+u+v, u-v)


def level(native):
    return 2*sum(native[:3])+3*sum(native[3:5])


def realizable(native):
    a,b,c,u,v,length = native
    if min(native)<0:
        return False
    if length == 0:
        return a+b+c+u+v == 0
    return max(a,b,c) <= length <= a+b+c+u+v


def slots(cutoff):
    yield (0,0,0,0,0,0)
    for a in range(cutoff//2+1):
        for b in range(cutoff//2-a+1):
            for c in range(cutoff//2-a-b+1):
                remaining=cutoff-2*(a+b+c)
                for u in range(remaining//3+1):
                    for v in range((remaining-3*u)//3+1):
                        total=a+b+c+u+v
                        for length in range(max(1,a,b,c),total+1):
                            yield (a,b,c,u,v,length)


def actions(native):
    a,b,c,u,v,length = native
    outgoing=[]
    for order in range(min(a,b,c)+1):
        target=(a-order,b-order,c-order,u+order,v+order,length+1-order)
        if realizable(target):
            outgoing.append((order,target))
    incoming=[]
    for order in range(min(u,v)+1):
        source=(a+order,b+order,c+order,u-order,v-order,length-1+order)
        if realizable(source):
            incoming.append((order,source))
    return incoming,outgoing


def corrections(native):
    a,b,c,u,v,length=native
    for order in range(min(a,b,c)+1):
        corrected=(a-order,b-order,c-order,u+order,v+order,length-order)
        if realizable(corrected):
            yield order,corrected


def alphabet(cutoff):
    for a,b,c in itertools.product(range(2),repeat=3):
        for u in range(cutoff//3+1):
            for v in range(cutoff//3+1):
                native=(a,b,c,u,v,1)
                if 0<level(native)<=cutoff:
                    yield native


def independent_letter_check(cutoff):
    """Exhaustive additive alphabet check of the closed-form envelope.

    This is a charge-realizability check before finite-colour exterior and
    invariant relations; all such zero spaces remain visible in the census.
    """
    reached={(0,0,0,0,0,0)}
    letters=list(alphabet(cutoff))
    for current_level in range(cutoff+1):
        current=[x for x in reached if level(x)==current_level]
        for row in current:
            for letter in letters:
                if current_level+level(letter)<=cutoff:
                    reached.add(tuple(x+y for x,y in zip(row,letter)))
    expected=set(slots(cutoff))
    assert reached==expected, {'missing':list(expected-reached)[:4], 'extra':list(reached-expected)[:4]}
    return len(expected)


def gzip_json_lines(path, rows):
    # Stable compressed bytes: no local path or wall-clock timestamp in header.
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9) as zipped:
            for row in rows:
                zipped.write((json.dumps(row,separators=(',',':'))+'\n').encode())


def build(output):
    output.mkdir(parents=True,exist_ok=True)
    scope=json.loads((ROOT/'research-roadmap/scope.json').read_text())
    cutoff=scope['core_maximum_target_level_inclusive']
    groups=scope['core_full_jet_groups']
    rows=sorted(slots(cutoff),key=lambda x:(level(x),x))
    native_ids={x:i for i,x in enumerate(rows)}
    assert len(native_ids)==len(rows)
    physical=defaultdict(list)
    totals=Counter()
    max_in=max_out=0
    for index,row in enumerate(rows):
        assert realizable(row)
        assert level(row)<=cutoff
        physical[twice_physical(row)].append(index)
        totals[level(row)]+=1
        incoming,outgoing=actions(row)
        for order,other in incoming+outgoing:
            assert other in native_ids and level(other)==level(row)
        for order,other in corrections(row):
            assert other in native_ids and twice_physical(row)==twice_physical(other)
        max_in=max(max_in,max((r for r,_ in incoming),default=0))
        max_out=max(max_out,max((r for r,_ in outgoing),default=0))
    physical_rows=sorted(physical)
    physical_ids={row:i for i,row in enumerate(physical_rows)}
    native_file=output/'native_slots.jsonl.gz'
    physical_file=output/'physical_blocks.jsonl.gz'
    gzip_json_lines(native_file, (list(row)+[physical_ids[twice_physical(row)]] for row in rows))
    gzip_json_lines(physical_file, (list(row)+[physical[row]] for row in physical_rows))
    # Explicit expansion is performed, checked, and hashed, without duplicating
    # the identical kinematic envelope eleven times in the repository.
    expansion=hashlib.sha256()
    expanded_count=0
    for group_index,group in enumerate(groups):
        for native_index,row in enumerate(rows):
            encoded=json.dumps([group_index,native_index,*row],separators=(',',':'))+'\n'
            expansion.update(encoded.encode()); expanded_count+=1
    assert expanded_count==len(groups)*len(rows)
    check=independent_letter_check(12)
    result={
        'schema_version':1,'scope_version':scope['scope_version'],
        'source_commit':'19cd8718b752b04104d6ff0becee42d580b27ca3',
        'source':'supplement/supplement.tex equations E-sector, E-charges and E-grading',
        'groups':groups,'maximum_target_level_inclusive':cutoff,
        'native_row_fields':['theta1','theta2','theta3','n_plus','n_minus','letter_length','physical_row_id'],
        'physical_row_fields':['twice_energy','twice_R1','twice_R2','twice_R3','twice_J_left','twice_J_right','native_row_ids'],
        'native_slots_per_group':len(rows),'physical_charge_blocks_per_group':len(physical_rows),
        'expanded_native_group_slots':expanded_count,
        'expanded_physical_group_blocks':len(groups)*len(physical_rows),
        'expanded_native_group_rows_sha256':expansion.hexdigest(),
        'native_slots_by_level':dict(sorted(totals.items())),
        'maximum_kinematically_allowed_incoming_order':max_in,
        'maximum_kinematically_allowed_outgoing_order':max_out,
        'all_classical_quantum_and_correction_neighbours_included':True,
        'independent_alphabet_check_level':12,'independent_alphabet_check_slots':check,
        'vacuum_included':True,
        'zero_spaces':'Length-one spaces vanish for each listed simple gauge algebra; additional zero spaces require invariant-character calculations.',
        'unknown_dimensions':'All other gauge-invariant dimensions remain explicitly uncomputed in this manifest.',
        'coverage_percentage':None,
        'counting_unit':'native cochain slot and complete physical charge block are distinct; neither count is a number of BPS classes',
        'storage':'Two explicit row lists plus the eleven-element group list. The complete Cartesian expansion has been executed and hashed; --expand emits every indexed group/native row.',
        'sha256':{file.name:hashlib.sha256(file.read_bytes()).hexdigest() for file in (native_file,physical_file)}
    }
    (output/'coverage_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def expand(output):
    manifest=json.loads((output/'coverage_manifest.json').read_text())
    with gzip.open(output/'native_slots.jsonl.gz','rt') as stream:
        native=[json.loads(line) for line in stream]
    for group_index,_group in enumerate(manifest['groups']):
        for native_index,row in enumerate(native):
            print(json.dumps([group_index,native_index,*row[:6]],separators=(',',':')))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'coverage')
    parser.add_argument('--expand',action='store_true')
    args=parser.parse_args()
    (expand if args.expand else build)(args.output)
