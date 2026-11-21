"""Verify the two precomputed binary descendants and record their ordered profile hashes."""
from pathlib import Path
import argparse, hashlib, json
from inputs import HERE, P, PARENT, PLAN, PACKAGE, ROOT, sha
BINARY={'S20':('36f2a8adda8623e95168cd6cfb5ceb7f44303e0520fd51132935574a2268dd3a',2770228), 'Tp21':('95e719290895aff772d020459c03e96096da985ca0e86cfbf2730a2c93a4b21e',3411015)}
def desc(s):return tuple(sorted(k[0]*5+k[1] for k,n in s['profile'] for _ in range(n)))
def filename(s):return 'p_'+'_'.join(map(str,desc(s)))+'.bin'
def receipt(folder):
    assert __debug__;read=lambda f:json.loads(f.read_bytes());pin=read(folder/'input.json');raw=folder/'controls/controls.json';ctl=read(raw);closure=read(folder/'letter_closure.json');assert sha(ROOT/'data/SU3_j30_parent_p1000033/parent.npy')==pin['parent_sha256']==PARENT and sha(HERE.parent/'su3_j30_classical_ranks/joint_plan.json')==pin['plan_sha256']==PLAN
    assert list(pin['recipe_sha256'])==[sha(PACKAGE/'second_basis/shard0'/f'targets{i:02d}/independent_seeds.pkl') for i in range(37)] and sha(folder/'alphabet.txt')==pin['alphabet_sha256'] and sha(folder/'parent.input.bin')==pin['input_sha256']
    assert closure['all_letter_identities']==1666 and sum(i['nonzero_ordered_word_coefficients'] for i in closure['checks'].values())==9094 and closure['lowering_sha256']==sha(HERE/'lower.cpp') and ctl['prime']==pin['prime']==closure['parent_field']==P
    result={'prime':P,'lengths':[8,9],'charge_conjugation':-1,'parent_and_recipes':pin,'controls':ctl,'letter_closure':closure,'descendants':{}}
    for op,(dig,n) in BINARY.items():
        plan=read(HERE/(op+'_plan.json'));assert (plan['source_length'],plan['target_length'],plan['parity'],plan['prime'])==(8,9,-1,P);q=folder/('parent.'+op+'.bin');assert sha(q)==dig
        counts={};blocks=[]
        for line in (folder/(op+'.counts')).read_text().splitlines():
            *p,m=map(int,line.split());assert len(p)==9 and tuple(p) not in counts;counts[tuple(p)]=m
        assert set(counts)=={desc(s) for s in plan['targets']} and sum(counts.values())==n
        for s in plan['targets']:
            f=folder/op/filename(s);blocks.append(dict(target=s['index'],terms=counts[desc(s)],bytes=f.stat().st_size,sha256=sha(f)))
        assert q.stat().st_size==sum(b['bytes'] for b in blocks)
        result['descendants'][op]=dict(mu=plan['mu'],terms=n,binary_sha256=dig,bytes=q.stat().st_size,plan_sha256=sha(HERE/(op+'_plan.json')),source_dimension=plan['source_dimension'],target_dimension=plan['target_dimension'],source_profiles=len(plan['sources']),target_profiles=len(plan['targets']),profile_blocks=blocks)
    result['scope']='Modular closed classical descendants. Incoming tests are a separate calculation; neither vector is identified as a higher supercharge obstruction.';return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('output',type=Path);x=p.parse_args();x.output.write_text(json.dumps(receipt(x.folder),indent=2)+'\n')
