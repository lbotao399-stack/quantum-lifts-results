"""Derive the small parent-extraction input manifest from the preceding committed receipts."""
import hashlib,json
from pathlib import Path
from zipfile import ZipFile

def make(repo):
    here=repo/'data';load=lambda p:json.loads(p.read_bytes());h=lambda b:hashlib.sha256(b).hexdigest();j=load(here/'N4_SU3_j30_joint_rank_20261119.json');r=load(here/'N4_SU3_j30_incoming_rank_20261119.json');plan=load(repo/'algorithms/su3_j30_classical_ranks/joint_plan.json')
    with ZipFile(here/'N4_SU3_j30_target_frames_p1000033.zip') as z:
        frames=[[h(z.read(f'e2_t{i:02d}/{name}.json')) for name in ('manifest','result')] for i in range(37)]
    return dict(joint=j['stages'],incoming=r['stages'],frames=frames,recipes=[s['recipe_sha256'] for s in plan['sources']],original_plan_sha256=r['input_plan_sha256'],matrices_sha256=r['matrix_receipt_sha256'])

if __name__=='__main__': print(json.dumps(make(Path(__file__).resolve().parents[2]),separators=(',',':')))
