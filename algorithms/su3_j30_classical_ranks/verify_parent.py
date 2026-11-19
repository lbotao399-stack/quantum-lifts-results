"""Offline receipt audit of the length-nine parent. No Slurm, pickle, or native code needed.

Checks recorded minor determinants, dimensions, hashes, and rank jumps. Full
matrix regeneration uses parent.py, closure_frames.py, and outgoing_frames.py.
"""
import argparse
import ast
import hashlib
import json
import struct
from pathlib import Path
from zipfile import ZipFile

from parent_pins import make as make_pins

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data'
BUNDLE = DATA / 'SU3_j30_parent_p1000033'
PRIME = 1000033
HEADER_HASH = '5759f1e92afbf02b7255c3ab249c039330d90fb1387f71c122ed677579eb29c3'

def sha(data): return hashlib.sha256(data).hexdigest()
def read(file): return json.loads(file.read_bytes())

def check_npy(file, digest, dimension, expected_nonzero):
    raw = file.read_bytes()
    assert sha(raw) == digest and raw[:8] == b'\x93NUMPY\x01\x00'
    offset = 10 + int.from_bytes(raw[8:10], 'little')
    head = ast.literal_eval(raw[10:offset].decode('ascii').strip())
    assert head == {'descr':'<u4','fortran_order':False,'shape':(dimension,)}
    assert len(raw) == offset + 4 * dimension
    values = [v for v, in struct.iter_unpack('<I', raw[offset:])]
    assert max(values) < PRIME and sum(v != 0 for v in values) == expected_nonzero

def witness(record, dimension, width):
    assert record['prime'] == PRIME and record['saturated'] is True
    w = record['witness']; rows = w['rows']; cols = w['columns']
    assert record['dimension'] == record['rank'] == w['rank'] == dimension == len(rows) == len(cols)
    assert record['width'] == width and len(set(cols)) == dimension and all(0 <= c < width for c in cols)
    assert len({json.dumps(row,separators=(',',':')) for row in rows}) == dimension
    n = 0; det = 1
    for order, factor in w['panel_orders_and_determinants']:
        assert 0 < order <= dimension and 0 < factor < PRIME
        n += order; det = det * factor % PRIME
    assert n == dimension and det == w['determinant_in_selected_order'] != 0

def run():
    assert __debug__ and sha((BUNDLE/'receipt_hashes.json').read_bytes()) == HEADER_HASH
    meta = read(BUNDLE/'receipt_hashes.json'); parent = read(BUNDLE/'result.json')
    assert sha((BUNDLE/'result.json').read_bytes()) == meta['parent']['result.json']
    inc = read(DATA/'N4_SU3_j30_incoming_rank_20261119.json'); joint = read(DATA/'N4_SU3_j30_joint_rank_20261119.json')
    pins = make_pins(ROOT); plan = read(Path(__file__).with_name('joint_plan.json'))
    assert parent['prime'] == inc['prime'] == joint['prime'] == plan['prime'] == PRIME
    assert parent['recipe_vector_sha256'] == meta['parent']['parent.npy'] and parent['evaluation_sha256'] == meta['parent']['evaluated.npy']
    check_npy(BUNDLE/'parent.npy',parent['recipe_vector_sha256'],parent['parent_dimension'],parent['parent_nonzero'])
    check_npy(BUNDLE/'evaluated.npy',parent['evaluation_sha256'],parent['evaluations'],parent['evaluation_nonzero'])
    out = parent['outgoing_and_Cartan']; stages = {v['stage']:v for v in joint['stages']}; assert len(stages) == len(out) == 51
    assert [(s['key'],s['coordinates']) for s in out] == [(s['stage'],s['coordinates']) for s in joint['stages']]
    assert all(s['chunks'] > 0 and s['nonzero'] == 0 for s in out)
    preds = {v['index']:v for v in read(DATA/'N4_SU3_j30_source_order_20261119.json')['profiles'] if v['equation']==2}
    cart = {'c%02d'%s['index']:s for s in plan['sources'] if s['cartan_dimension']}
    dims = {'q%02d'%s['index']:s['dimension'] for s in plan['targets']}
    existing = {k for k in dims if stages[k]['rank_increment']==dims[k]}
    assert len(existing)==9 and len(preds)==57 and len(cart)==32
    for name,kind in [('cartan_predecessor_minors.zip','cartan_predecessor'),('outgoing_minors.zip','outgoing')]:
        hashes=meta[kind]; assert sha((BUNDLE/name).read_bytes())==meta['archives'][name]
        expected=(set(cart)|{'p%02d'%i for i in preds}) if kind=='cartan_predecessor' else set(dims)-existing
        assert set(hashes)==expected
        with ZipFile(BUNDLE/name) as z:
            assert z.testzip() is None and sorted(z.namelist())==[k+'.json' for k in sorted(expected)]
            for k in sorted(expected):
                raw=z.read(k+'.json'); assert sha(raw)==hashes[k]; v=json.loads(raw); assert v['key']==k
                i=int(k[1:]); w=v['witness']
                if k[0]=='p':
                    s=preds[i]; assert v['recipe_sha256']==s['seeds_sha256'];witness(v,s['dimension'],v['width'])
                    assert set(w['rows'])==set(range(s['dimension']))
                    if k=='p08': assert raw==(BUNDLE/'local_p08.json').read_bytes() and sha(raw)==meta['local_p08_sha256'] and w['determinant_in_selected_order']==122908
                elif k[0]=='c':
                    s=cart[k];assert v['recipe_sha256']==pins['recipes'][i] and v['stage_sha256']==stages[k]['stage_sha256']
                    witness(v,s['cartan_dimension'],stages[k]['coordinates']);assert all(0<=row<s['dimension'] for row in w['rows'])
                else:
                    assert v['plan']==joint['plan'] and v['stage_sha256']==stages[k]['stage_sha256']
                    witness(v,dims[k],stages[k]['coordinates']); m=v['image_rank']; assert v['increment']==stages[k]['rank_increment']<=m<=dims[k]
                    sources=[s['index'] for s in plan['sources'] if i in s['targets']];assert sources==v['sources']
                    assert sum(plan['sources'][s]['dimension'] for s in sources)==v['image_rows']
                    assert all(s in sources and 0<=j<plan['sources'][s]['dimension'] for s,j in w['rows'][:m]);assert all(label=='trace' and sum(map(len,word))==10 for label,word in w['rows'][m:])
    checks=parent['incoming_checks']; assert parent['nonboundary_mod_prime'] and 0<len(checks)<=len(inc['stages'])
    for i,(s,z) in enumerate(zip(checks,inc['stages'])):
        jump=(i==len(checks)-1);assert s['target']==z['target'] and s['old_increment']==z['rank_increment'] and s['inconsistent']==jump and s['augmented_increment']==s['old_increment']+int(jump)
    n=sum(s['dimension'] for s in plan['sources']);rank=sum(s['rank_increment'] for s in inc['stages']);both=sum(s['rank_increment'] for s in joint['stages'])
    assert n==parent['parent_dimension']==63983 and rank==inc['incoming_rank_over_prime']==parent['incoming_rank']==32144 and both==joint['joint_measured_rank']==parent['joint_rank']==31838
    assert parent['free']==n-both and parent['augmented_rank']==rank+1 and n-both-rank==1
    return {'prime':PRIME,'cochains':n,'incoming_rank':rank,'augmented_rank':rank+1,'outgoing_and_Cartan_rank':both,'outgoing_injective_profiles':len(dims),'outgoing_zero_coordinates':sum(s['coordinates'] for s in out if s['key'][0]=='q'),'Cartan_injective_profiles':len(cart),'Cartan_zero_coordinates':sum(s['coordinates'] for s in out if s['key'][0]=='c'),'complete_predecessor_profiles':len(preds),'Cartan_invisible_dimension_mod_prime':n-both-rank,'Cartan_invisible_dimension_rational_upper_bound':1,'rational_second_obstruction_decided':False}

if __name__=='__main__': print(json.dumps(run(),indent=2))
