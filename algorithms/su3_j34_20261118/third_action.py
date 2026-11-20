"""Direct Q3 on the published even SU(3) BMN source, before either lower correction."""
import argparse,hashlib,json,os,sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import combinations,permutations
from math import lcm
from pathlib import Path
from second_action import read_source,save_rows,cartan,SOURCE_SHA,HERE
from verify import add_scaled
sys.path.insert(0,str(HERE/'third_action_rules'))
from qn_sp_bmn_exact import Qn_bmn_multi,theta_rule
from q3d_sp_exact import Q3d_multi
from fortuity.full import words as fw
import su_trace_adapter_v1 as alg
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def hits(profile):
    eligible=active=0
    for chosen in combinations(range(len(profile)),4):
        if sum(len(profile[i]) for i in chosen)<9:continue
        eligible+=1
        for path in permutations(chosen):
            if path>path[::-1]:continue
            if theta_rule(tuple(profile[i] for i in path)):active+=1
    return eligible,active
@lru_cache(24000)
def fliptrace(w):return fw.canon_word(tuple((b,a,S) for a,b,S in w))
def swap(w):
    z=[];sg=1
    for tr in w:
        n,t=fliptrace(tr);sg*=n;z.append(t)
    n,t=fw.canon_multi(tuple(z));assert sg*n in (1,-1);return sg*n,t
def audit(image,sign):
    assert all(image.get(sw,0)==sign*sg*val for w,val in image.items() for sg,sw in [swap(w)])
def worker(args):
    i,w,c,verify=args;image=Qn_bmn_multi(w,3,3,'SU');audit(image,-1)
    if verify:
        check=Q3d_multi(w,3,'SU',both_paths=True);assert check==image
    return i,{t:v*c for t,v in image.items()},len(image),verify
def run(source,output,workers):
    assert workers>=1 and __debug__;output.mkdir(parents=True,exist_ok=True);s=read_source(source);counts=Counter();profiles=Counter();active=[];dead=[]
    for i,(w,c) in enumerate(sorted(s.items())):
        p=tuple(sorted(L[2] for tr in w for L in tr));profiles[p]+=1;quads,paths=hits(p);counts[quads,paths]+=1
        if paths:active.append((i,w,c))
        else:dead.append((i,w,c))
    assert len(s)==6270 and len(active)+len(dead)==6270
    for i in (0,len(dead)//2,len(dead)-1):
        _,w,_=dead[i];assert not Q3d_multi(w,3,'SU',both_paths=True)
    checks={active[k][0] for k in range(0,len(active),max(1,len(active)//8))}|{active[-1][0]};image={};controls={};nonzero=[];input_terms=0
    with ProcessPoolExecutor(workers) as pool:
        for i,d,n,ck in pool.map(worker,((i,w,c,i in checks) for i,w,c in active),chunksize=2):
            if n:nonzero.append(i)
            if ck:controls[str(i)]=n
            input_terms+=n;add_scaled(image,d,1)
    assert set(controls)==set(map(str,checks));audit(image,-1);assert alg.conjugate(image)==image and all(sum(map(len,t))==9 and fw.sector_of(t)==(5,5,16,10,8) for t in image)
    sha=save_rows(output/'Q3.jsonl.gz',([str(image[w]),w] for w in sorted(image)));d=cartan(image);(output/'Q3_Cartan.json').write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
    result=dict(source_sha256=SOURCE_SHA,source_words=len(s),source_letter_profiles=len(profiles),indexing='zero-based in sorted coloured trace-word order',eligibility=[dict(eligible_quadruples=k[0],active_paths=k[1],source_words=n) for k,n in sorted(counts.items())],surviving_theta_filter=len(active),annihilated_before_colour=len(dead),nonzero_individual_images=len(nonzero),distinct_terms_before_weighted_sum=input_terms,independent_both_orientation_checks=controls,eligible_zero_controls=3,total_path_hits=sum(n*k[1] for k,n in counts.items()),terms=len(image),common_denominator=str(lcm(*(w.denominator for w in image.values()))),length=9,doubled_charges=[5,5,16,10,8],charge_conjugation=1,derivative_exchange=-1,formal_terms_after_diagonal_ordering=d['formal_diagonal_trace_terms'],generic_Cartan_coefficients=len(d['terms']),output_uncompressed_sha256=sha,rules={f.name:digest(f) for f in sorted((HERE/'third_action_rules').glob('*.py'))},scope='Direct Q3 O0 on this representative. The third obstruction requires the full lower descent.')
    (output/'direct.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,default=HERE/'O34_C_even.json.gz');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=1);args=ap.parse_args();run(**vars(args))
