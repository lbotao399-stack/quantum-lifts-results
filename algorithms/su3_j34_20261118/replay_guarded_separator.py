"""Rebuild the rational functional on all restricted parent basis elements.

Use the upstream Q_multi and parity_project on the original three-colour words.
For each output word, evaluate only the monomial coefficients in the supplied
functional. Neither the elimination matrix nor the solver's pivot lists is read.
"""
import argparse,hashlib
from collections import Counter,defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path
from su_trace_adapter_v1 import parity_project
import replay_guarded as v
r,s=v.r,v.s

def run(files,receipt,exe,work,out):
    rec=s.read(receipt)['receipt'];q=files/'pair';cert=q/'internal.result.json';func=q/'internal.functional.json.gz';fr=q/'parent_frames';job=rec['results']['pair/internal']
    for p in (cert,func):assert job['files'][str(p.relative_to(files))]==[p.stat().st_size,r.sha(p)]
    z=s.read(cert);assert z['kind']=='exact_separator' and rec['base']==z['base']==v.BASE_DIGEST and rec['rhs']==z['rhs']==s.INPUT
    W=defaultdict(dict)
    for t,k,c in s.read(func):
        k=int(k);assert t in z['targets'] and k not in W[t] and Fraction(c);W[t][k]=Fraction(c)
    D=lcm(*(c.denominator for row in W.values() for c in row.values()));wi={t:{k:int(a*D) for k,a in row.items()} for t,row in W.items()};ts=set(z['targets']);x,ok,bad=v.base315()
    import parent_frontier as pf
    F,B=pf.incidence();guards=set(z['solved_guards']);all_g={g for t in ts for g in B[bytes.fromhex(t)]};allowed=sorted(g.hex() for g in all_g if {k.hex() for k in F[g]}.isdisjoint(ok-guards))
    assert allowed==z['parents'] and set(z['failed_requested'])<=bad and guards<=ok and len(ts)==4
    _,_,rr=v.before(x,ts);wanted={t:set(rr[t]) for t in wi};sources={};records={};work.mkdir(exist_ok=True,parents=True)
    for g in allowed:
        path=fr/(g+'.json.gz');base=s.read(path);assert r.sha(path)==rec['frames'][g]['sha256']
        det=r.E.check(base,exe,work);rows=[];assert det==rec['frames'][g]['det']
        for raw in base['words']:
            w=s.expanded(s.unpack_hex(raw));image=Counter()
            for p,c in parity_project({w:1},1).items():
                for u,k in s.Q_multi(p).items():
                    wz=s.packed(u);t=r.ch(wz)
                    if t in wi:image[wz]+=int(c*k)
            clean={w:k for w,k in image.items() if k};rows.append(clean)
            for w in clean:wanted[r.ch(w)].add(w)
        records[g]=dict(dimension=len(rows),independent_minor_mod_prime=det);sources[g]=rows
        print('parent frame',g,len(rows),flush=True)
    values={};word_count={}
    for t,ws in sorted(wanted.items()):
        words=sorted(ws);word_count[t]=len(words)
        for i in range(0,len(words),64):
            chunk=words[i:i+64]
            for w,p in zip(chunk,r.E.poly_iter(chunk,exe,work),strict=True):values[w]=sum(wi[t].get(k,0)*c for k,c in p)
    for g,rows in sources.items():
        tests=[sum(k*values[w] for w,k in row.items()) for row in rows];assert not any(tests);records[g]['annihilated_columns']=len(tests)
    pair=-sum(k*values[w] for t,row in rr.items() if t in wi for w,k in row.items())/D
    assert isinstance(pair,Fraction) and pair and pair==Fraction(z['pairing']);P=z['prime'];res=pair.numerator*pow(pair.denominator,-1,P)%P;assert sum(d['dimension'] for d in records.values())==z['parent_coordinates']
    s.write(out,dict(job=rec['job'],source=z['base'],rhs=s.INPUT,targets=sorted(ts),parent_frames=records,annihilated_parent_coordinates=z['parent_coordinates'],target_coordinates=z['target_coordinates'],denominator=str(D),nonzero_functional_coefficients=sum(map(len,W.values())),pairing_sign='negative of Q0 O1 + Q1 O0',pairing=str(pair),prime=P,residue=res,trace_words_per_supported_content=word_count))
    print('annihilator',res,z['parent_coordinates'],flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);h=v.HERE/'joint_70/guarded';a.add_argument('--files',type=Path,default=h);a.add_argument('--receipt',type=Path,default=h/'a8_5665888.json')
    for n in ('exe','work','out'):a.add_argument('--'+n,type=Path,required=True)
    assert __debug__;run(**vars(a.parse_args()))
