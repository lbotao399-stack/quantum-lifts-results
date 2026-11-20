"""Apply the exact antisymmetric projector to the 318-equation correction.

For an orbit with exactly one old zero residual, the antisymmetric residual is
nonzero at both members. For two nonzero old residuals neither cancellation nor
noncancellation is assumed. Evaluate a representative of every undecided orbit.
"""
import argparse
from fractions import Fraction
from pathlib import Path
import derivative_exchange as e
v,r,s=e.v,e.r,e.s
def run(files,receipt,original,exchange,columns,exterior,work,out,export=None):
    e.control();proof=s.read(exchange);assert proof['rhs']==dict(uncompressed_sha256=s.INPUT,formal_words=517468,negative_exchange_mismatches=0)
    _,ds,_=v.read_inputs(files,receipt);x,ok,bad=v.base315();singles=r.add_rows(*(ds[k] for k in v.SINGLES));x=r.add_rows(x,singles);old=s.read(original)['variants']['three'];assert old['ordered_vector_sha256']==r.vector_sha(x) and not old['lost'];ok|=set(old['gained']);bad-=set(old['gained'])
    raw=e.exchange_vector(x);assert e.exchange_vector(raw)==x;projected=r.add_rows(x,{w:-c for w,c in raw.items()});projected={w:c/2 for w,c in projected.items()};assert e.exchange_vector(projected)=={w:-c for w,c in projected.items()}
    md=r.add_rows(projected,{w:-c for w,c in x.items()});touch=v.changed(md);universe=ok|bad;undecided=set();good=set();fail=set();reason={}
    for t in universe:
        u=e.flip(t)
        if t not in touch:
            (good if t in ok else fail).add(t);reason[t]='unchanged_support'
        elif t in ok and u in ok:good.add(t);reason[t]='both_old_residuals_zero'
        elif (t in ok)!=(u in ok):fail.add(t);reason[t]='one_old_residual_nonzero'
        elif t==u and all(e.TABLE[k]==k for k in bytes.fromhex(t)):good.add(t);reason[t]='exchange_is_identity'
        else:undecided.add(min(t,u))
    work.mkdir(exist_ok=True,parents=True);rows=v.before(projected,undecided);nr,b,rr=rows;tags={('odd',t):rr[t] for t in sorted(undecided)};ev={}
    for seed in (994471,665941):
        r.probe.controls(exterior,seed,work);todo={k:a for k,a in tags.items() if k not in ev};ev.update(r.witnesses(todo,exterior,work,seed))
    evidence={}
    for (name,t),row in tags.items():
        k=(name,t)
        if k in ev:nonzero=True;evv=dict(exterior_nonzero=ev[k])
        elif not row:nonzero=False;evv=dict(formal_terms=0)
        else:
            assert bytes.fromhex(t).count(0)>=4,('unclassified',t);p=r.full({'residual':row},columns,work)['residual'];nonzero=bool(p['nonzero_monomials']);evv=dict(full_polynomial=p)
        evidence[t]=evv
        for u in {t,e.flip(t)}:
            assert u not in (good if nonzero else fail);(fail if nonzero else good).add(u);reason[u]=t
    assert good|fail==universe and not good&fail and {e.flip(t) for t in good}==good
    digest=r.vector_sha(projected);doc=dict(base=old['ordered_vector_sha256'],rhs=s.INPUT,projector='(1 - derivative exchange)/2',ordered_vector_sha256=digest,words=len(projected),delta_words=len(md),affected_contents=len(touch),solved=len(good),failed=len(fail),inhomogeneous=len(good&set(nr)),homogeneous=len(good-set(nr)),unresolved=0,gained=sorted(good-ok),lost=sorted(ok-good),witness_orbits=evidence,reasons=reason);s.write(out,doc)
    if export:
        export.mkdir(exist_ok=True,parents=True);s.write(export/'odd.vector.json.gz',r.word_rows(projected));s.write(export/'odd.census.json',dict(vector_sha=digest,rhs_sha=s.INPUT,exact_zero=sorted(good),mathematically_nonzero=sorted(fail),unresolved=[]))
    print('negative exchange',len(projected),len(good),len(fail),doc['inhomogeneous'],doc['homogeneous'],len(good-ok),len(ok-good),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);h=v.HERE/'joint_70/guarded';a.add_argument('--files',type=Path,default=h);a.add_argument('--receipt',type=Path,default=h/'a8_5665888.json');a.add_argument('--original',type=Path,default=h/'verified.json');a.add_argument('--exchange',type=Path,default=h/'exchange_dimensions.json')
    for n in ('columns','exterior','work','out'):a.add_argument('--'+n,type=Path,required=True)
    a.add_argument('--export',type=Path);assert __debug__;run(**vars(a.parse_args()))
