"""Solve selected first-correction equations together on a specified rational base.

Parents incident on several targets share one coefficient. Protect either every
previous identity (modulo explicitly retested ones), or use *all* parents of the
selected targets. In the latter case a zero on the selected equations carries
no assertion on the other 500-content equations. A rational dual witness is
normalized against the original signed residual, not its cleared denominator.
"""
import argparse
from fractions import Fraction
from math import lcm
from pathlib import Path
import inhomogeneous_block as single
j,s=single.j,single.s

def prepare(targets,guards,rhs,vector,census,out):
    ts=sorted(set(targets)|set(guards)); assert targets and len(ts)==len(targets)+len(guards)
    old=s.read(census); ok=set(old['exact_zero']); bad=set(old['mathematically_nonzero']);F,B=single.front.front.incidence()
    assert not old['unresolved'] and not ok&bad and ok|bad=={k.hex() for k in B} and set(targets)<=bad and set(guards)<=ok
    digest,counts,res=single.front.residuals(rhs,vector,ts);assert digest==old['vector_sha'] and old['rhs_sha']==s.INPUT
    gs=sorted({g for t in ts for g in B[bytes.fromhex(t)]}); H={g.hex():{k.hex() for k in F[g]} for g in gs};keep=sorted(g for g,t in H.items() if not (t & (ok-set(guards))))
    dims={g:j.E.parity_frame_dims(bytes.fromhex(g))[0] for g in list(H)+ts}; assert all(dims.values()); assert all(0<len(bytes.fromhex(k))-bytes.fromhex(k).count(0)<=7 for k in dims)
    z=dict(prime=j.P,base=digest,rhs=s.INPUT,targets=ts,failed_requested=sorted(targets),solved_guards=sorted(guards),passes_at_start=len(ok),groups={g:sorted(set(ts)&v) for g,v in H.items()},internal=keep,dimensions=dims,scope='Internal parents avoid all other previous successes. All includes every parent of each selected content, including retested successes.')
    out.mkdir(parents=True,exist_ok=True);den=lcm(*(v.denominator for row in res.values() for v in row.values()));z['denominator']=str(den)
    for t in ts:s.write(out/(t+'.rhs.json.gz'),[[s.hexes(w),str(int(-v*den))] for w,v in sorted(res[t].items())]);
    z['formal_rhs_words']={t:counts[t] for t in ts};s.write(out/'graph.json',z);return z

def solve(out,frames,scope):
    z,A,y,pool,keys=j.assemble(out,scope,frames);history=[]
    def progress(v):history.append(v);s.write(out/(scope+'.progress.json'),history);print(scope,v,flush=True)
    rec,co=single.certify(A,y,progress);D=int(z['denominator']);gs=sorted(set(g for g,w in pool));rec.update(base=z['base'],rhs=z['rhs'],scope=scope,targets=z['targets'],failed_requested=z['failed_requested'],solved_guards=z['solved_guards'],parents=gs,parent_coordinates=len(pool),target_coordinates=len(keys),matrix_entries=int(A.nnz)); dest=out/(scope+'.result.json')
    if rec['kind']=='exact_separator':
        # Keys include the target content, so distinct summands cannot alias.
        dd=lcm(*(v.denominator for v in co.values()));c=[int(co.get(i,0)*dd) for i in range(len(y))];assert all(v==0 for v in j.multiply(A,c));p=Fraction(sum(a*b for a,b in zip(c,y)),D*dd);assert p==Fraction(rec['integer_rhs_pairing'])/D and p;s.write(out/(scope+'.functional.json.gz'),[[t,str(k),str(co[i])] for i,(t,k) in enumerate(keys) if i in co]);rec.update(nonzero_functional=len(co),pairing=str(p));
    elif rec['kind']=='exact_solution':
        delta={};terms=[]
        for i,a in sorted(co.items()):
            g,w=pool[i];v=a/D;terms.append([g,w,str(v)]);j.l.p.axpy(delta,j.E.projected(s.unpack_hex(w)),v)
        s.write(out/(scope+'.frame_terms.json.gz'),terms);j.l.save_vector(out/(scope+'.delta.json.gz'),delta);rec.update(frame_terms=len(terms),delta_words=len(delta))
    s.write(dest,rec);return rec

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='cmd',required=True);a=subs.add_parser('prepare');a.add_argument('--targets',nargs='+',required=True);a.add_argument('--guards',nargs='*',default=[])
    for n in ('rhs','vector','census','out'):a.add_argument('--'+n,type=Path,required=True)
    b=subs.add_parser('solve');b.add_argument('--out',type=Path,required=True);b.add_argument('--frames',type=Path,required=True);b.add_argument('--scope',choices=['internal','all'],required=True);x=vars(p.parse_args());assert __debug__;print({'prepare':prepare,'solve':solve}[x.pop('cmd')](**x).get('kind','prepared'),flush=True)
