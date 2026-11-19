"""Test all rational affine combinations of two first corrections, not just sampled parameters.
For each content two verified exterior evaluations with nonzero two-by-two
minor modulo p prove independence of the rational residual polynomials. Emit
both evaluations so the determinant can be recomputed independently.
"""
import argparse
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import full_content_probes as f
s,l,P=f.s,f.l,f.P

def run(rhs,va,vb,pa,pb,exe,out,work,seeds):
    a,b=s.read(pa),s.read(pb); records={};counts={};unresolved=[];seeds=list(dict.fromkeys(seeds));assert seeds and len(a['exact_zero'])+len(a['mathematically_nonzero'])==len(b['exact_zero'])+len(b['mathematically_nonzero'])==500 and not a['unresolved'] and not b['unresolved'] and a['rhs_sha']==b['rhs_sha']
    # First test for an identically constant residual using the exact word map;
    # equality of sampled residues modulo a prime could not replace this check.
    delta=l.load_vector(vb);l.p.axpy(delta,l.load_vector(va),-1);qc={}
    for w,v in delta.items():
        for u,k in f.cen.r.Q_multi(s.expanded(w)).items():l.p.axpy(qc,{s.packed(u):k},v)
    movable={f.cen.r.content(w).hex() for w in qc}
    with TemporaryDirectory(dir=work) as name:
        root=Path(name);pair=[];reference_checks={}
        for idx,(path,z) in enumerate(((va,a),(vb,b))):
            d=root/str(idx);d.mkdir(); files={}
            def emit(g,w,x):
                if g not in files:files[g]=(d/g.hex()).open('w')
                x=Fraction(x);assert x.denominator%P;files[g].write(f'{x.numerator*pow(x.denominator,-1,P)%P} {len(w)} '+ ' '.join(s.hexes(w))+'\n')
            try:c=f.cen.inventory(rhs,path,emit)
            finally:
                for fp in files.values():fp.close()
            assert c['source_sha256']==z['vector_sha'] and c['rhs_sha256']==z['rhs_sha'];rows=d/'rows'
            with rows.open('wb') as to:
                import shutil
                for g in sorted(files):
                    to.write(f"{g.hex()} {sum(l.f.ODD[k] for k in g)} {c['contents'][g.hex()][2]}\n".encode())
                    with (d/g.hex()).open('rb') as src:shutil.copyfileobj(src,to)
            side={}
            for seed in seeds:
                if idx==0:reference_checks[seed]=f.controls(exe,seed,d)
                vv=f.read_data(rows,exe,seed)
                for t,r in vv.items():
                    for mask,x in r.items():
                        if x:side[(t,seed,mask)]=x
            assert set(t for t,seed,mask in side)==set(z['mathematically_nonzero']);pair.append(side)
        A,B=pair;allkeys=sorted(set(A)|set(B)); by={}
        for t,u,k in allkeys:by.setdefault(t,[]).append((u,k,A.get((t,u,k),0),B.get((t,u,k),0)))
        az=set(a['exact_zero']);bz=set(b['exact_zero']);allvars=set(c['contents']);assert az|set(a['mathematically_nonzero'])==bz|set(b['mathematically_nonzero'])==allvars
        for t in sorted(allvars):
            if t in az and t in bz:kind='all parameters';w=[]
            elif t in az:kind='first endpoint only';w=[]
            elif t in bz:kind='second endpoint only';w=[]
            else:
                r=by[t];x=r[0];y=next((q for q in r[1:] if (x[2]*q[3]-q[2]*x[3])%P),None)
                if t not in movable:assert x[2]==x[3] and x[2] and y is None;kind='nonzero constant as formal words';w=[x]
                elif y:kind='no parameter';w=[x,y,(x[2]*y[3]-y[2]*x[3])%P]
                else:
                    kind='not decided by the two-column test';w=[x];unresolved.append(t)
            records[t]=[kind,w];counts[kind]=counts.get(kind,0)+1
    n=counts.get('all parameters',0);outrec=dict(rhs=a['rhs_sha'],vectors=[a['vector_sha'],b['vector_sha']],prime=P,seeds=seeds,reference_checks=reference_checks,entries=records,nonzero_formal_difference_words=len(qc),different_formal_contents=sorted(movable),scope_counts=counts,unresolved=unresolved,endpoint_success=[len(az),len(bz)],other_parameter_lower_bound=n,other_parameter_upper_bound=n+len(unresolved),global_max_if_endpoint_certified=(max(len(az),len(bz)) if max(len(az),len(bz))>=n+len(unresolved) else None),row_order=['seed','exterior_mask','first_residual','second_residual'],interpretation='A nonzero minor excludes every rational affine combination of the two fixed corrections for this content. A zero minor modulo one prime gives no rank upper bound over Q.');s.write(out,outrec);print(counts,'generic bound',n,n+len(unresolved),flush=True);return outrec
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');p.add_argument('--vectors',type=Path,nargs=2,required=True);p.add_argument('--certificates',type=Path,nargs=2,required=True);p.add_argument('--exe',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seeds',type=int,nargs='+',default=[817234,419071]);a=p.parse_args();assert __debug__;a.work.mkdir(exist_ok=True,parents=True);run(a.rhs,*a.vectors,*a.certificates,a.exe,a.out,a.work,a.seeds)
