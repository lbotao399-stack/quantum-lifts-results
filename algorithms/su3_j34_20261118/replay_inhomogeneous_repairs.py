"""Reconstruct both A8 corrections and classify their sum from the original RHS.

Uses a packed cyclic/graded Leibniz differential independently of the producer's
upstream Q_multi. Support is computed from each original letter's splits. Zero
on a random assignment is never counted as equality. Only a complete polynomial
or a zero inherited with disjoint support certifies a solved equation.
"""
import argparse, gzip, hashlib, json
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import formal_trace_blocks as f
import full_content_probes as probe
import even_content_frames as E
s=probe.s
HERE=Path(__file__).resolve().parent
BASE=HERE/'joint_70/support_cut'
TARGETS=('0000000000010105050716','0000000001010106060d0f')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ch(w):return bytes(sorted(b''.join(w))).hex()
def add(x,w,a):
    if a:
        x[w]=x.get(w,0)+a
        if not x[w]:del x[w]
def add_rows(*rows):
    out={}
    for row in rows:
        for w,c in row.items():add(out,w,c)
    return out
def load_words(p):
    data=s.read(p);row={s.unpack_hex(w):Fraction(c) for w,c in data};assert len(row)==len(data) and all(row.values());return row
def word_rows(row):return [[s.hexes(w),str(c)] for w,c in sorted(row.items())]
def vector_sha(x):return hashlib.sha256(json.dumps(word_rows(x),separators=(',',':')).encode()).hexdigest()
def splits(g):
    v=bytes.fromhex(g);return {bytes(sorted(v[:i]+pair+v[i+1:])).hex() for i,k in enumerate(v) for _,pair in f.SPLITS[k]}
def project_terms(path):
    row={};seen=set();terms=s.read(path)
    for g,raw,c in terms:
        w=s.unpack_hex(raw);v=Fraction(c);assert g==ch(w) and v and w not in seen;seen.add(w)
        for u,k in E.projected(w).items():add(row,u,k*v)
    return row,len(seen)
def action(vec,targets):
    rows={t:{} for t in targets};possible={}
    for w,c in vec.items():
        g=ch(w)
        if g not in possible:possible[g]=bool(splits(g)&targets)
        if not possible[g]:continue
        for u,k in f.q(w).items():
            t=ch(u)
            if t in rows:add(rows[t],u,k*c)
    return rows

def full(rows,exe,work):
    # Export a word only once even if it occurs on both sides or on multiple vectors.
    words=sorted(set().union(*rows.values()));scaled={};acc={};den={}
    for key,row in rows.items():
        den[key]=lcm(*(c.denominator for c in row.values()));scaled[key]={w:int(c*den[key]) for w,c in row.items()};acc[key]=Counter()
    for start in range(0,len(words),64):
        chunk=words[start:start+64]
        for w,stream in zip(chunk,E.poly_iter(chunk,exe,work),strict=True):
            col=list(stream)
            for key,cc in scaled.items():
                if w in cc:
                    n=cc[w];a=acc[key]
                    for j,k in col:a[j]+=n*k
    out={}
    for key,a in acc.items():
        nz=[j for j in a if a[j]];j=min(nz,default=None);out[key]=dict(formal_words=len(rows[key]),collected_monomials=len(a),nonzero_monomials=len(nz),first_key=j,first_coefficient=str(Fraction(a[j],den[key])) if j is not None else None)
    return out

def witnesses(rows,exe,work,seed):
    with TemporaryDirectory(dir=work) as folder:
        path=Path(folder)/'samples';ids={};res={}
        with path.open('w') as st:
            for i,((variant,t),row) in enumerate(rows.items()):
                if not row:continue
                ids[str(i)]=(variant,t);st.write(f'{i} {sum(f.ODD[k] for k in bytes.fromhex(t))} {len(row)}\n');den=lcm(*(v.denominator for v in row.values()));inv=pow(den,-1,probe.P)
                for w,v in row.items():st.write(f'{(int(v*den)*inv)%probe.P} {len(w)} '+ ' '.join(s.hexes(w))+'\n')
        if ids:
            for i,vals in probe.read_data(path,exe,seed).items():
                found=next(((k,v) for k,v in sorted(vals.items()) if v),None)
                if found:res[ids[i]]=[seed,*found]
        return res

def run(files,receipt_file,columns,exterior,work,out):
    work.mkdir(exist_ok=True,parents=True); rec=s.read(receipt_file)['receipt'];assert rec['job']=='5665806' and set(rec['targets'])==set(TARGETS)
    rhs=HERE/'negative_Q1.jsonl.gz';vecfile=BASE/'preserve_296.vector.json.gz';cfile=BASE/'preserve_296.census.json';assert rec['inputs']=={p.name:sha(p) for p in (rhs,vecfile,cfile)};old=s.read(cfile);base=load_words(vecfile);ok=set(old['exact_zero']);bad=set(old['mathematically_nonzero']);assert len(ok)==313 and len(bad)==187 and not old['unresolved'] and not ok&bad and vector_sha(base)==old['vector_sha'];docs={};vec={};touch={}
    for t in TARGETS:
        d=files/t;data=rec['targets'][t];res=s.read(d/'result.json');assert res==data['result'] and res['kind']=='exact_solution' and res['base']==old['vector_sha'] and res['target']==t
        for name in ('result.json','delta.json.gz','frame_terms.json.gz'):
            path=d/name;assert [path.stat().st_size,sha(path)]==data['files'][name]
        delta,nt=project_terms(d/'frame_terms.json.gz');assert delta==load_words(d/'delta.json.gz') and nt==res['frame_terms'] and len(delta)==res['delta_words'];x=add_rows(base,delta);assert len(x)==res['vector_words']
        # gzip encodings depend on the Python version; compare ordered JSON words if a raw file is supplied.
        gp=sorted({ch(w) for w in delta});assert gp==res['changed_parents'] and set(gp)<=set(res['allowed']);touch[t]=set.union(*(splits(g) for g in gp));assert touch[t]==set(res['affected_contents']) and touch[t]<=bad
        docs[t]=res;vec[t]=delta
    a,b=TARGETS;assert a not in touch[b] and b not in touch[a];union=touch[a]|touch[b];assert len(union)==22
    bb={t:{} for t in union};nr=Counter();hd=hashlib.sha256();n=0
    for w,c in s.stream_rows(rhs,hd):
        t=ch(w);nr[t]+=1;n+=1
        if t in bb:assert w not in bb[t];bb[t][w]=Fraction(c)
    assert n==517468 and len(nr)==174 and hd.hexdigest()==s.INPUT==old['rhs_sha']; im=action(base,union);rr={t:add_rows(im[t],{w:-c for w,c in bb[t].items()}) for t in union}; acts={t:action(vec[t],touch[t]) for t in TARGETS};both=add_rows(base,vec[a],vec[b]); assert len(both)<=len(base)+len(vec[a])+len(vec[b]);print('reconstructed words:',len(base),*(len(add_rows(base,vec[t])) for t in TARGETS),len(both),flush=True)
    pending={};ans={};detail={};vs=(*TARGETS,'combined');whole={a:add_rows(base,vec[a]),b:add_rows(base,vec[b]),'combined':both};touch['combined']=union
    for t in TARGETS:
        final=add_rows(rr[t],acts[t][t]);check=full({'rhs':bb[t],'before':rr[t],'after':final},columns,work);stored=docs[t];assert check['rhs']==stored['physical_before']['negative_Q1'] and check['rhs']['nonzero_monomials'] and check['before']==stored['physical_before']['Q0_O1_plus_Q1_O0'] and check['before']['nonzero_monomials'];assert check['after']==stored['physical_after'] and check['after']['nonzero_monomials']==0;detail[t]=check
        for v in (t,'combined'):ans[v,t]={'exact_polynomial':'after','nonzero':False}
    for v in vs:
        for t in sorted(touch[v]):
            if (v,t) not in ans:pending[v,t]=add_rows(rr[t],*(acts[k][t] for k in TARGETS if (v=='combined' or v==k) and t in acts[k]))
    seed=817234;probe.controls(exterior,seed,work);sample=witnesses(pending,exterior,work,seed)
    for key,row in pending.items():
        if key in sample:ans[key]={'exterior_nonzero':sample[key],'nonzero':True}
        elif not row:ans[key]={'formal_terms':0,'nonzero':False}
        else:
            _,t=key;assert bytes.fromhex(t).count(0)>=4,('unclassified; expand at a different diagonal scalar or use more exterior witnesses',key);p=full({'r':row},columns,work)['r'];ans[key]={'exact_polynomial':p,'nonzero':bool(p['nonzero_monomials'])}
    summary={}
    for v in vs:
        affected=touch[v];gains={t for t in affected if not ans[v,t]['nonzero']};success=ok|gains;fail=(bad-affected)|{t for t in affected if ans[v,t]['nonzero']};assert len(success|fail)==500 and success.isdisjoint(fail);summary[v]=dict(vector_words=len(whole[v]),ordered_vector_sha256=vector_sha(whole[v]),affected_contents=sorted(affected),gained=sorted(gains),lost=[],inhomogeneous_passes=len(success&set(nr)),homogeneous_passes=len(success-set(nr)),fixed_vector_zero=len(success),fixed_vector_nonzero=len(fail),unresolved=0,changed_content_evidence={t:ans[v,t] for t in sorted(affected)})
    result=dict(a8_job=rec['job'],rhs_sha=s.INPUT,base_vector_sha=old['vector_sha'],base_census_sha256=sha(cfile),sign='Q0 O1 + Q1 O0',exterior_prime=probe.P,shared_affected_contents=sorted(touch[a]&touch[b]),full_traceless_three_colour=detail,variants=summary);s.write(out,result);print([(k,x['fixed_vector_zero'],x['inhomogeneous_passes'],x['homogeneous_passes']) for k,x in summary.items()],flush=True);return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);h=HERE/'joint_70/inhomogeneous';ap.add_argument('--files',type=Path,default=h);ap.add_argument('--receipt-file',type=Path,default=h/'a8_5665806.json');ap.add_argument('--columns',type=Path,required=True);ap.add_argument('--exterior',type=Path,required=True);ap.add_argument('--work',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);assert __debug__;run(**vars(ap.parse_args()))
