"""Complete the last C4 boundary allowing closed corrections one level above."""
from pathlib import Path
import pickle,json,time,resource
import certify_c4 as C,certify_nonzero as N,exact_product as E
ROOT=Path(__file__).resolve().parent
C.CHECK_PATH=ROOT/'final_boundary_checks.json'
START=time.time()

def save(name,p):
    with (ROOT/name).open('wb') as f:pickle.dump(p,f)

def q(n,p):return E.Q0(p) if n==0 else E.T.qn(n,p)

def main():
    d=C.load('c4_hidden_product_exact.pkl')
    full=C.load('full_products_partial.pkl')
    P,L=full['P'],full['L'];bases=C.load('completion_bases_mod.pkl')
    prim=[E.sub(d['C2_squared_primitives'][0],E.sc(d['C0_C4_primitives'][0],3)),
          E.add(E.sub(d['C2_squared_primitives'][1],E.sc(d['C0_C4_primitives'][1],3)),d['boundary_Q1']),
          d['boundary_Q0']]
    p3path=ROOT/'last_boundary_primitive3.pkl'
    if p3path.exists():prim.append(C.load(p3path.name))
    else:
        target=E.sub(E.sub(P[3],E.sc(L[3],3)),E.add(*(q(3-i,p) for i,p in enumerate(prim))))
        source=[C.poly(p) for p in bases[(1,3,6)]]
        co,_=C.rational_solve([E.Q0(p) for p in source],target)
        prim.append(E.add(*(E.sc(p,c) for p,c in zip(source,co))))
        save(p3path.name,prim[-1])
    save('last_boundary_primitives_partial.pkl',prim)
    pieces=[]
    for i,p in enumerate(prim):
        name='last_boundary_Q%d_primitive%d.pkl'%(4-i,i)
        C.say('START '+name+' input terms '+str(len(p)))
        if (ROOT/name).exists():value=C.load(name)
        else:
            value=q(4-i,p);save(name,value)
        pieces.append(value)
        C.say('DONE '+name+' output terms '+str(len(value))+' seconds '+str(time.time()-START))
    target=E.sub(E.sub(P[4],E.sc(L[4],3)),E.add(*pieces))
    save('last_boundary_target_exact.pkl',target)
    C.check('last boundary target Q0 closed',E.Q0(target))
    source=[C.poly(p) for p in bases[(0,3,6)]]
    space=C.load('last_boundary_space_mod.pkl')
    closed=[C.poly(p) for p in space['closed']]
    for p in closed:C.check('last-stage closed correction',E.Q0(p))
    columns=[E.Q0(p) for p in source]+[E.Q1(p) for p in closed]
    ix=N.indices(columns)
    co,_=C.rational_solve([columns[i] for i in ix],target)
    u=E.add(*(E.sc(source[i],c) for i,c in zip(ix,co) if i<len(source)))
    z=E.add(*(E.sc(closed[i-len(source)],c) for i,c in zip(ix,co) if i>=len(source)))
    C.check('full last boundary equation',E.sub(target,E.add(E.Q0(u),E.Q1(z))))
    C.check('closed correction preserves previous boundary equation',E.Q0(z))
    prim[3]=E.add(prim[3],z);prim.append(u)
    full['boundary']=prim
    save('full_hidden_product_exact.pkl',full)
    summary=dict(passed=True,seconds=time.time()-START,maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 runtime_scope='This assembler process only; saved Qn images were computed in earlier workers.',
                 computation_logs=['complete.log','finish.log','precompute.log'],
                 boundary_terms=list(map(len,prim)),last_closed_correction_terms=len(z),
                 relation='C2 star C2 - 3 C0 star C4_source = Qpert boundary',
                 article_normalization='C4_article = -C4_source gives C2 star C2 + 3 C0 star C4_article = Qpert boundary')
    (ROOT/'full_completion_summary.json').write_text(json.dumps(summary,indent=2))
    C.say(json.dumps(summary))

if __name__=='__main__':main()
