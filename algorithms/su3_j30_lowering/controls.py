"""Compare native operators term by term to the upstream word maps at the parent's prime."""
import argparse, json, random, subprocess, sys
from pathlib import Path
from inputs import P, PLAN, LETTERS, HERE, PACKAGE, RUNTIME, sha, key, unpack, records, write_records, write_alphabet


def run(up,exe,out):
    assert __debug__; out.mkdir(parents=True,exist_ok=True); alpha=out/'alphabet.txt';write_alphabet(alpha)
    code=up/'code'; scripts=up/'scripts/oneloop/su-threshold-classes'
    pin={'build_cand.py':'7edf415f3a37fd312a51ab28924fc05575a35cfd8a6f17e4153e362ac1added5','words.py':'87dc52706e6b3ec4efca6dccfd9b76eb823b8e91142524503046c046f5fc131d','letters.py':'0e68c21fc553197aa737ea673294647a9e9cf57a51a38fff4dd226dccea846be'}
    for f in (scripts/'build_cand.py',code/'fortuity/full/words.py',code/'fortuity/full/letters.py'):assert sha(f)==pin[f.name]
    sys.path[:0]=[str(code),str(scripts)];import build_cand as ref
    from fortuity.full.words import canon_word, canon_multi, word_parity, Q_expr
    ref.P=P; ops={'S20':ref._S20(3),'Tp21':ref._Tp21(3)}
    canon=lambda a:{w:int(v)%P for w,v in a.items() if v%P}
    def at(w):
        sign=1;tr=[]
        for t in w:
            s,l=canon_word(t)
            if not s:return 0,None
            sign*=s;tr.append(l)
        s,l=canon_multi(tr);return s*sign,l
    def baseline(words,op):
        out={}
        for w,c in words.items():
            for fac,f in ops[op]:
                for z,k in f(w).items():
                    s=1
                    if op=='Tp21':s,z=at(tuple(tuple((b,a,I) for a,b,I in t) for t in z))
                    if s:out[z]=(out.get(z,0)+c*k*s*fac)%P
        return canon(out)
    def native(x,mode='raw'):
        inp=out/'control.bin'; stem=out/'result';write_records(inp,{key(w):c for w,c in x.items()});subprocess.run([exe,alpha,inp,stem,mode],check=True,capture_output=True)
        return {op:{unpack(k):c for k,c in records(Path(str(stem)+'.'+op+'.bin'))} for op in ops}
    rng=random.Random(2026112001);pool=[l for l in LETTERS if sum(l[:2])<=1];checked=odd_pairs=0;supports={op:0 for op in ops};anticomm={op:0 for op in ops}
    while checked<18:
        shapes=(3,3) if checked<9 else (2,3) if checked<14 else (4,)
        s,w=at(tuple(tuple(rng.choice(pool) for _ in range(n)) for n in shapes))
        if not s or (checked<9 and not all(word_parity(t) for t in w)):continue
        x={w:s};qx=canon(Q_expr(x));
        if not qx:continue
        want={op:baseline(x,op) for op in ops};qwant={op:canon(Q_expr(z)) for op,z in want.items()}
        if not all(qwant.values()):continue
        dx=native(x);dqx=native(qx);assert dx==want
        for op in ops:
            assert dqx[op]==baseline(qx,op)==canon({z:-c for z,c in qwant[op].items()});supports[op]+=len(dx[op]);anticomm[op]+=len(dqx[op])
        odd_pairs+=int(len(w)==2 and all(word_parity(t) for t in w));checked+=1
    # Also replay the extremal projection with the independent Python implementation, then the public maps.
    import pickle
    sys.path.insert(0,str(RUNTIME));import su_trace_adapter_v1 as alg
    pf=HERE.parent/'su3_j30_classical_ranks/joint_plan.json';assert sha(pf)==PLAN;plan=json.loads(pf.read_bytes());tested=[]
    for i in (0,12,35):
        f=PACKAGE/'second_basis/shard0'/f'targets{i:02d}/independent_seeds.pkl';assert sha(f)==plan['sources'][i]['recipe_sha256'];w=pickle.loads(f.read_bytes())[0];x=alg.parity_project({w:1},-1);px=alg.projection.project(x,(4,4,14,8,8),P)
        want={op:baseline(px,op) for op in ops};assert native(x,'parent')==want;tested.append([i,len(px),*[len(want[op]) for op in ops]])
    assert all(any(row[k] for row in tested) for k in (2,3))
    return dict(prime=P,public_commit='4256f207986029540060fb44dd494746828f65ce',pins=pin,nonzero_controls=checked,two_odd_trace_controls=odd_pairs,word_coefficients_compared=supports,nonzero_exact_image_coefficients_compared=anticomm,projected_recipe_control_counts=tested,control_projection_column_order=['parent_profile','projected','S20','Tp21'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('upstream',type=Path);p.add_argument('executable',type=Path);p.add_argument('work',type=Path);args=p.parse_args();res=run(args.upstream.resolve(),args.executable.resolve(),args.work);(args.work/'controls.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
