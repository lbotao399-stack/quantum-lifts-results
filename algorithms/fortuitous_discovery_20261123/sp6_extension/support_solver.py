"""Find a universal Q0 preimage using only the support of a given defect.

All successful outputs are verified over Q by recomputing the native Q0.
Failure of a restricted ansatz is never a non-existence certificate.
"""
from inverse_support import *
from flint import nmod_mat
from fractions import Fraction
import argparse
import math
import numpy as np

PRIME=2147483647


def use_highest_weights(defect):
    key=block(next(iter(defect)))
    active=[w for w in key if w[0] or w[1]]
    return len(active)==1 and active[0][:2]==(3,3) and not pf.spin_action(defect,True)


def spin_orbit(word):
    locations=[(t,i,l) for t,tr in enumerate(word) for i,l in enumerate(tr) if l[0]+l[1]]
    result=set()
    for choices in product(*(range(l[0]+l[1]+1) for _,_,l in locations)):
        if sum(choices)!=3:continue
        traces=[list(tr) for tr in word]
        for (t,i,l),a in zip(locations,choices):traces[t][i]=(a,l[0]+l[1]-a,l[2])
        sign,v=pf.m.so_canon(tuple(tuple(tr) for tr in traces))
        if sign:result.add(v)
    return result


def source_expressions(sources,highest):
    if not highest:return [{w:1} for w in sorted(sources)]
    representatives={pf.m.shape(w,'spin'):w for w in sources}
    result=[]
    for key,word in sorted(representatives.items()):
        result.extend(pf.m.singlets(sorted(spin_orbit(word))))
    return result


def rational_reconstruct(value,prime=PRIME):
    if not value:return Fraction(0)
    r,s=prime,int(value);t,u=0,1;bound=math.isqrt((prime-1)//2)
    while abs(s)>bound:
        q=r//s;r,s=s,r-q*s;t,u=u,t-q*u
    if not u or abs(u)>bound or math.gcd(s,u)!=1:return None
    answer=Fraction(s,u)
    return answer if (answer.numerator-value*answer.denominator)%prime==0 else None


def solve_block(defect,max_iterations=4,max_sources=12000):
    divisor=math.gcd(*defect.values())
    original=defect
    defect={w:c//divisor for w,c in defect.items()}
    highest=use_highest_weights(defect)
    sources=set().union(*(parents(w) for w in defect))
    records=[]
    for iteration in range(max_iterations):
        source=source_expressions(sources,highest)
        if len(source)>max_sources:
            return None,dict(reason='source cap',sources=len(source),raw_parents=len(sources),iterations=records)
        echelon=pf.Echelon(PRIME)
        for w in sorted(defect):echelon.words[w]=len(echelon.words)
        images=[];selected=[]
        for i,expression in enumerate(source):
            row=pf.m.q(expression);images.append(row)
            if echelon.add(row):selected.append(i)
        residual=np.zeros(len(echelon.words),dtype=np.int64)
        for w,c in defect.items():residual[echelon.words[w]]=c%PRIME
        for pivot,(indices,values) in sorted(echelon.rows.items()):
            factor=int(residual[pivot])
            if factor:residual[indices]=(residual[indices]-factor*values)%PRIME
        inverse={i:w for w,i in echelon.words.items()}
        missing={inverse[i] for i in np.flatnonzero(residual)}
        record=dict(iteration=iteration,sources=len(source),raw_parents=len(sources),
                    highest_weight_reduction=highest,targets=len(inverse),
                    rank=len(selected),residual_terms=len(missing))
        records.append(record)
        if not missing:
            columns=[inverse[i] for i in sorted(echelon.rows)]
            matrix=nmod_mat([[images[j].get(w,0)%PRIME for j in selected] for w in columns],PRIME)
            rhs=nmod_mat([[defect.get(w,0)%PRIME] for w in columns],PRIME)
            solution=matrix.solve(rhs)
            coefficients=[rational_reconstruct(int(c)) for c in solution.entries()]
            if any(c is None for c in coefficients):
                return None,dict(reason='rational reconstruction failed',iterations=records)
            denominator=math.lcm(*(c.denominator for c in coefficients))
            integral=Counter()
            for j,c in zip(selected,coefficients):
                if c:
                    for w,d in source[j].items():integral[w]+=int(c*denominator)*d
            integral=clean(integral)
            if pf.m.q(integral)!={w:c*denominator for w,c in defect.items()}:
                return None,dict(reason='exact native Q verification failed',iterations=records)
            return {w:Fraction(c*divisor,denominator) for w,c in integral.items()},dict(
                passed=True,exact_rational=True,operator_terms=len(integral),
                denominator=denominator,iterations=records)
        old=len(sources)
        sources.update(w for v in missing for w in parents(v))
        if len(sources)==old:
            return None,dict(reason='restricted support does not span defect',iterations=records)
    return None,dict(reason='iteration cap',iterations=records)


def main():
    cache_native_canonicalization()
    ap=argparse.ArgumentParser()
    ap.add_argument('--defect',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--max-sources',type=int,default=12000)
    ap.add_argument('--seconds',type=int,default=300)
    args=ap.parse_args();args.output.mkdir(exist_ok=True)
    defect=pickle.loads(args.defect.read_bytes())
    groups=defaultdict(dict)
    for word,c in defect.items():groups[block(word)][word]=c
    ordered=sorted(groups.items(),key=lambda kv:(len(kv[1]),kv[0]))
    start=time.monotonic();operator=Counter();records=[]
    for i,(key,part) in enumerate(ordered):
        path=args.output/f'block_{i}.pkl'
        digest=__import__('hashlib').sha256(pickle.dumps(part,protocol=4)).hexdigest()
        if path.exists():
            old,answer,record=pickle.loads(path.read_bytes());assert old==digest
            if answer is None and record.get('reason')=='source cap' and record['sources']<=args.max_sources:
                answer,record=solve_block(part,max_sources=args.max_sources)
                record.update(trace_weights=key,defect_terms=len(part))
                path.write_bytes(pickle.dumps((digest,answer,record),protocol=4))
        else:
            answer,record=solve_block(part,max_sources=args.max_sources)
            record.update(trace_weights=key,defect_terms=len(part))
            path.write_bytes(pickle.dumps((digest,answer,record),protocol=4))
        records.append(record)
        if answer is not None:
            for w,c in answer.items():operator[w]+=c
        if i%50==0 or len(part)>500 or answer is None or i+1==len(ordered):
            print(json.dumps(dict(block=i,total=len(ordered),passed=answer is not None,
                                  defect_terms=len(part),seconds=time.monotonic()-start,
                                  final_iteration=record.get('iterations',[])[-1:])),flush=True)
            (args.output/'progress.json').write_text(json.dumps(dict(completed=i+1,total=len(ordered),
                solved=sum(bool(r.get('passed')) for r in records),records=records),indent=2)+'\n')
        if time.monotonic()-start>args.seconds:
            print('Checkpoint time cap; rerun to continue.',flush=True);return
    passed=all(r.get('passed') for r in records)
    result=dict(passed=passed,blocks=len(records),solved=sum(bool(r.get('passed')) for r in records),
                operator_terms=len(operator),seconds=time.monotonic()-start,records=records)
    if passed:
        denominator=math.lcm(*(c.denominator for c in operator.values()))
        integral={w:int(c*denominator) for w,c in operator.items() if c}
        assert pf.m.q(integral)=={w:c*denominator for w,c in defect.items()}
        (args.output/'universal_preimage.pkl').write_bytes(pickle.dumps(dict(
            polynomial=integral,denominator=denominator),protocol=4))
        result.update(denominator=denominator,exact_native_Q_recomputed=True)
    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
