"""Full source Q0 equality, partitioned by its output bare-ghost color.

All source coefficients, q0field images, canonical signs and odd-prefix signs
are preserved. Splitting disjoint output monomials bounds peak memory. Integer
denominators are cleared exactly; zero accumulator entries are removed at once.
"""
from fractions import Fraction as F
from math import lcm
import gc

def verify(theory,polynomial,expected,require_zero=True,verbose=False):
    source_den=1;target_den=1;letters=set()
    for word,c in polynomial.items():
        assert not any(u[0][0]=='c' and u[1:]==(0,0) for u in word)
        source_den=lcm(source_den,c.denominator);letters.update(word)
    for c in expected.values():target_den=lcm(target_den,c.denominator)
    raw={u:theory.q0field(u) for u in letters};field_den=1
    colors=sorted({v[0] for image in raw.values() for word in image for v in word if v[0][0]=='c' and v[1:]==(0,0)})
    color_groups=15 if len(polynomial)>1000000 else 4
    colors={color:i%color_groups for i,color in enumerate(colors)}
    def group(word):
        bare=[u for u in word if u[0][0]=='c' and u[1:]==(0,0)]
        assert len(bare)<=1
        return colors[bare[0][0]] if bare else color_groups
    for image in raw.values():
        for c in image.values():field_den=lcm(field_den,c.denominator)
    denominator=lcm(source_den*field_den,target_den)
    source_factor=denominator//(source_den*field_den)
    images=[{} for _ in range(color_groups+1)]
    for u,image in raw.items():
        for word,c in image.items():
            images[group(word)].setdefault(u,[]).append((word,c.numerator*(field_den//c.denominator)))
    records=[]
    for channel in range(color_groups+1):
        residual={word:-c.numerator*(denominator//c.denominator) for word,c in expected.items() if group(word)==channel}
        peak=len(residual)
        for word,c in polynomial.items():
            value=c.numerator*(source_den//c.denominator)*source_factor;prefix=1
            for i,u in enumerate(word):
                for image,co in images[channel].get(u,()):
                    key,sign=theory.canon(word[:i]+image+word[i+1:])
                    if sign:
                        result=residual.get(key,0)+value*prefix*sign*co
                        if result:residual[key]=result
                        else:residual.pop(key,None)
                if theory.parity(u):prefix=-prefix
            peak=max(peak,len(residual))
        records.append({'channel':channel,'peak_accumulator_terms':peak,'residual_terms':len(residual)})
        if verbose:print('full Q0 output channel',records[-1],flush=True)
        if require_zero:assert not residual,('full Q0 mismatch',records[-1])
        del residual;gc.collect()
    return records

if __name__=='__main__':
    import os,sys,pickle,json
    from pathlib import Path
    os.environ['SUN']='4';root=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(root/'rational_c4'))
    import exact_product as E
    tests=[E.mul(E.v('be12',1,0),E.v('ch1',0,1)),
           E.mul(E.v('be21',0,1),E.v('ch2',1,0),E.v('ce12',1,0))]
    saved=pickle.load(open(Path(__file__).with_name('k2_full.pkl'),'rb'))
    tests.extend([saved['leading'],saved['first']])
    for p in tests:verify(E.T,p,E.T.q0(p))
    negative=verify(E.T,tests[0],{},require_zero=False)
    assert sum(r['residual_terms'] for r in negative)>0
    result={'status':'PASS','full_source_Q0_comparisons':len(tests),'noninvariant_bare_ghost_sectors_tested':True,'negative_control':negative}
    Path(__file__).with_name('full_q0_stream_gates.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
