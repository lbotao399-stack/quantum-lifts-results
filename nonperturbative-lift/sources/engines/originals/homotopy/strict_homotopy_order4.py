"""Triangular, rotation-equivariant extension of the b-derivative homotopy.

Only local work files are written. Full rational equations are checked after
each solve. An unsuccessful obstruction solve is saved, never promoted.
"""
from invariant_covariants import *
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
from math import factorial
import exact_star as STAR

C.CHECK_PATH=ROOT/'strict_homotopy_order4_checks.json'
OUT=ROOT/'order4';OUT.mkdir(exist_ok=True)

def local_check(name,residual):
    if residual:
        write('failed_residual',dict(name=name,residual=residual))
        raise AssertionError((name,len(residual)))
    print(name,'PASS',flush=True)
C.check=local_check

for kind,partner in [('q','p'),('r','s'),('s','r')]:
    E.S.UP[kind]=True
    for color in E.S.labels:
        name=kind+color
        T.info[name]=dict(par=0,wt=0,kind='probe')
        T.partners[name]=partner+color
        T.base[name]={};T.wt[name]=0
T.fields=tuple(sorted(T.info))

initial=pickle.load(open(ROOT/'strict_h3_bare_triples.pkl','rb'))
C0,x=initial['C0'],initial['x']
tables={1:{},2:dict(initial['h2_pairs']),3:{},4:{}}
solved={1:set(),2:set(),3:set(),4:set()}
for key in tables[2]:solved[2].add(tuple(sorted(u[1:] for u in key)))

def write(name,value):
    path=OUT/(name+'.pkl');tmp=path.with_suffix('.tmp')
    with open(tmp,'wb') as file:pickle.dump(value,file)
    tmp.replace(path)

def sector_name(n,jets):
    return 'h%d_'%n+'_'.join('%d-%d'%j for j in jets)

def encoding(jets):
    groups=sorted(set(jets))
    assert len(groups)<=4
    mapping=dict(zip(groups,'pqrs'))
    poly={():F(1)}
    for jet in jets:
        p=mapping[jet]
        factor=E.add(*(E.mul(T.v('b'+color,*jet),T.v(p+color)) for color in E.S.labels))
        poly=E.mul(poly,factor)
    return poly,mapping

def install(n,jets,value,mapping=None):
    if mapping is None:mapping=encoding(jets)[1]
    reverse={p:j for j,p in mapping.items()}
    local={}
    for word,c in value.items():
        ps=[u for u in word if u[0][0] in reverse]
        ghosts=tuple(u for u in word if u[0][0] not in reverse)
        assert len(ps)==n and all(u[0][0]=='c' for u in ghosts)
        key=tuple(sorted(('b'+u[0][1:],*reverse[u[0][0]]) for u in ps))
        multiplicity=1
        for kind in reverse:
            group=[u for u in ps if u[0][0]==kind]
            count=Counter(group)
            multiplicity*=factorial(len(group))
            for k in count.values():multiplicity//=factorial(k)
        coefficient=c/multiplicity
        assert ghosts not in local.setdefault(key,{})
        local[key][ghosts]=coefficient
    tables[n].update(local);solved[n].add(tuple(sorted(jets)))

for m in range(5):
    for jet,value in pickle.load(open(ROOT/('strict_eta_jet%d.pkl'%m),'rb')).items():
        install(1,(jet,),value)
install(3,((0,0),)*3,initial['value'])

def h(n,poly):
    if n==0:return E.mul(x,poly)
    out={}
    for word,c in poly.items():
        bs=[i for i,u in enumerate(word) if u[0][0]=='b']
        for positions in combinations(bs,n):
            key=tuple(word[i] for i in positions)
            jets=tuple(sorted(u[1:] for u in key))
            if jets not in solved[n]:raise KeyError(('unsolved',n,jets))
            coefficient=tables[n].get(key,{})
            if coefficient:
                rest=tuple(u for i,u in enumerate(word) if i not in positions)
                out=E.add(out,E.mul(coefficient,{rest:c}))
    return out

@lru_cache(maxsize=65536)
def l_coefficient(n,key):
    return STAR.star_n(n,C0,{key:F(1)})

def L(n,poly):
    if n==0:return E.mul(C0,poly)
    out={}
    for word,c in poly.items():
        bs=[i for i,u in enumerate(word) if u[0][0]=='b']
        for positions in combinations(bs,n):
            key=tuple(word[i] for i in positions)
            coefficient=l_coefficient(n,key)
            if coefficient:
                rest=tuple(u for i,u in enumerate(word) if i not in positions)
                out=E.add(out,E.mul(coefficient,{rest:c}))
    return out

def Q(n,poly):
    return E.Q0(poly) if n==0 else T.qn(n,poly)

def tagged(a,b):
    return {(0,w):c for w,c in a.items()}|{(1,w):c for w,c in b.items()}

def fast_q_gate():
    sys.path.insert(0,str(ROOT.parent/'su3_rational'))
    import fast_wick
    old,new=fast_wick.install(T)
    cases=[]
    for n in range(1,4):
        for c in E.S.labels:
            poly=E.mul(x,*([T.v('b'+c)]*n))
            T.wick=old;T.selected.cache_clear();a=T.qn(n,poly)
            T.wick=new;T.selected.cache_clear();b=T.qn(n,poly)
            assert a==b
            cases.append(dict(order=n,color=c,terms=len(a)))
    write('fast_wick_gate',cases)
    print('Exact fast-Q gate PASS',flush=True)

def install_moment_projection():
    folder=ROOT.parent/'rational_c4'/'m6_lift'
    certificate=json.loads((folder/'adapter_gates.json').read_text())
    assert certificate['passed'] and certificate['full_compact_images']==[1,2,3]
    sys.path.insert(0,str(folder))
    import fast_selected
    original,optimized=fast_selected.install(T)
    # A fresh full-image comparison on the current source x and mixed colors.
    poly=E.mul(x,T.v('bh1'),T.v('be12'),T.v('be21'))
    original.cache_clear();optimized.cache_clear();T.selected=original
    a=T.qn(3,poly)
    original.cache_clear();optimized.cache_clear();T.selected=optimized
    b=T.qn(3,poly)
    assert a==b and a
    write('moment_projection_gate',dict(order=3,terms=len(a),passed=True))
    print('Exact active-output moment gate PASS',len(a),flush=True)

def extend_eta():
    eta={}
    for m in range(7):
        path=ROOT/('strict_eta_jet%d.pkl'%m)
        if path.exists():
            jets=pickle.load(open(path,'rb'));eta.update(jets)
            for jet,value in jets.items():install(1,(jet,),value)
            continue
        poly,mapping=encoding(((m,0),))
        target=E.sub(E.sc(E.mul(C0,L(1,poly)),2),E.add(Q(1,E.mul(x,poly)),h(1,Q(0,poly))))
        source=basis(3+m,3,'p',True)
        images=[Q(0,p) for p in source];selected=N.indices(images)
        print('eta',m,'source',len(source),'rank',len(selected),'target',len(target),flush=True)
        co,_=C.rational_solve([images[i] for i in selected],target)
        value=E.add(*(E.sc(source[i],c) for i,c in zip(selected,co)))
        jets={(m,0):value};lower=value;fac=1
        for j in range(1,m+1):
            lower=E.S.Jm(lower);fac*=m-j+1;jets[(m-j,j)]=E.sc(lower,F(1,fac))
        for jet,v in jets.items():install(1,(jet,),v)
        for jet,v in jets.items():
            poly,_=encoding((jet,))
            assert E.add(Q(0,v),h(1,Q(0,poly)),Q(1,E.mul(x,poly)))==E.sc(E.mul(C0,L(1,poly)),2)
        with open(path,'wb') as file:pickle.dump(jets,file)
        eta.update(jets)
        print('eta',m,'PASS terms',len(value),'seconds',time.time()-START,flush=True)

def solve_sector(n,jets):
    name=sector_name(n,jets);path=OUT/(name+'.pkl')
    poly,mapping=encoding(jets)
    if path.exists():
        result=pickle.load(open(path,'rb'));install(n,jets,result['value'],mapping)
        print(name,'readback',len(result['value']),flush=True);return
    target_path=OUT/(name+'_target.pkl')
    if target_path.exists():
        target=pickle.load(open(target_path,'rb'))
    else:
        target=E.add(*(L(i,L(n-i,poly)) for i in range(n+1)))
        print(name,'L-square',len(target),flush=True)
        for r in range(1,n+1):
            target=E.sub(target,E.add(Q(r,h(n-r,poly)),h(n-r,Q(r,poly))))
            print(name,'Q',r,'target',len(target),flush=True)
        target=E.sub(target,h(n,Q(0,poly)))
        write(name+'_target',target)
    # Solve the affine Q0 equation first. Rotation is imposed on the completed
    # fixed-jet-degree map by Reynolds projection, not by greedy descendants.
    source=list(candidates(n+2+sum(j[0] for j in jets),n+2+sum(j[1] for j in jets),
                           ''.join(mapping[j] for j in jets)))
    images=[Q(0,p) for p in source]
    selected=N.indices(images)
    print(name,'candidates',len(source),'Q0 rank',len(selected),'target',len(target),flush=True)
    co,minor=C.rational_solve([images[i] for i in selected],target)
    value=E.add(*(E.sc(source[i],c) for i,c in zip(selected,co)))
    assert Q(0,value)==target
    install(n,jets,value,mapping)
    assert h(n,poly)==value
    result=dict(value=value,target=target,selected_source=[source[i] for i in selected],
                coefficients=co,minor=minor,jets=jets,arity=n,seconds=time.time()-START)
    write(name,result)
    print(name,'PASS terms',len(value),'seconds',time.time()-START,flush=True)

def rotate_map(values,keys,raising):
    out={}
    for key in keys:
        value=(E.Jp if raising else E.S.Jm)(values.get(key,{}))
        for i,u in enumerate(key):
            count=u[2] if raising else u[1]
            if not count:continue
            changed=(u[0],u[1]+(1 if raising else -1),u[2]+(-1 if raising else 1))
            other=tuple(sorted(key[:i]+(changed,)+key[i+1:]))
            value=E.sub(value,E.sc(values.get(other,{}),count))
        if value:out[key]=value
    return out

def project_total(n,total,sectors):
    path=OUT/('h%d_total%d_projected.pkl'%(n,total))
    if path.exists():
        projected=pickle.load(open(path,'rb'))['coefficients']
        for key in list(tables[n]):
            if sum(u[1]+u[2] for u in key)==total:del tables[n][key]
        tables[n].update(projected)
        for jets in sectors:solved[n].add(jets)
        print('Projected map readback',n,total,flush=True);return
    keys=sorted({tuple(sorted(('b'+c,*j) for c,j in zip(colors,jets)))
                 for jets in sectors for colors in product(E.S.labels,repeat=n)})
    values={key:tables[n][key] for key in keys if tables[n].get(key)}
    for j in range(1,n+3+total):
        raised=rotate_map(values,keys,True)
        if not raised:break
        lowered=rotate_map(raised,keys,False)
        values={key:E.sub(values.get(key,{}),E.sc(lowered.get(key,{}),F(1,j*(j+1)))) for key in keys}
        values={key:v for key,v in values.items() if v}
        print('Reynolds',n,total,'spin',j,'terms',sum(map(len,values.values())),flush=True)
    assert not rotate_map(values,keys,True)
    assert not rotate_map(values,keys,False)
    for key in keys:tables[n].pop(key,None)
    tables[n].update(values)
    # Check every stored affine equation after projection, not just covariance.
    for jets in sectors:
        poly,_=encoding(jets)
        target=pickle.load(open(OUT/(sector_name(n,jets)+'_target.pkl'),'rb'))
        assert Q(0,h(n,poly))==target,('projected Q0',n,jets)
    write('h%d_total%d_projected'%(n,total),dict(coefficients=values,total=total,arity=n,
                                               Jplus=True,Jminus=True,Q0=True))
    print('Fixed total jet degree projected and checked',n,total,'keys',len(values),flush=True)

def main():
    # Decode and operator-normalization checks against stored certificates.
    for n in (1,2,3):
        poly,_=encoding(((0,0),)*n)
        ref=pickle.load(open(ROOT/('strict_h%d_bare_%s.pkl'%(n,'pairs' if n==2 else 'triples')),'rb'))['value'] if n>1 else initial['eta'][(0,0)]
        assert h(n,poly)==ref
    test=E.mul(T.v('bh1',1,0),T.v('be12'),T.v('ch1',0,1))
    for n in (1,2):assert L(n,test)==STAR.star_n(n,C0,test)
    fast_q_gate();install_moment_projection();extend_eta()
    nstop=int(os.environ.get('HOMOTOPY_STOP_ARITY','4'))
    for n,limit in ((2,4),(3,2),(4,0)):
        if n>nstop:break
        for total in range(limit+1):
            jets=[(r,s) for r in range(total+1) for s in range(total+1-r)]
            sectors=[j for j in combinations_with_replacement(jets,n) if sum(sum(k) for k in j)==total]
            sectors.sort(key=lambda j:(-sum(k[0]-k[1] for k in j),j))
            if all(j in solved[n] for j in sectors):continue
            for j in sectors:
                if j not in solved[n]:solve_sector(n,j)
            project_total(n,total,sectors)
        write('h%d_coefficients'%n,tables[n])
    write('all_coefficients',dict(tables=tables,solved=solved,C0=C0,x=x))
    print('EXTENSION COMPLETED THROUGH ARITY',nstop,'seconds',time.time()-START,flush=True)

if __name__=='__main__':main()
