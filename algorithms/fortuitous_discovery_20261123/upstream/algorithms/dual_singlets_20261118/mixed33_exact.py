"""Rational rank and cycle quotient, with independently checkable sparse certificates.

Matrices have source vectors in rows. In a dependence relation the coefficient of
its last input row is nonzero. Distinct last rows certify independence over Q.
The verifier performs direct integer combinations, not Gaussian elimination.
"""
import argparse, gzip, hashlib, json, math
from pathlib import Path

HERE=Path(__file__).resolve().parent

def read(path):
    data=path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix=='.gz' else data)

def save(x,path):
    b=(json.dumps(x,separators=(',',':'))+'\n').encode()
    path.write_bytes(gzip.compress(b,mtime=0) if path.suffix=='.gz' else b)

def sha(path):return hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()

def linear(x,a,y,b):
    r={i:a*c for i,c in x.items()}
    for i,c in y.items():r[i]=r.get(i,0)+b*c
    return {i:c for i,c in r.items() if c}

def divide(d,g):return {i:c//g for i,c in d.items()}

def exact(rows):
    pivots={};deps=[];inds=[]
    for j,x in enumerate(rows):
        r=x.copy();v={j:1}
        while r and (key:=min(r)) in pivots:
            e,h=pivots[key];t,s=e[key],r[key];g=math.gcd(t,s);t//=g;s//=g
            r=linear(r,t,e,-s);v=linear(v,t,h,-s)
            g=math.gcd(*r.values(),*v.values());r=divide(r,g);v=divide(v,g)
        if r:
            pivots[min(r)]=r,v;inds.append(j)
        else:deps.append(v)
    return inds,sorted(pivots),deps

def pairs(x):return sorted(x.items())

def apply(v,rows):
    ans={}
    for j,c in v.items():
        for k,d in rows[j].items():ans[k]=ans.get(k,0)+c*d
    return {k:c for k,c in ans.items() if c}

def integer_rows(raw,n):
    assert all(len(dict(r))==len(r) and all(type(k)==type(c)==int and 0<=k<n and c!=0 for k,c in r) for r in raw)
    return list(map(dict,raw))

def key(w):return tuple(tuple((i,j,tuple(s)) for i,j,s in t) for t in w)

def input_rows(inc,out):
    words=sorted(set(map(key,inc['target_dictionary']))|set(map(key,out['source_dictionary'])));idx={w:i for i,w in enumerate(words)}
    def lift(x,ws):
        ids=[idx[key(w)] for w in ws];return [{ids[i]:c for i,c in r.items()} for r in integer_rows(x,len(ids))]
    return lift(inc['images'],inc['target_dictionary']),lift(out['sources'],out['source_dictionary']),integer_rows(out['images'],len(out['target_dictionary']))

def build(ip,op):
    inc,out=read(ip),read(op);incoming,sources,outgoing=input_rows(inc,out)
    I,J,dep=exact(incoming);a,b,ker=exact(outgoing)
    closed=[apply(v,sources) for v in ker];sel,piv,_=exact([incoming[i] for i in I]+closed);co=[j-len(I) for j in sel if j>=len(I)]
    return dict(prime_free=True,incoming_plan=sha(ip),outgoing_plan=sha(op),incoming_shape=[len(incoming),len(inc['target_dictionary'])],incoming_nnz=sum(map(len,incoming)),incoming_rank=len(I),incoming_null_relations=[pairs(x) for x in dep],incoming_basis_rows=I,incoming_pivot_columns=J,outgoing_shape=[len(outgoing),len(out['target_dictionary'])],outgoing_rank=len(a),outgoing_basis_rows=a,outgoing_columns=b,closed_singlet_recipes=[pairs(x) for x in ker],formal_quotient_dimension=len(co),formal_class_recipes=[pairs(ker[j]) for j in co])

def check(ip,op,cert):
    inc,out=read(ip),read(op);incoming,sources,outgoing=input_rows(inc,out)
    assert sha(ip)==cert['incoming_plan'] and sha(op)==cert['outgoing_plan']
    for rows,recipes,idkey,rank in [(incoming,cert['incoming_null_relations'],'incoming_basis_rows',225),(outgoing,cert['closed_singlet_recipes'],'outgoing_basis_rows',85)]:
        rs=list(map(dict,recipes));last=[max(x) for x in rs];assert len(rs)==len(set(last))==len(rows)-rank and set(range(len(rows)))-set(last)==set(cert[idkey]);assert all(not apply(x,rows) for x in rs)
        # Independent exact lower bound using a rational/integer minor, formed from the original rows.
        from flint import fmpz_mat
        # For incoming the columns refer to the union dictionary ordered in input_rows().
        cols=cert['incoming_pivot_columns' if rows is incoming else 'outgoing_columns'];det=int(fmpz_mat([[rows[i].get(j,0) for j in cols] for i in cert[idkey]]).det());assert len(cols)==rank and det
        print('rank',rank,'integer minor',det,'relations',len(rs))
    classes=[apply(dict(v),sources) for v in cert['formal_class_recipes']];assert all(not apply(dict(v),outgoing) for v in cert['formal_class_recipes']);rs=[incoming[i] for i in cert['incoming_basis_rows']]+classes
    from flint import fmpz_mat
    ii,jj,_=exact(rs);assert len(ii)==len(rs)==227 and int(fmpz_mat([[r.get(j,0) for j in jj] for r in rs]).det())
    assert len(sources)-85-225==cert['formal_quotient_dimension']==len(classes)==2

if __name__=='__main__':
    arg=argparse.ArgumentParser(description=__doc__);arg.add_argument('--write',action='store_true');args=arg.parse_args();assert __debug__
    ip=HERE/'mixed33_incoming.json.gz';op=HERE/'mixed33/outgoing.json.gz';dst=HERE/'mixed33/formal_certificate.json.gz'
    if args.write:save(build(ip,op),dst)
    check(ip,op,read(dst))
