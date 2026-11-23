"""The complete first-identity space one rank above the Sp(6) candidate."""
from fast_paths import *
from flint import fmpz_mat
import hashlib,json,math,pickle,time
import abelian_obstruction as ab

def main():
    cache_native_canonicalization();start=time.monotonic()
    rank=7;weight=(3,3,11,0,0);length=16
    from trace_descent import dimensions
    records=[]
    for character_rank in (7,17):
        value=sum(dimensions(weight,length,2*character_rank,'Sp').values())
        records.append(dict(rank=character_rank,native_weight=weight,length=length,
                            dimension=value,seconds=time.monotonic()-start))
    assert records[1]['dimension']-records[0]['dimension']==4
    (HERE/'next_rank_character.json').write_text(json.dumps(dict(
        records=records,next_rank_identity_dimension=4),indent=2)+'\n')
    basis=sorted({p for content in pf.m.contents(weight,length)
                  for p in pc.paired_content(tuple(sorted(content)))})
    images=[pc.raising_pairs(p) for p in basis]
    words=sorted({w for row in images for w in row});ids={w:i for i,w in enumerate(words)}
    raising=fmpz_mat(len(words),len(basis))
    for j,row in enumerate(images):
        for w,c in row.items():raising[ids[w],j]=c
    kernel,dimension=raising.nullspace();assert dimension==4
    rows=[];recipes=[];expand=lru_cache(maxsize=256)(expand_paths)
    for j in range(dimension):
        co=[int(kernel[i,j]) for i in range(len(basis))];d=math.gcd(*co)
        recipe=[(basis[i],c//d) for i,c in enumerate(co) if c]
        row=Counter()
        for paths,c in recipe:
            for w,n in expand(paths,rank).items():row[w]+=c*n
        row=pf.primitive(clean(row));assert not pf.spin_action(row,True)
        rows.append(row);recipes.append(recipe)
        print('next rank relation',j+1,len(row),'seconds',time.monotonic()-start,flush=True)
    rowcert=ab.rank_certificate(rows,1000003);assert rowcert['rank']==4
    outgoing=[pf.m.q(row) for row in rows]
    qcert=ab.rank_certificate(outgoing,1000003)
    qcoordinates=ab.rational_rank_upper(outgoing,qcert)
    minor=fmpz_mat([[outgoing[i].get(w,0) for w in qcert['columns']] for i in qcert['rows']])
    determinant=int(minor.det());assert determinant
    restrictions=[ab.restrict(row) for row in rows]
    acert=ab.rank_certificate(restrictions,1000003)
    aminor=fmpz_mat([[restrictions[i].get(w,0) for w in acert['columns']] for i in acert['rows']])
    adet=int(aminor.det());assert adet
    combined=[{**{('Q0',w):c for w,c in q.items()},**{('Cartan',w):c for w,c in a.items()}}
              for q,a in zip(outgoing,restrictions)]
    ccert=ab.rank_certificate(combined,1000003)
    cminor=fmpz_mat([[combined[i].get(w,0) for w in ccert['columns']] for i in ccert['rows']])
    cdet=int(cminor.det());assert cdet
    target=HERE/'next_rank_relations.pkl';target.write_bytes(pickle.dumps((rows,recipes),protocol=4))
    result=dict(group='Sp(7)',native_weight=weight,length=length,
        pair_space=len(basis),raised_space=len(words),highest_weight_relations=4,
        Q0_rank=qcert['rank'],integer_Q0_minor_determinant=str(determinant),
        exact_Q0_image_coordinates=qcoordinates,
        Q0_minor_columns=[pf.m.encode(w) for w in qcert['columns']],
        no_nonzero_closed_identity=qcert['rank']==4,
        stable_abelianization_rank=acert['rank'],
        abelianization_integer_minor_determinant=str(adet),
        abelianization_minor_columns=[pf.m.encode(w) for w in acert['columns']],
        exact_fortuitous_dimension_at_source_length_15=0 if acert['rank']==4 else None,
        combined_Q0_abelianization_rank=ccert['rank'],
        combined_integer_minor_determinant=str(cdet),
        combined_minor_columns=[(kind,pf.m.encode(w)) for kind,w in ccert['columns']],
        no_nonzero_closed_and_abelian_zero_identity=ccert['rank']==4,
        character_completeness_receipt='next_rank_character.json',
        polynomial_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        seconds=time.monotonic()-start)
    (HERE/'next_rank_relations.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('Q0_minor_columns','combined_minor_columns')}),flush=True)

if __name__=='__main__':main()
