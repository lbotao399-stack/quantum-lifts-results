"""Write the complete mathematical certificate as deterministic JSON archives."""
import argparse,gzip,hashlib,json,pickle
from fractions import Fraction
from pathlib import Path


def read(path):return pickle.loads(path.read_bytes())
def entries(row):return [[int(j),Fraction(c).numerator,Fraction(c).denominator] for j,c in sorted(row.items()) if c]
def save(path,value):
    raw=(json.dumps(value,separators=(",",":"),sort_keys=True)+"\n").encode()
    path.write_bytes(gzip.compress(raw,mtime=0) if path.suffix==".gz" else raw)
    print(json.dumps(dict(file=path.name,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())),flush=True)


def build(work,out):
    out.mkdir(parents=True,exist_ok=True)
    plans={length:read(work/f"stable{length}.pkl") for length in (6,7,8)}
    all_words=sorted({word for plan in plans.values() for row in plan["sources"] for word in row})
    word_ids={word:i for i,word in enumerate(all_words)}
    letters=sorted({letter for word in all_words for trace in word for letter in trace});letter_ids={letter:i for i,letter in enumerate(letters)}
    encode=lambda row:entries({word_ids[word]:coefficient for word,coefficient in row.items()})
    maps={length:read(work/f"map{length}.pkl") for length in (6,7)}
    generator=read(work/"generator.pkl");relation=read(work/"relation_test.pkl")
    if (work/"polynomial.pkl").exists():polynomial=read(work/"polynomial.pkl")
    else:
        from trace_algebra import combine,differential
        row=combine(generator["candidate"],[plans[7]["sources"][j] for j in maps[7]["source_basis"]])
        polynomial=dict(polynomial=row,residual=differential(row))
    basis=dict(schema=1,letters=[[a,b,sum(1<<t for t in theta)] for a,b,theta in letters],
               words=[[[letter_ids[letter] for letter in trace] for trace in word] for word in all_words],
               cochains={str(length):[encode(row) for row in plans[length]["sources"]] for length in (6,7,8)})
    save(out/"basis.json.gz",basis)
    encoded_maps={str(length):dict(source_basis=maps[length]["source_basis"],target_basis=maps[length]["target_basis"],
                                  stable=[entries(row) for row in maps[length]["stable"]],finite=[entries(row) for row in maps[length]["finite"]]) for length in (6,7)}
    save(out/"maps.json.gz",encoded_maps)
    frames={}
    for length in (6,7,8):
        blocks=[]
        for file in sorted((work/f"frame{length}").glob("frame_*.pkl")):
            block=read(file)
            blocks.append(dict(skeleton=block["skeleton"],indices=block["indices"],alphabet=[letter_ids[letter] for letter in block["alphabet"]],
                               points=block["points"],seed=block["seed"],degree=block["degree"],columns=block["columns"],basis_rows=block["basis_rows"],
                               determinant_mod_prime=block["determinant_mod_prime"],prime=block["prime"],bounds=block["bounds"],
                               minor=block["minor"].tolist(),coordinates=[entries({j:Fraction(c) for j,c in enumerate(row) if c!="0"}) for row in block["coordinates"]]))
        frames[str(length)]=blocks
    save(out/"frames.json.gz",frames)
    g=dict(schema=1,weight_flavour_then_dotted=[3,0,0,4,4],length=7,charge_conjugation=-1,dotted_spin=0,
           physical_doubled_charges_J1_J2_q1_q2_q3=[4,4,14,8,8],level=30,
           differential_convention="Divided jets, Q(Psi)=Psi^2. This Q is minus the pinned undivided native Q after factorial jet rescaling.",
           representative=entries(generator["candidate"]),dual=entries(generator["dual"]),
           full_polynomial=encode(polynomial["polynomial"]),full_stable_residual=encode(polynomial["residual"]),
           stable_products=[entries(row) for row in generator["stable_products"]],
           stable_proper=[entries(row) for row in generator["stable_proper"]],
           finite_products=[entries(row) for row in generator["finite_products"]],
           finite_proper=[entries(row) for row in generator["finite_proper"]],
           closed_span_independent_indices=generator["selected"],dual_columns=generator["dual_columns"],
           complete_source_relations=3095,complete_target_relations=5021,
           stable_outgoing_rank=3348,stable_closed_span_rank=1643,stable_closed_finite_image_rank=858,
           target_relations_in_stable_image=2311,image_of_all_source_relations=2310,fortuitous_quotient=1)
    save(out/"generator.json.gz",g)
    save(out/"rank_certificates.json",{key:relation[key] for key in ("stable_outgoing","stable_closed_span","stable_closed_finite_image")})
    archived=work/"algorithms/su3_j30_classical_ranks/instanton_probe/outgoing_products.json"
    if not archived.exists():archived=Path(__file__).resolve().parent/"certificate/proper_product_recipes.json"
    save(out/"proper_product_recipes.json",json.loads(archived.read_text()))
    products=read(work/"products.pkl")
    save(out/"complete_graviton_recipes.json",products["complete_recipes"])
    save(out/"summary.json",dict(source_stable=4991,source_finite=1896,source_relations=3095,target_stable=6609,target_finite=1588,target_relations=5021,
                                 representative_finite_terms=len(generator["candidate"]),representative_trace_terms=len(polynomial["polynomial"]),
                                 dual_nonzero_coordinates=len(generator["dual"]),stable_residual_terms=len(polynomial["residual"]),
                                 incoming_rank=825,outgoing_rank=1037,cohomology=34,proper_products=33,ring_generator_quotient=1,
                                 stable_outgoing_rank=3348,stable_kernel=1643,stable_cycles_finite_image=858,source_relations_closed=785,
                                 image_of_source_relations=2310,target_relations_in_stable_image=2311,fortuitous_quotient=1,
                                 candidate_pairing=1,boundary_pairings_nonzero=0,proper_product_pairings_nonzero=0,stable_graviton_pairings_nonzero=0))
    from maps import sparse_combine,finite_coordinates
    source={maps[7]["source_basis"][i]:c for i,c in generator["candidate"].items()}
    residual=sparse_combine(source,maps[7]["stable"])
    target_basis,coordinates=finite_coordinates(work/"frame8",6609)
    assert not sparse_combine(residual,coordinates)
    nonbasis=set(range(6609))-set(target_basis)
    relation_coefficients={j:c for j,c in residual.items() if j in nonbasis}
    reconstructed=dict(relation_coefficients)
    for j,c in relation_coefficients.items():
        for k,v in coordinates[j].items():reconstructed[target_basis[k]]=reconstructed.get(target_basis[k],0)-c*v
    assert {j:c for j,c in reconstructed.items() if c}==residual
    save(out/"closure_relation.json",dict(stable_image=entries(residual),target_relation_coefficients=entries(relation_coefficients),
         target_relation_convention="For every nonbasis stable row j, relation_j = row_j - sum_k coordinates[j,k] finite_basis_row_k.",exact_finite_residual=0))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();build(args.work,args.output)
