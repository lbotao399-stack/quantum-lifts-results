"""Verify the archived B05 certificate; full source reconstruction is reproduce.py."""
import argparse,gzip,hashlib,json,resource,time
from fractions import Fraction
from pathlib import Path
from flint import nmod_mat
from characters import singlets
from maps import sparse_combine
from extract import modular,rational
from relation_test import rank_receipt

HERE=Path(__file__).resolve().parent


def read(path):return json.loads(gzip.decompress(path.read_bytes()) if path.suffix==".gz" else path.read_text())
def rows(values):return [{j:Fraction(n,d) for j,n,d in row} for row in values]
def row(values):return {j:Fraction(n,d) for j,n,d in values}


def manifest_check(directory=None):
    manifest=read(HERE/"MANIFEST.json")
    for entry in manifest["files"]:
        relative=Path(entry["path"])
        path=(directory/Path(*relative.parts[1:]) if directory is not None and relative.parts[0]=="certificate" else HERE/relative)
        assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==entry["sha256"],entry["path"]
    root=HERE.parents[1]
    for entry in manifest["repository_dependencies"]:
        path=root/entry["path"]
        assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==entry["sha256"],entry["path"]
    return len(manifest["files"]),len(manifest["repository_dependencies"])


def finite_coordinates(frames,total):
    result=[None]*total;basis=[];offset=0
    for block in frames:
        basis.extend(block["indices"][j] for j in block["basis_rows"])
        for j,values in zip(block["indices"],block["coordinates"]):result[j]={offset+k:c for k,c in row(values).items()}
        offset+=len(block["basis_rows"])
    assert all(value is not None for value in result)
    assert all(result[j]=={i:Fraction(1)} for i,j in enumerate(basis))
    return basis,result


def verify(directory,rational_ranks=False,check_manifest=True):
    start=time.monotonic();manifest=manifest_check(directory) if check_manifest else None
    basis=read(directory/"basis.json.gz");frames=read(directory/"frames.json.gz");maps=read(directory/"maps.json.gz");g=read(directory/"generator.json.gz")
    finite={};coordinate={};frame_count=0
    for length,stable,expected in ((6,2091,1169),(7,4991,1896),(8,6609,1588)):
        assert len(basis["cochains"][str(length)])==singlets(length,length)["dotted_singlets"]==stable
        assert singlets(length,3)["dotted_singlets"]==expected
        chosen,coordinate[length]=finite_coordinates(frames[str(length)],stable)
        assert len(chosen)==expected
        indices=[]
        for block in frames[str(length)]:
            indices.extend(block["indices"]);minor=nmod_mat(block["minor"],block["prime"])
            assert minor.nrows()==minor.ncols()==len(block["basis_rows"])
            assert int(minor.det())==block["determinant_mod_prime"]!=0
            assert len(block["columns"])==minor.nrows() and len(set(block["columns"]))==minor.nrows()
            frame_count+=1
        assert sorted(indices)==list(range(stable))
        if length in (6,7):assert chosen==maps[str(length)]["source_basis"]
        if length in (7,8):assert chosen==maps[str(length-1)]["target_basis"]
    incoming=rows(maps["6"]["finite"]);outgoing=rows(maps["7"]["finite"])
    stable_incoming=rows(maps["6"]["stable"]);stable_outgoing=rows(maps["7"]["stable"])
    proper=rows(g["finite_proper"]);products=rows(g["finite_products"]);stable_products=rows(g["stable_products"])
    candidate=row(g["representative"]);dual=row(g["dual"])
    pair=lambda value:sum(c*dual.get(j,0) for j,c in value.items())
    assert all(pair(value)==0 for value in incoming+proper+products) and pair(candidate)==1
    assert not sparse_combine(candidate,outgoing)
    for value in stable_incoming+stable_products:assert not sparse_combine(value,stable_outgoing)
    closure=read(directory/"closure_relation.json");residual=row(closure["stable_image"])
    lifted={maps["7"]["source_basis"][j]:c for j,c in candidate.items()}
    assert sparse_combine(lifted,stable_outgoing)==residual
    assert residual and not sparse_combine(residual,coordinate[8])
    target_basis=maps["7"]["target_basis"];relation_coefficients=row(closure["target_relation_coefficients"])
    reconstructed=dict(relation_coefficients)
    for j,c in relation_coefficients.items():
        assert j not in target_basis
        for k,v in coordinate[8][j].items():reconstructed[target_basis[k]]=reconstructed.get(target_basis[k],0)-c*v
    assert {j:c for j,c in reconstructed.items() if c}==residual
    ranks={"incoming":modular(incoming,1896).rank(),"outgoing":modular(outgoing,1588).rank(),
           "boundaries_and_proper_products":modular(incoming+proper,1896).rank(),
           "closed_span_with_representative":modular(incoming+proper+[candidate],1896).rank()}
    assert ranks==dict(incoming=825,outgoing=1037,boundaries_and_proper_products=858,closed_span_with_representative=859)
    saved=read(directory/"rank_certificates.json")
    for key,values in (("stable_outgoing",stable_outgoing),("stable_closed_span",stable_incoming+stable_products),
                       ("stable_closed_finite_image",incoming+products)):
        receipt=rank_receipt(values,label=key)
        assert receipt==saved[key]
    assert saved["stable_outgoing"]["rank"]+saved["stable_closed_span"]["rank"]==4991
    assert saved["stable_closed_finite_image"]["rank"]==858
    if rational_ranks:
        assert rational(incoming,1896).rank()==825 and rational(outgoing,1588).rank()==1037
    return dict(verified=True,manifest_counts=manifest,complete_finite_frames=frame_count,
                complete_source_relations=3095,complete_target_relations=5021,finite_ranks=ranks,
                representative_pairing=1,all_boundary_pairings_zero=True,all_proper_product_pairings_zero=True,
                full_stable_kernel=1643,stable_kernel_finite_image=858,
                complete_relation_quotient="2311 - 2310 = 1",rational_rank_recalculation=rational_ranks,
                source_rebuilt=False,scope="Archived exact coefficients, complete frame minors, exact map identities and all decisive ranks; source rebuilding is reproduce.py.",
                seconds=time.monotonic()-start,maximum_resident_set_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory",type=Path,nargs="?",default=HERE/"certificate")
    parser.add_argument("--rational-ranks",action="store_true");parser.add_argument("--skip-manifest",action="store_true")
    parser.add_argument("--receipt",type=Path);args=parser.parse_args();assert __debug__
    result=verify(args.directory,args.rational_ranks,not args.skip_manifest)
    if args.receipt:args.receipt.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
