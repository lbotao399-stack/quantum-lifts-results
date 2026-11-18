"""Read certified profiles once, without double counting migrated checkpoints."""
from sp_full_profiles import dimensions
from integer_j30_boundary_system import HERE,digest
from projected_checkpoint import atomic,json_bytes
import json


def inventory():
    exact=dimensions((5,5,10,10,10),9,2);records=[];missing=[]
    for degree in (1,3,5):
        profiles=sorted((p,n) for (d,p),n in exact.items() if d==degree)
        for number,(profile,dimension) in enumerate(profiles):
            stem='sp2_j30_first_correction_degree%d_profile%d'%(degree,number)
            for suffix in ('_stratified','_bounded','_compact','_fullnative','_native',''):
                folder=HERE/(stem+suffix);path=folder/'rational_result.json'
                if not path.exists():continue
                record=json.loads(path.read_text())
                assert record['status']=='passed' and record['degree']==degree and record['profile_number']==number
                assert record['rational_basis_dimension']==record['fresh_evaluation_rank']==record['exact_character_dimension']==dimension
                assert record['profile']==json.loads(json.dumps(profile))
                assert digest(folder/'rational_basis.pkl')==record['rational_basis_sha256']
                assert digest(folder/'rational_recipe.json')==record['rational_recipe_sha256']
                assert digest(folder/'rational_evaluation/manifest.json')==record['fresh_evaluation_manifest_sha256']
                records.append(dict(degree=degree,number=number,dimension=dimension,folder=folder.name,
                                    certificate_sha256=digest(path),rational_basis_sha256=record['rational_basis_sha256']))
                break
            else:missing.append(dict(degree=degree,number=number,dimension=dimension))
    return dict(status='complete' if not missing else 'partial',sector=[5,5,10,10,10],letters=9,derivatives=2,
                rational_profile_dimension=sum(r['dimension'] for r in records),previous_rational_degree7_dimension=552,
                source_dimension=30148,records=records,missing=missing,
                scope='Complete rational correction cochains in the listed field-content profiles. No d1 or higher obstruction is claimed.')


if __name__=='__main__':
    result=inventory();atomic(HERE/'correction_profile_inventory.json',json_bytes(result))
    print(json.dumps(result),flush=True)
