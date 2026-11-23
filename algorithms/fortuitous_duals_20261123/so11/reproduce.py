#!/usr/bin/env python3
"""Rebuild the complete scalar SO11 relation quotient from pinned sources.

Requires git, Node.js, C++17, Python>=3.12 and python-flint0.9.
Never changes the supplied checkout or writes to GitHub.
"""
from pathlib import Path
import argparse,base64,hashlib,io,json,os,shutil,subprocess,sys,tarfile,time

HERE=Path(__file__).resolve().parent
DOT='e1d184d388c588c0a5a825a19ecbca4e24427242'
CHAR='9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9'
REVIEW='ac8bdb18447d32cbdd0a67806a87fdd3d9c1df1d'

def git(repo,*args,**kw):
    return subprocess.run(['git','-C',str(repo),*args],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,**kw).stdout

def fetch(repo,commit,paths,destination):
    try:git(repo,'cat-file','-e',commit+'^{commit}')
    except subprocess.CalledProcessError:git(repo,'fetch','--depth','1','--filter=blob:none','origin',commit)
    objects=[]
    for path in paths:
        lines=git(repo,'ls-tree','-r',commit,'--',path).decode().splitlines();assert len(lines)==1
        header,name=lines[0].split('\t');mode,kind,sha=header.split();assert kind=='blob' and name==path
        objects.append((name,sha))
    git(repo,'fetch','--no-tags','--no-write-fetch-head','origin','--stdin',input=('\n'.join(sha for _,sha in objects)+'\n').encode())
    receipts=[]
    for path,sha in objects:
        b=git(repo,'cat-file','blob',sha);p=destination/Path(path).name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
        receipts.append({'commit':commit,'path':path,'git_blob':sha,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
    return receipts

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--work',type=Path,required=True);a=ap.parse_args()
    start=time.monotonic();repo=a.repo.resolve();work=a.work.resolve();assert work!=HERE
    work.mkdir(parents=True,exist_ok=True)
    for p in HERE.iterdir():
        if p.is_file() and p.suffix in ['.py','.js','.json']:shutil.copy2(p,work/p.name)
    inputs=work/'inputs';inputs.mkdir(exist_ok=True)
    manifest=fetch(repo,DOT,['reviews/dot-orthogonal-d3/so11-j36-replay.tar.gz.b64'],inputs)
    manifest+=fetch(repo,CHAR,['algorithms/dual_singlets_20261118/b2_schur_complex.py','algorithms/dual_singlets_20261118/orthogonal_glueball_sources.py','algorithms/dual_singlets_20261118/orthogonal_glueball_first_pair.py'],inputs/'characters')
    manifest+=fetch(repo,REVIEW,['algorithms/so11_source_comparison_20261123/so11_derivative_review.py'],inputs)
    archive=base64.b64decode((inputs/'so11-j36-replay.tar.gz.b64').read_bytes())
    assert hashlib.sha256(archive).hexdigest()=='d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b'
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz') as tf:tf.extractall(inputs,filter='data')
    dot=inputs/'so11-j36'
    for line in (dot/'MANIFEST.sha256').read_text().splitlines():
        digest,name=line.split(maxsplit=1);assert hashlib.sha256((dot/name.lstrip('*')).read_bytes()).hexdigest()==digest
    native=work/'probe_trials';native.mkdir(exist_ok=True)
    original=(dot/'bmngram.cpp').read_text()
    for n,name in [(15,'bmngram'),(14,'bmngram14')]:
        code=original if n==15 else original.replace('constexpr int n=15,m=12;','constexpr int n=14,m=12;')
        cpp=native/(name+'.cpp');cpp.write_text(code)
        subprocess.run(['c++','-O3','-std=c++17',str(cpp),'-o',str(native/name)],check=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',SO11_REVIEW_PATH=str(inputs/'so11_derivative_review.py'),SO11_PFAFFIAN_PATH=str(inputs/'characters/orthogonal_glueball_first_pair.py'))
    os.environ.update(env);sys.path.insert(0,str(work))
    from probe_q0_relations import source,q_inputs
    reference=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    selected=reference['source_relation_basis_indices'];qselected=reference['terminal_basis_Q0_source_indices']
    recorded_specs={r['index']:r['source_spec'] for r in reference['input_receipts']}
    assert all(source(i)==recorded_specs[i] for i in selected)
    src=work/'source_factor5';src.mkdir(exist_ok=True);image=work/'factor5_rank';image.mkdir(exist_ok=True)
    history=[{'index':i,'source':source(i),'independent':i in selected,'image_independent':i in qselected} for i in range(max(selected)+1)]
    for step,i in enumerate(selected,1):
        spec=source(i)
        inp='5 0\n'+' '.join(str('ABCx'.index(c)) for c in spec['species'])+'\n'+' '.join(['0']*13+['1'])+'\n'+' '.join(map(str,spec['endpoints']))+'\n'
        raw=src/f'input_{i}.raw';raw.write_text(subprocess.check_output([str(native/'bmngram14')],input=inp,text=True))
        result=json.loads(subprocess.check_output(['node',str(work/'project_gram_source.js'),str(dot),str(raw)],text=True))
        (src/f'projected_{i}.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
        pair=[]
        for j,(sign,data) in enumerate(q_inputs(spec,5)):
            file=image/f'input_{i}_{j}.raw';file.write_text(subprocess.check_output([str(native/'bmngram')],input=data,text=True));pair.append(file)
        separately=json.loads(subprocess.check_output(['node',str(work/'project_gram_pair.js'),str(dot),*map(str,pair)],text=True))
        assert separately['projected']==result['Q0_projected']
        print(f'Fresh predecessor {step}/21: source {i}, exact Q0/projector comparison passed',flush=True)
    (src/'rank.json').write_text(json.dumps({'rows':history},indent=2)+'\n')
    target=work/'terminal_factor5';target.mkdir(exist_ok=True);(target/'zero.raw').write_text('! 0\n')
    for label,script in [('cycle','seed_cycle.py'),('path3','seed_3path.py'),('path3-alt','seed_3path_alt.py')]:
        data=subprocess.check_output([sys.executable,str(dot/script)],text=True).splitlines();data[0]='5 0'
        raw=target/(label+'.raw');raw.write_text(subprocess.check_output([str(native/'bmngram')],input='\n'.join(data)+'\n',text=True))
        result=subprocess.check_output(['node',str(work/'project_gram_pair.js'),str(dot),str(raw),str(target/'zero.raw')],text=True)
        (target/(label+'.json')).write_text(result)
    commands=[['check_relation_characters.py','--character-source',str(inputs/'characters')],['certify_relation_quotient.py'],['verify_scalar_certificate.py'],['solve_scalar_preimages.py'],['verify_q0_relation_probes.py'],['codex_scalar_return.py'],['solve_codex_preimage.py'],['verify_scalar_preimages.py']]
    for command in commands:
        log=subprocess.check_output([sys.executable,str(work/command[0]),*command[1:]],text=True,env=env,cwd=work)
        (work/(Path(command[0]).stem+'.log')).write_text(log);print(command[0]+': passed',flush=True)
    actual=json.loads((work/'complete_scalar_relation_quotient.json').read_text())
    for key in ['source_basis_probe_matrix','terminal_basis_probe_matrix','Q0_relation_matrix_on_probes','dual_coefficients','dot_seed_pairings','dot_combination_pairing']:
        assert actual[key]==reference[key],key
    receipt={'passed':True,'fresh_predecessor_circuits':21,'fresh_terminal_circuits':3,'all21_source_Q0_projection_comparisons':True,
             'source_relation_dimension':21,'terminal_relation_dimension':17,'rational_Q0_rank':16,'fortuitous_quotient_dimension':1,
             'dual_support':6,'Dot_scalar_combination_pairing':'1/6','Codex_scalar_ancestor_pairing':'-189/128',
             'combined_scalar_preimage_found':True,'combined_scalar_predecessor_terms':16,
             'input_files':manifest,'seconds':time.monotonic()-start}
    (work/'FULL_REPLAY_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
