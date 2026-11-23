#!/usr/bin/env python3
"""Materialize immutable upstream inputs without changing a Git checkout."""
from pathlib import Path
import argparse,hashlib,json,subprocess

HERE=Path(__file__).resolve().parent
DUAL='9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9'
RUNTIME='082be7070d207e2cdecfd0a91b03d5e014fd303d'

def git(repo,*args,input=None):
    return subprocess.run(['git','-C',str(repo),*args],input=input,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE/'inputs')
    args=parser.parse_args();repo=args.repo.resolve();out=args.output.resolve()
    requests=[]
    for commit,prefix in [(DUAL,'algorithms/dual_singlets_20261118'),(RUNTIME,'algorithms/su3_independent_prime_20261118/runtime')]:
        try:git(repo,'cat-file','-e',commit+'^{commit}')
        except subprocess.CalledProcessError:git(repo,'fetch','--depth','1','--filter=blob:none','origin',commit)
        for line in git(repo,'ls-tree','-r',commit,prefix).decode().splitlines():
            header,name=line.split('\t',1);mode,kind,sha=header.split()
            if kind!='blob':continue
            if commit==DUAL and len(Path(name).parts)!=3 and '/sp3_omega/' not in name:continue
            requests.append((commit,name,sha))
    git(repo,'fetch','--no-tags','--no-write-fetch-head','origin','--stdin',input=('\n'.join(sorted({r[2] for r in requests}))+'\n').encode())
    manifest=[]
    for commit,name,sha in requests:
        data=git(repo,'cat-file','blob',sha);path=out/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        manifest.append({'commit':commit,'path':name,'git_blob':sha,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    (out/'INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'files':len(manifest),'bytes':sum(r['bytes'] for r in manifest),'output':str(out)}))

if __name__=='__main__':main()
