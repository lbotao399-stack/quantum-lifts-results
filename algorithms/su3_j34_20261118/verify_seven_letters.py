"""Test the seven-entry polynomial keys against direct integer matrix/exterior products.
The exterior determinant below is the signed coefficient of the top ordered wedge.
Only the values in the selected random assignments are point checks.
"""
import argparse
from functools import lru_cache
from itertools import permutations
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib, subprocess
from math import prod
import finite_trace_identity as ft
import replay_formal_blocks as s

# The low-content side of the 70 equations; traces are evaluated in the supplied order.
WORDS=[('0000000001040505060810',),('000005','0001040506081000'),('00000000','01040505060810'),('000005','000005','0104060810'),('000105','000405','0000060810')]
def read_polys(exe,odd):
    with TemporaryDirectory() as d:
        d=Path(d);inp=d/'in';out=d/'out';inp.write_text(str(len(WORDS))+'\n'+''.join(f"{len(r)} {' '.join(r)}\n" for r in WORDS));subprocess.run([str(exe.resolve()),str(inp),odd,str(out)],check=True);raw=out.read_bytes();lines=raw.splitlines();assert len(lines)==len(WORDS)
    data=[]
    for line in lines:
        n,*x=map(int,line.split());row=dict(zip(x[::2],x[1::2]));assert len(x)==2*n and len(row)==n and all(row.values());data.append(row)
    return data,hashlib.sha256(raw).hexdigest()
def matrix_value(row,M,content):
    k=content.count(0);n=sum(ft.ODD[i] for i in content);perm=[(r,ft.sign(r)) for r in permutations(range(n))];x,y=M[0,-1][0][0],M[0,-1][1][1];result=0
    @lru_cache(None)
    def wedge(ids):
        return sum(sg*prod(M[l,h][a][b] for (l,a,b),h in zip(ids,r)) for r,sg in perm)
    for key,c in row.items():
        raw=key>>3;e=key&7;letters=[]
        while raw:
            v=raw&255;assert v>=8;letters.append(v);raw>>=8
        letters.reverse();assert len(letters)==7 and letters==sorted(letters) and e<=k
        ids=[(v//8,(v%8)//3,(v%8)%3) for v in letters];assert bytes(sorted([0]*k+[l for l,i,j in ids]))==content
        even=prod(M[l,-1][a][b] for l,a,b in ids if not ft.ODD[l]);o=tuple(t for t in ids if ft.ODD[t[0]]);result+=c*x**e*y**(k-e)*even*wedge(o)
    return result
def run(exe):
    odd=bytes(i for i,v in enumerate(ft.ODD) if v).hex();data,digest=read_polys(exe,odd);controls=[];target=bytes.fromhex('0000000001040505060810')
    for word,row in zip(WORDS,data):
        w=tuple(map(bytes.fromhex,word));assert bytes(sorted(b''.join(w)))==target
        vals=[]
        for seed in (3407,7164):
            M,n=ft.assignment(w,3,seed);x,y=(2,-3) if seed==3407 else (3,5);M[0,-1]=((x,0,0),(0,y,0),(0,0,-x-y));a=ft.exterior_value(w,M,n);b=matrix_value(row,M,target);assert a==b,(word,seed,a,b);vals.append(a)
        controls.append(dict(word=word,monomials=len(row),values=vals))
    assert not data[3] and all(all(v!=0 for v in z['values']) for j,z in enumerate(controls) if j!=3)
    return dict(content=target.hex(),ordered_column_sha256=digest,controls=controls,seeds=[3407,7164],diagonal_eigenvalues=[[2,-3,1],[3,5,-8]],method='Independent integer trace evaluation with full exterior permutations versus packed polynomial keys; fourth word squares the same odd trace.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--exe',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__;res=run(a.exe);s.write(a.out,res);print([(u['monomials'],u['values']) for u in res['controls']],flush=True)
