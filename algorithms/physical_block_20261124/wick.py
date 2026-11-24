"""Physical free radial Wick form on exact complete native trace spaces.

Colour generators are Hermitian and tr(T_A T_B)=delta_AB. Use precisely the
Chang--Feng--Lin--Tao physical superfield norm (arXiv:2306.04673v2, (3.13),
App. E), for native Taylor coefficients rescaled by (4*pi/g)*2**(a+b).
The common rescaling gives Qphysical = (g/(4*pi))*Qnative at leading order.
This is a fixed-order Wick calculation, not higher-order renormalized matching.
"""
from native import *
import ctypes, subprocess, numpy as np, time
from fractions import Fraction
from flint import fmpq_mat, fmpq

def library():
    target=HERE/'wick_native.so';source=HERE/'wick.cpp'
    if not target.exists() or target.stat().st_mtime<source.stat().st_mtime:
        subprocess.run(['c++','-O3','-std=c++17','-shared','-fPIC',str(source),'-o',str(target)],check=True)
    lib=ctypes.CDLL(str(target))
    lib.wick_gram.argtypes=[ctypes.c_int]*3+[ctypes.c_void_p]*5
    lib.wick_gram.restype=ctypes.c_int
    lib.wick_cache_size.restype=ctypes.c_uint64
    return lib

LIB=library()

def letter_norm(letter):
    a,b,s=letter;n=a+b;k=len(s)
    assert n+k>=1
    return Fraction(math.factorial(n+k-1),math.factorial(a)*math.factorial(b))

def profile(word):return tuple(sorted(l for tr in word for l in tr))

def colour_gram(words,size):
    if not words:return np.zeros((0,0),dtype=np.int64)
    n=sum(map(len,words[0]));alphabet=sorted(set(profile(words[0])));ids={l:i for i,l in enumerate(alphabet)}
    assert all(profile(w)==profile(words[0]) for w in words)
    letters=[];successors=[];maps=[]
    for w in words:
        letters.append([ids[l] for tr in w for l in tr])
        successor=[];start=0;blocks=[]
        for tr in w:
            successor+=list(range(start+1,start+len(tr)))+[start]
            blocks.append(list(range(start,start+len(tr))));start+=len(tr)
        remap=[0]*n
        for j,i in enumerate(i for block in sorted(blocks,key=lambda b:(len(b),b)) for i in block):remap[i]=j
        successors.append(successor);maps.append(remap)
    arrays=[np.ascontiguousarray(a,dtype=np.int32) for a in (letters,successors,maps,[int(is_ferm(l)) for l in alphabet])]
    out=np.empty((len(words),len(words)),dtype=np.int64)
    result=LIB.wick_gram(n,size,len(words),*(a.ctypes.data for a in arrays),out.ctypes.data)
    assert result==0,result
    return out

def gram(rows,size=4,label=''):
    """Entire Gram matrix, with every raw word retained before exact relations."""
    start=time.monotonic();by_profile=defaultdict(set)
    for row in rows:
        for w in row:by_profile[profile(w)].add(w)
    result=fmpq_mat(len(rows),len(rows));receipts=[]
    for step,(content,wordset) in enumerate(sorted(by_profile.items())):
        words=sorted(wordset);wid={w:i for i,w in enumerate(words)}
        coefficients=defaultdict(dict)
        for i,row in enumerate(rows):
            r={wid[w]:c for w,c in row.items() if w in wordset}
            if r:coefficients[i]=r
        raw=colour_gram(words,size)
        factor=Fraction(1,2**len(content))
        for letter in content:factor*=letter_norm(letter)
        indices=sorted(coefficients)
        a=fmpz_mat([[coefficients[i].get(j,0) for j in range(len(words))] for i in indices])
        b=a*fmpz_mat(raw.tolist())*a.transpose()
        for i,ri in enumerate(indices):
            for j,rj in enumerate(indices):
                if b[i,j]:result[ri,rj]+=fmpq(int(b[i,j])*factor.numerator,factor.denominator)
        receipts.append({'profile':[[a,b,list(s)] for a,b,s in content],
                         'words':len(words),'rows':len(indices),'factor':str(factor),
                         'raw_sha256':hashlib.sha256(raw.astype('<i8').tobytes()).hexdigest()})
        if (step+1)%20==0 or step+1==len(by_profile):
            print(label,'profiles',step+1,'/',len(by_profile),'cache',LIB.wick_cache_size(),'seconds',time.monotonic()-start,flush=True)
    assert result==result.transpose()
    return result,receipts

def export_matrix(matrix):return [[str(matrix[i,j]) for j in range(matrix.ncols())] for i in range(matrix.nrows())]
def import_matrix(value):return fmpq_mat(value)
