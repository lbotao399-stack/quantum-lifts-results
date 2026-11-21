"""Factor complete SU(3) first images by a shared trace and expand in the scalar chart."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess

from su3_subregular_image import ARCHIVES, ODD, THETA_DEGREE, entries, write


def letter_code(letter):
    z, w, theta = letter
    return (2*z+w)*8+sum(1<<i for i in theta)


def parity(trace):
    return sum((1-(k&1))+THETA_DEGREE[(k>>1)%4] for k in trace)%2


def prepare(level, path):
    name, n, denominator, expected = ARCHIVES[level]
    archive = Path(__file__).resolve().parents[1]/name
    groups = defaultdict(Counter)
    alphabet = set()
    digest = hashlib.sha256()
    count = 0
    with gzip.open(archive, "rb") as stream:
        for line in stream:
            digest.update(line)
            c, raw = json.loads(line)
            traces = [bytes(map(letter_code, t)) for t in raw]
            letters = b"".join(traces)
            assert len(letters) == 11 and sum(k&1 for k in letters) == 7
            alphabet.update(letters)
            pos = max(range(len(traces)), key=lambda i:(len(traces[i]), traces[i]))
            longest = traces[pos]
            sign = (-1)**(parity(longest)*sum(map(parity, traces[:pos])))
            rest = tuple(traces[:pos]+traces[pos+1:])
            numerator = Fraction(c)*denominator*sign
            assert numerator.denominator == 1
            groups[longest][rest] += numerator.numerator
            count += 1
    if count != n or digest.hexdigest() != expected:
        raise ValueError("Incomplete or unexpected source archive")
    with path.open("w") as stream:
        print(*map(int,ODD), file=stream)
        print(len(alphabet), file=stream)
        for k in sorted(alphabet):
            print(k, file=stream)
            for poly in entries(k):
                print(len(poly), file=stream)
                for word,c in sorted(poly.items()):
                    print(c, len(word), *word, file=stream)
        print(len(groups), file=stream)
        for word,sub in sorted(groups.items()):
            print(word.hex(), len(sub), file=stream)
            for rest,c in sorted(sub.items()):
                print(c, len(rest), *(t.hex() for t in rest), file=stream)
    return dict(level=level, full_image_words=n, full_image_sha256=expected, denominator=denominator,
                sign=1 if level==32 else -1, trace_groups=len(groups),
                grouped_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def collect(files):
    result = Counter()
    for path in files:
        with path.open("rb") as stream:
            while chunk := stream.read(24*60000):
                for packed,n in struct.iter_unpack("<16sq",chunk):
                    size = packed[15]
                    assert 0 < size <= 15 and packed[size:15] == bytes(15-size)
                    word = packed[:size]
                    assert word == bytes(sorted(word)) and n
                    result[word] += n
    return [[w.hex(),c] for w,c in sorted(result.items()) if c]


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", required=True, type=int, choices=(32,34))
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--degree", type=int, choices=range(1,16), default=15)
    p.add_argument("--roots", type=int, choices=range(5), default=4)
    args = p.parse_args()
    if not __debug__ or args.workers<1:
        p.error("Enable assertions and select a positive worker count")
    work = args.work.resolve()
    work.mkdir(exist_ok=True, parents=True)
    inputs = work/f"level-{args.level}.txt"
    prepared = work/f"level-{args.level}.json"
    if prepared.exists():
        metadata = json.loads(prepared.read_bytes())
        assert metadata["level"]==args.level and metadata["full_image_sha256"]==ARCHIVES[args.level][3]
        assert metadata["grouped_sha256"]==hashlib.sha256(inputs.read_bytes()).hexdigest()
    else:
        metadata = prepare(args.level, inputs)
        write(prepared, metadata)
    exe = work/"expand"
    subprocess.run(["c++","-O3","-std=c++17",str(Path(__file__).with_name("subregular_trace_expansion.cpp")),
                    "-o",str(exe)],check=True)
    groups = metadata["trace_groups"]
    chunk = (groups+args.workers-1)//args.workers
    def expand_range(start):
        end=min(start+chunk, groups)
        dest=work/f"level-{args.level}-{args.degree}-{args.roots}-{start}-{end}.bin"
        subprocess.run([str(exe),str(inputs),str(dest),*map(str,(args.degree,args.roots,start,end))],check=True)
        return dest
    with ThreadPoolExecutor(args.workers) as pool:
        files=list(pool.map(expand_range, range(0,groups,chunk)))
    rows=collect(files)
    degrees=Counter(len(w)//2 for w,_ in rows)
    data=dict(**metadata,degree=args.degree,root_limit=args.roots,
              nonzero_by_degree={str(i):degrees[i] for i in range(args.degree+1)},coefficients=rows)
    if args.roots==1:
        data["linear_roots"]=True
    write(args.out,data)
    print(args.level,args.degree,args.roots,len(rows),data["nonzero_by_degree"],flush=True)
