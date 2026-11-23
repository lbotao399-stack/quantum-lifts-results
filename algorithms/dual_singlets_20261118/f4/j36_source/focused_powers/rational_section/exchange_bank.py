"""Generate the six by six ordered two-exchange networks from Dot's pinned compiler."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from zipfile import ZipFile

from other_crossing import PACKET_SHA, ROOT, COLOURS

PAIRS = tuple(combinations(range(4), 2))


def scalar_remainder(terms):
    """Alternate the actual ordered scalar insertions in each of the three flavours."""
    result = Counter()
    for encoded, weight in terms["rawResidual"].items():
        insertions = json.loads(encoded)
        for order in permutations(range(3)):
            inversions = sum(x > y for i, x in enumerate(order) for y in order[i + 1 :])
            relabelled = [[order[x] for x in edge] for edge in insertions]
            key = json.dumps(relabelled, separators=(",", ":"))
            result[key] += (-1) ** inversions * weight
    literal = {k: v for k, v in result.items() if v}
    assert literal == terms["altResidual"]
    return literal


def compiler(original):
    m = re.search(r'^const raw=Function\((".*?")\)\(\);', original)
    src = json.loads(m[1])
    old = 'const EX=process.argv[2]==="B"?[[1,2],[0,1]]:[[0,1],[1,2]];'
    assert src.count(old) == 1
    src = src.replace(old, "const EX=JSON.parse(process.argv[2]);")
    old = """const D=[3,4].map(r=>g.chords.map((_,drop)=>{
   const ID=Q.filter(x=>x[2]===r).map(x=>x[0]);
   return det(g.K.map(X=>X.filter((_,j)=>j!==drop)),ID)
  }));"""
    new = """const D=[3,4].map(r=>{
   const ids=Q.filter(x=>x[2]===r).map(x=>x[0]);
   if(new Set(ids).size<ids.length)return g.chords.map(()=>0);
   const sign=sgn(ids), sorted=[...ids].sort((a,b)=>a-b), key=sorted.join(',');
   if(!g.restrictions)g.restrictions=g.chords.map((_,i)=>g.K.map(x=>x.filter((_,j)=>j!==i)));
   if(!g.minors)g.minors=new Map();
   if(!g.minors.has(key))g.minors.set(key,g.restrictions.map(k=>det(k,sorted)));
   return g.minors.get(key).map(v=>sign*v);
  });"""
    assert src.count(old) == 1
    return (
        original[: m.start(1)]
        + json.dumps(src.replace(old, new))
        + original[m.end(1) :]
    )


def build(item):
    pair, compiler_path = item
    # Use the unreduced graph bank. Only literally identical graphs are combined.
    sys.path.insert(0, str(compiler_path.parent))
    from build_bank import combine, dis, overlap

    axes = [PAIRS[i] for i in pair]
    proc = subprocess.run(
        ["node", str(compiler_path), json.dumps(axes)], capture_output=True, check=True
    )
    data = json.loads(proc.stdout)
    banks = {}
    for side, ward in data["W"].items():
        remaining = {k: r for k, v in ward.items() if (r := scalar_remainder(v))}
        gs = combine(list(dis(data["raw"]["disjoint"])) + list(overlap(ward)))
        banks[side] = dict(
            ordinary=remaining,
            banks=[
                {k: g[k] for k in ("E", "d", "cp", "cm", "w", "cubics")} for g in gs
            ],
        )
    return dict(indices=pair, axes=axes, banks=banks)


def generate(archive, output, workers):
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == PACKET_SHA
    with ZipFile(archive) as z, tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        original = z.read(f"{ROOT}/{COLOURS}/typed.js").decode()
        (folder / "generic.js").write_text(compiler(original))
        (folder / "build_bank.py").write_bytes(
            z.read(f"{ROOT}/{COLOURS}/build_bank.py")
        )
        with ProcessPoolExecutor(workers) as pool:
            all_banks = list(
                pool.map(
                    build,
                    [(ab, folder / "generic.js") for ab in product(range(6), repeat=2)],
                )
            )
        for label, i, j in [("A", 0, 3), ("B", 3, 0)]:
            old = json.loads(z.read(f"{ROOT}/{COLOURS}/bank-{label}-reduced.json"))
            trial = all_banks[6 * i + j]["banks"]
            for side in ("left", "right"):
                assert trial[side]["ordinary"] == old[side]["ordinary"] == {}
            if label == "B":
                before = json.loads(z.read(f"{ROOT}/{COLOURS}/bank-B.json"))
                for side in ("left", "right"):
                    assert trial[side]["banks"] == [
                        {k: g[k] for k in ("E", "d", "cp", "cm", "w", "cubics")}
                        for g in before[side]["banks"]
                    ]
        import gzip

        output.write_bytes(
            gzip.compress(
                json.dumps(all_banks, separators=(",", ":")).encode(), mtime=0
            )
        )
        for item in all_banks:
            print(
                item["indices"],
                [
                    (side, len(data["banks"]), len(data["ordinary"]))
                    for side, data in item["banks"].items()
                ],
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    generate(args.packet, args.output, args.workers)
