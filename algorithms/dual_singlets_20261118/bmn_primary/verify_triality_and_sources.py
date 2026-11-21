"""Verify the SO(8) triality indices and the G2 source exclusions over Q."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from bmn_primary_index import cochains
from torus_index import ROOTS, prime_moduli
from verify import receipt


def triality():
    receipt("7658707")
    folder = HERE/"7658707"
    components = []
    for name in ("D4","D4_2","D4_3"):
        d = json.loads((folder/(name+".json")).read_bytes())
        assert (d["group"],d["maximum_level"],d["maximum_first_part"],d["index_type"]) == (name,144,28,"primary")
        _,coordinate,shift = ROOTS[name]
        assert d["grid"] > 72*coordinate+shift
        primes = list(prime_moduli(d["grid"],int(d["coefficient_bound"])))
        from math import prod
        assert prod(primes) == int(d["modulus"]) > 2*int(d["coefficient_bound"])
        assert set(primes)==set(map(int,d["modular_sha256"]))
        components.append({tuple(r[:3]):r[3] for r in d["rows"]})
    expected = {(a,b,c) for a in range(2,29) for b in range(2,a+1)
                for c in range(2,b+1) if a+b+c<=72}
    assert all(set(t)==expected for t in components)
    previous = json.loads((HERE/"7658467/D4.json").read_bytes())
    for a,b,c,n in previous["rows"]:
        if (a,b,c) in expected:
            assert components[0][a,b,c] == n
    ranks = {}
    for w in expected:
        identity,swap,cycle = (table[w] for table in components)
        values = (identity+3*swap+2*cycle,identity-3*swap+2*cycle,identity-cycle)
        denominators = (6,6,3)
        assert all(v%d==0 for v,d in zip(values,denominators))
        ranks[w] = tuple(v//d for v,d in zip(values,denominators))
    assert len(ranks)==3571 and ranks[3,3,3] == (0,0,1)
    assert all(min(v)>=0 for v in ranks.values())
    return len(ranks)


def odd_sources():
    dims = cochains(18,"G2")[1]
    weights = ((6,5,4),(6,6,4),(6,6,6))
    open_lengths = {w:{length for length,table in dims.items() if length%2 and w in table} for w in weights}
    for job in ("5691136","5691142","5691420","5691422","5691423","5691483"):
        receipt(job)
        for path in (HERE/job).glob("chain_6_*.json"):
            d=json.loads(path.read_bytes())
            w=tuple(d["weight"])
            for r in d["blocks"]:
                length=r["length"]
                assert r["incoming"]+r["outgoing"] == r["dimension"] == dims[length][w][0]
                open_lengths[w].remove(length)
    assert open_lengths == {(6,5,4):set(), (6,6,4):{11,13}, (6,6,6):{11,13,15}}
    return {str(w):sorted(v) for w,v in open_lengths.items()}


if __name__ == "__main__":
    print("SO(8) triality partitions",triality())
    print("G2 outstanding odd lengths",odd_sources())
