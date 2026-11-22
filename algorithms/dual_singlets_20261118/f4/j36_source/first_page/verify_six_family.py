"""Check the individual A8 hashes and the characteristic-zero completion of the missing weight."""
from hashlib import sha256
import json
from pathlib import Path

from colour import Colour
from one_fermion import orbits
from six_invariants import run

HERE = Path(__file__).resolve().parent
FIELDS = ("columns", "classical_rank_lower", "boundary_rank_lower", "known_cycle_rank_lower", "unresolved_upper")


def verify():
    data = json.loads((HERE/"six_nonzero_weights.json").read_bytes())
    all_weights = dict(orbits(Colour(),5))
    blocks = data["blocks"]
    seen = {tuple(r["weight"]) for r in blocks}
    assert len(seen) == len(blocks) == 54 and set(all_weights)-seen == {(0,0,0,0)}
    for row in blocks:
        w = tuple(row["weight"])
        assert row["orbit"] == all_weights[w] and row["unresolved_upper"] == 0
        assert row["classical_rank_lower"]+row["known_cycle_rank_lower"] == row["columns"]
        packet = dict(flavour_degrees=[2,2,2], flavour_singlet=True, prime=1000003,
                      blocks=[row], total={f:row["orbit"]*row[f] for f in FIELDS})
        name = "weight-"+",".join(map(str,w))+".json"
        assert sha256((json.dumps(packet,indent=2)+"\n").encode()).hexdigest() == data["sha256"][name]
    totals = {f:sum(r["orbit"]*r[f] for r in blocks) for f in FIELDS}
    proof = run()
    assert json.loads(json.dumps(proof)) == json.loads((HERE/"six_invariants.json").read_bytes())
    assert proof["determinant"] == 292863838464 and totals["columns"] == 45503928
    return dict(nonzero_weights=len(seen), nonzero_weight_totals=totals,
                zero_weight_columns=52*(52*53*51*50//8)-totals["columns"],
                exceptional_invariant_kernel=0)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(verify(),indent=2))
