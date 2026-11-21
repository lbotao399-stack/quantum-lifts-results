"""Replay the exceptional degree-eleven scalar determinants and the product exclusion."""
import gzip
import json
from pathlib import Path

from e6_minuscule_443 import compute
from verify import cochains, read
from verify_scalar_stage import check_minor


HERE = Path(__file__).resolve().parent / "scalar_443"
WEIGHT = (4, 4, 3)


def stored_minor(matrix, minor, prime):
    n = minor["rank"]
    assert len(matrix) == n and all(len(row) == n for row in matrix)
    assert all(len(set(minor[key])) == len(minor[key]) == n for key in ("rows", "columns"))
    return check_minor(matrix, dict(minor, rows=list(range(n)), columns=list(range(n))), prime)


def verify():
    groups = ("F4", "E6", "E7", "E8")
    level = {g:cochains(11, g)[1][11] for g in groups}
    sizes = {g:level[g][WEIGHT] for g in groups}
    assert sizes == {"F4":[14], "E6":[14, 23], "E7":[12], "E8":[9]}
    for prime in (1000003, 1000033):
        f4 = read(HERE/f"f4_{prime}.json")["boundaries"]
        assert stored_minor(f4["values"], f4["minor"], prime) == sizes["F4"][0]
        e6 = read(HERE/f"e6_casimir_{prime}.json")["with_casimir"]
        assert stored_minor(e6["values"], e6["minor"], prime) == sizes["E6"][0]
        e6_trace = read(HERE/f"e6_trace_{prime}.json")
        odd, = [r for r in e6_trace["results"] if tuple(r["weight"]) == WEIGHT and r["outer"] == -1]
        assert stored_minor(odd["minor_matrix"], odd["boundary_minor"], prime) == 22
        augmented = read(HERE/f"e6_augmented_{prime}.json")
        assert augmented["prime"] == prime and augmented["seed"] == e6_trace["seed"]
        assert augmented["points"] == e6_trace["points"] and augmented["boundary_minor"] == odd["boundary_minor"]
        minor = augmented["augmented_minor"]
        assert minor["rows"][-1] == odd["parents"]
        assert [augmented["scalar_values"][k] for k in minor["columns"]] == augmented["augmented_matrix"][-1]
        assert stored_minor(augmented["augmented_matrix"], minor, prime) == sizes["E6"][1]
        adjoint = read(HERE/f"adjoint_{prime}.json")
        row, = [r["groups"] for r in adjoint["results"] if tuple(r["weight"]) == WEIGHT]
        for g in ("E7", "E8"):
            assert stored_minor(row[g]["minor_matrix"], row[g]["boundary_minor"], prime) == sizes[g][0]
        extra = json.loads(gzip.decompress((HERE/f"degree11_adj_{prime}.json.gz").read_bytes()))
        complete = {w for w in level["E7"] if sum(w) == 11 and w[2] >= 2}
        assert {tuple(r["weight"]) for r in extra["results"]} == complete - {WEIGHT}
        for r in extra["results"]:
            w = tuple(r["weight"])
            for g, result in r["groups"].items():
                rank = stored_minor(result["minor_matrix"], result["boundary_minor"], prime)
                assert rank + result["cohomology_upper"] == result["dimension"] == level[g][w][0]
                if g in ("E7", "E8") or w == (7, 2, 2):
                    assert rank == result["dimension"]
        extra = json.loads(gzip.decompress((HERE/f"degree11_e6_{prime}.json.gz").read_bytes()))
        for r in extra["results"]:
            w = tuple(r["weight"])
            rank = stored_minor(r["minor_matrix"], r["boundary_minor"], prime)
            assert rank + r["cohomology_upper"] == r["dimension"] == level["E6"][w][r["outer"] == -1]
            if w == (5, 3, 3) and r["outer"] == -1:
                assert rank == 15
    cartan = compute()
    assert json.loads(json.dumps(cartan)) == read(HERE/"e6_minuscule_443.json")
    assert len(cartan["all_odd_subshapes"]) == 16
    assert [w for w in cartan["all_odd_subshapes"] if w[-1]] == [[4, 4, 3, 1]]
    return dict(weight=WEIGHT, scalar_cohomology={"F4":[0], "E6":[0,1], "E7":[0], "E8":[0]},
                E6_Cartan_value=cartan["evaluations"][0]["scalar_443"])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with certificate assertions enabled")
    print(json.dumps(verify(), indent=2))
