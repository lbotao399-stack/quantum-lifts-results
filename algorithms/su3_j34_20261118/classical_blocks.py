"""Resolve the classical maps adjacent to the level-34 obstruction targets by letter content."""
import argparse
import json
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNTIME = HERE.parent / "su3_independent_prime_20261118/runtime/work/sp3-native.nosync"
sys.path.insert(0, str(RUNTIME))
from su3_charge_dimensions_v1 import dimensions, parity_parts
from g2_full_classical_screen_v1 import contents
from sp_full_profiles import profile
from fortuity.full.letters import Q_letter, is_ferm


@lru_cache(None)
def resolve(mu, length, parity):
    multiplicity = Counter()
    rho = (2, 1, 0)
    for spin in (0, 1):
        for order in permutations(range(3)):
            shift = tuple(mu[2+i] + rho[i] - rho[order[i]] for i in range(3))
            weight = (mu[0]+spin, mu[1]-spin) + shift
            if min(weight) < 0:
                continue
            sign = (-1)**(spin + sum(order[i] > order[j] for i in range(3) for j in range(i+1, 3)))
            for content in contents(weight, length):
                counts = Counter(content)
                even, odd = [tuple(sorted((n for l,n in counts.items() if is_ferm(l) == f), reverse=True))
                             for f in (False, True)]
                multiplicity[sum(odd), profile(content)] += sign * parity_parts(even, odd)[parity == -1]
    totals = Counter()
    records = []
    for (degree, pattern), size in sorted(multiplicity.items()):
        assert size >= 0
        if size:
            totals[degree] += size
            records.append(dict(index=len(records), fermions=degree, profile=pattern, dimension=size))
    assert dict(totals) == dimensions(mu, length, parity)
    return records


def adjacent(mu, length, parity):
    source = [r.copy() for r in resolve(mu, length, parity)]
    target = resolve(mu, length+1, parity)
    source_ids = {r["profile"]:r["index"] for r in source}
    target_ids = {r["profile"]:r["index"] for r in target}
    edges = defaultdict(set)
    for content in contents(mu, length):
        i = source_ids.get(profile(content))
        if i is None:
            continue
        for pos, letter in enumerate(content):
            for _, left, right in Q_letter(letter):
                t = target_ids.get(profile(content[:pos]+content[pos+1:]+(left,right)))
                if t is not None:
                    edges[i].add(t)
    for r in source:
        r["targets"] = sorted(edges[r["index"]])
    vertices = {(0,r["index"]) for r in source} | {(1,t["index"]) for t in target}
    graph = defaultdict(set)
    for i, js in edges.items():
        for j in js:
            graph[0,i].add((1,j))
            graph[1,j].add((0,i))
    components = []
    while vertices:
        todo = [min(vertices)]
        visited = set()
        while todo:
            v = todo.pop()
            if v not in visited:
                visited.add(v)
                todo.extend(graph[v]-visited)
        vertices -= visited
        a, b = [sorted(i for side,i in visited if side == chosen) for chosen in (0,1)]
        components.append(dict(sources=a, targets=b, source_dimension=sum(source[i]["dimension"] for i in a),
                               target_dimension=sum(target[i]["dimension"] for i in b)))
    return dict(mu=mu, parity=parity, source_length=length, target_length=length+1,
                sources=source, targets=target, components=components,
                source_dimension=sum(r["dimension"] for r in source),
                target_dimension=sum(r["dimension"] for r in target),
                nonzero_block_entries_upper_bound=sum(r["dimension"]*target[t]["dimension"]
                                                     for r in source for t in r["targets"]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", required=True, type=int, choices=(1,2,3))
    parser.add_argument("--output", required=True, type=Path)
    args=parser.parse_args()
    r=args.order
    mu=(r,r,8-r,5-r,4-r)
    incoming, outgoing = [adjacent(mu, n, 1) for n in (11-r,12-r)]
    assert incoming["targets"] == [{k:v for k,v in row.items() if k!="targets"} for row in outgoing["sources"]]
    if r==1:
        rhs=json.loads((HERE/"first_action.json").read_text())
        profile_terms={tuple((tuple(k),n) for k,n in row["profile"]):row["terms"] for row in rhs["profiles"]}
        targets={t["profile"]:t["index"] for t in incoming["targets"]}
        assert profile_terms.keys() <= targets.keys()
        incoming["first_rhs_trace_terms"]={targets[p]:v for p,v in profile_terms.items()}
        assert sum(incoming["first_rhs_trace_terms"].values())==rhs["formal_trace_terms"]
    record=dict(order=r, incoming=incoming, outgoing=outgoing)
    args.output.write_text(json.dumps(record,indent=2)+"\n")
    for direction in ("incoming","outgoing"):
        b=record[direction]
        print(direction,b["source_dimension"],b["target_dimension"],
              len(b["sources"]),len(b["targets"]),len(b["components"]),b["nonzero_block_entries_upper_bound"])


if __name__=="__main__":
    main()
