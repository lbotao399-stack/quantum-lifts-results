"""Exact C-odd profile decomposition for the nonabelian second correction."""
from collections import defaultdict
from pathlib import Path
import hashlib
import json
import time

from su3_correction_profiles_v1 import profile_dimensions, contents, profile, Q_letter


def main():
    start = time.monotonic()
    mu = (2, 2, 5, 2, 2)
    source = profile_dimensions(mu, 8)
    target = profile_dimensions(mu, 9)
    assert sum(source.values()) == 49328
    assert sum(target.values()) == 63983
    source_by_profile = {p: (i, degree, dimension)
                         for i, ((degree, p), dimension) in enumerate(sorted(source.items()))}
    target_by_profile = {p: (i, degree, dimension)
                         for i, ((degree, p), dimension) in enumerate(sorted(target.items()))}
    edges = defaultdict(set)
    for content in contents(mu, 8):
        current = profile(content)
        if current not in source_by_profile:
            continue
        index = source_by_profile[current][0]
        for position, letter in enumerate(content):
            for _, left, right in Q_letter(letter):
                following = profile(content[:position] + content[position+1:] + (left, right))
                if following in target_by_profile:
                    edges[index].add(target_by_profile[following][0])
    graph = defaultdict(set)
    for index, targets in edges.items():
        for j in targets:
            graph['s%d' % index].add('t%d' % j)
            graph['t%d' % j].add('s%d' % index)
    # Include isolated profiles as well as the coupled components.
    remaining = {'s%d' % i for i in range(len(source))} | {'t%d' % i for i in range(len(target))}
    components = []
    while remaining:
        todo, visited = [min(remaining)], set()
        while todo:
            vertex = todo.pop()
            if vertex in visited:
                continue
            visited.add(vertex)
            todo.extend(graph[vertex] - visited)
        remaining -= visited
        ss = sorted(int(v[1:]) for v in visited if v.startswith('s'))
        tt = sorted(int(v[1:]) for v in visited if v.startswith('t'))
        components.append(dict(source_profiles=ss, target_profiles=tt,
            source_dimension=sum(n for i, _, n in source_by_profile.values() if i in ss),
            target_dimension=sum(n for i, _, n in target_by_profile.values() if i in tt)))
    result = dict(status='exact_integer_second_correction_profiles_and_structural_couplings',
        mu=mu, parity=-1, source_length=8, target_length=9,
        source_dimension=sum(source.values()), target_dimension=sum(target.values()),
        sources=[dict(index=i, fermions=d, dimension=n, profile=p, targets=sorted(edges[i]))
                 for p, (i, d, n) in sorted(source_by_profile.items(), key=lambda item: item[1])],
        targets=[dict(index=i, fermions=d, dimension=n, profile=p)
                 for p, (i, d, n) in sorted(target_by_profile.items(), key=lambda item: item[1])],
        connected_components=components, seconds=time.monotonic()-start,
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Exact dimensions and structural Q0 couplings only. Does not solve the correction or quotient by the first-order homogeneous freedom.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('status','source_dimension','target_dimension','connected_components','seconds')}))
    print(json.dumps(dict(source_profiles=len(source), target_profiles=len(target))))


if __name__ == '__main__':
    main()
