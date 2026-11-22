"""Test the SU(4) weights and the free endpoint-contraction bounds of the R-007 Ward channel."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import zipfile

from verify import HASHES, blob, endpoint

FLAVOURS = frozenset(range(4))
MATTER = frozenset(range(3))


def occupations(fields):
    return tuple(sum(i in field for field in fields) for i in range(4))


def reverse(fields):
    return [FLAVOURS - field for field in fields]


def source_profile():
    compressed = blob("dot-r007-exact-input/phase34-public-q0-q1-replay.zip")
    with zipfile.ZipFile(io.BytesIO(compressed)) as archive:
        data = archive.read(
            next(f for f in archive.namelist() if f.endswith("/dmk_o330_hw.json"))
        )
    assert (
        sha256(data).hexdigest()
        == "8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34"
    )
    terms = json.loads(data)["terms"]
    counts = Counter()
    for coefficient, word in terms:
        assert Fraction(coefficient)
        letters = [letter for tr in word for letter in tr]
        assert len(letters) == 10 and sum(x < 3 for x in letters) == 5
        fields = [MATTER - {x} if x < 3 else {x - 3} for x in letters]
        counts[occupations(fields), occupations(reverse(fields))] += 1
    assert counts == {((3, 6, 6, 0), (7, 4, 4, 10)): 6139}
    return next(iter(counts))


def radial_profile():
    inputs = {f: blob("dot-r007-t2-primary/" + f) for f in HASHES}
    assert all(sha256(data).hexdigest() == HASHES[f] for f, data in inputs.items())
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "replay.py"
        path.write_bytes(inputs["replay.py"])
        engine = SimpleNamespace(**runpy.run_path(str(path)))
        target = endpoint(engine, inputs)
    counts = Counter()
    for word, coefficient in target.items():
        assert coefficient
        fields = [set(letter[2]) | {3} for trace in word for letter in trace]
        counts[occupations(fields), occupations(reverse(fields))] += 1
    assert counts == {((5, 2, 2, 9), (4, 7, 7, 0)): 5080}
    return next(iter(counts))


def admissible(weights, modes_per_flavour):
    # Equal Cartan charge allows only a common integer shift of all four occupations.
    return list(range(-min(weights), modes_per_flavour - max(weights) + 1))


def main():
    source = source_profile()
    radial = radial_profile()
    current = ((1, 1, 1, 0), (0, 0, 0, 1))
    results = []
    for bare, target, insertion in zip(source, radial, current, strict=True):
        assert not admissible(target, 6)
        total = [a + b + c for a, b, c in zip(bare, target, insertion, strict=True)]
        assert len(set(total)) == 1
        results.append(
            dict(
                source=bare,
                radial_target=target,
                current=insertion,
                endpoint_current_minimum=total,
                minimum_direct_free_pairs=max(total) - 6,
                radial_allowed_common_shifts=admissible(target, 6),
            )
        )
    assert [r["minimum_direct_free_pairs"] for r in results] == [3, 5]
    # At one point the leading chiral gluino and scalar use two local exact
    # spinor combinations and 2(N-2) lifted variables of each flavour.
    assert max(source[0]) == 6 and max(source[0]) > 2 + 2 * (3 - 2)
    assert max(source[1]) > 6
    return dict(
        group="SU(3)",
        instanton_number_absolute_value=1,
        modes_per_flavour=6,
        leading_local_modes_per_flavour=4,
        source_words=6139,
        radial_words=5080,
        orientations=results,
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    print(json.dumps(main(), indent=2))
