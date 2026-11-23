"""Independent edge-contract and undefined-behavior audit for the new catalogue."""

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

import numpy as np


def representative(word):
    backwards = word[::-1]
    return min(*(word[i:] + word[:i] for i in range(len(word))),
               *(backwards[i:] + backwards[:i] for i in range(len(word))))


def encode(word):
    result = 0
    for letter in word:
        result = result * 4 + letter
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("optimizer", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    sys.path.insert(0, str(args.optimizer))
    import native_catalogue as native

    rng = random.Random(144221)
    tests = []
    for dtype in (np.int32, np.int64, np.uint32, np.uint64):
        raw = np.array([0, 1, (1 << 28) - 1] + [rng.randrange(1 << 28) for _ in range(97)], dtype=dtype)
        for values in (raw, raw.reshape(10, 10), raw[::3], raw[:0], raw[0].reshape(())):
            actual = native.canonical14(values)
            expected = np.array([encode(representative(tuple((int(v) >> (2 * (13 - j))) & 3 for j in range(14))))
                                 for v in values.ravel()], dtype=dtype).reshape(values.shape)
            assert np.array_equal(actual, expected)
            assert actual.shape == values.shape and actual.dtype == values.dtype
            tests.append({"dtype": str(values.dtype), "shape": values.shape, "elements": values.size})
    for bad in (np.array([-1], dtype=np.int64), np.array([1 << 28], dtype=np.uint64), np.array([1.5])):
        try:
            native.canonical14(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("bad canonical input accepted")

    contents = [(2,), (3,), (30,), (31,), (0, 2), (0, 0, 0, 2), (1, 1, 1), (2, 2), (3, 2, 1), (0, 2, 1, 1)]
    counts = []
    for content in contents:
        word = tuple(i for i, count in enumerate(content) for _ in range(count))
        if len(word) > 8:
            expected = [] if len(word) % 2 else [encode(word)]
        else:
            possible = set()
            for w in set(permutations(word)):
                rep = representative(w)
                reverse = w[::-1]
                if len(w) % 2 and min(w[i:] + w[:i] for i in range(len(w))) == min(reverse[i:] + reverse[:i] for i in range(len(w))):
                    continue
                possible.add(encode(rep))
            expected = sorted(possible)
        actual = native.fixed_content_catalogue(content).tolist()
        assert actual == expected
        counts.append({"content": content, "count": len(actual)})
    # Length 31 without a large permutation expansion: independently place the
    # two distinct nonzero letters on a cycle and retain signed bracelets.
    expected = set()
    for i in range(31):
        for j in range(31):
            if i == j:
                continue
            word = [0] * 31
            word[i], word[j] = 1, 2
            expected.add(encode(representative(tuple(word))))
    actual = native.fixed_content_catalogue((29, 1, 1)).tolist()
    assert actual == sorted(expected) and len(actual) == 15
    counts.append({"content": [29, 1, 1], "count": 15})

    source = (args.optimizer / "dihedral_catalogue.cpp").resolve()
    harness = '#include "' + str(source) + '"\n' + r'''
#include <cassert>
#include <climits>
int main() {
    std::size_t size = 123;
    const int invalid[][4] = {{INT_MAX,INT_MAX,0,0},{-1,4,0,0},{31,1,0,0},{0,0,0,0}};
    for (auto& c: invalid) { auto* p=dihedral_catalogue(4,c,&size); assert(!p && size==0); }
    const int edge[][4] = {{31,0,0,0},{30,0,0,0},{29,1,1,0},{0,0,0,2},{1,1,1,0}};
    for (auto& c: edge) { auto* p=dihedral_catalogue(4,c,&size); assert(p); dihedral_free(p); }
    std::uint32_t x = 8127331;
    for (int i=0;i<10000;++i) {
        x ^= x<<13; x ^= x>>17; x ^= x<<5;
        auto input=x&((1u<<28)-1); std::uint32_t output=0;
        dihedral_canonical14(&input,&output,1);
        assert(output<=input && output<(1u<<28));
    }
}
'''
    with tempfile.TemporaryDirectory(prefix="quantum-native-audit-") as tmp:
        code = Path(tmp) / "audit.cpp"
        binary = Path(tmp) / "audit"
        code.write_text(harness)
        subprocess.run(["c++", "-std=c++17", "-O1", "-fsanitize=undefined", "-fno-sanitize-recover=all", str(code), "-o", str(binary)], check=True)
        completed = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
        assert not completed.stderr
    result = {
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "python_wrapper_sha256": sha256((args.optimizer / "native_catalogue.py").read_bytes()).hexdigest(),
        "canonical_shape_dtype_cases": tests,
        "fixed_content_edge_cases": counts,
        "invalid_python_input_cases": 3,
        "ubsan_invalid_native_content_cases": 4,
        "ubsan_valid_content_cases": 5,
        "ubsan_canonical_cases": 10000,
        "undefined_behavior_sanitizer_passed": True,
        "scope": "edge contracts and bit arithmetic; full degree14 coordinate equality audited separately",
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"canonical_cases": len(tests), "content_cases": len(counts), "ubsan": "passed"}))


if __name__ == "__main__":
    main()
