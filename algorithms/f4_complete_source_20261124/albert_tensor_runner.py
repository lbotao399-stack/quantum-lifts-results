"""Integer replay of the fixed native Albert boundary and closed-product bank."""
import argparse
import ctypes
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np


class AlbertRunner:
    def __init__(self, inputs, library, prime):
        sys.path.insert(0, str(inputs))
        import albert_boundaries
        import albert_product_cycles
        self.boundaries, self.cycles = albert_boundaries, albert_product_cycles
        self.prime = prime
        self.library = ctypes.CDLL(str(library.resolve()))
        lib = self.library
        lib.f4_albert_create.argtypes = [ctypes.c_uint32] + [ctypes.c_void_p] * 4 + [ctypes.c_size_t]
        lib.f4_albert_create.restype = ctypes.c_void_p
        lib.f4_albert_run.argtypes = [ctypes.c_void_p] * 3
        lib.f4_albert_run.restype = ctypes.c_int
        lib.f4_albert_error.restype = ctypes.c_char_p
        lib.f4_albert_destroy.argtypes = [ctypes.c_void_p]
        lower, upper, exchange, scales = albert_boundaries.tensors()
        assert scales == (28, 252, 84)
        constants = []
        for raw in (lower, upper, exchange):
            integer = np.asarray(raw, dtype=np.int64, order="C")
            if not np.array_equal(raw, integer):
                raise ArithmeticError("nonintegral native Albert tensor")
            constants.append(integer)
        self.columns, self.coefficients = albert_product_cycles.highest()
        active = np.flatnonzero(self.coefficients[0])
        assert len(self.columns) == 96 and len(active) == 87
        self.columns = [self.columns[i] for i in active]
        self.coefficients = self.coefficients[:, active]
        terms = []
        for value, (omitted, gram, pair, last) in zip(self.coefficients[0], self.columns):
            terms.append([int(value), *omitted, *gram, albert_boundaries.PAIRS.index(pair), last])
        terms = np.asarray(terms, dtype=np.int64, order="C")
        self.context = lib.f4_albert_create(prime, *(v.ctypes.data for v in constants), terms.ctypes.data, len(terms))
        if not self.context:
            raise ArithmeticError(lib.f4_albert_error().decode())
        self.constants_sha256 = [sha256(v.astype("<i8").tobytes()).hexdigest() for v in constants]
        self.product_terms_sha256 = sha256(terms.astype("<i8").tobytes()).hexdigest()

    def evaluate(self, point, closed=False):
        scalars, _, marks = self.boundaries.samples(point, self.prime, closed)
        inputs = np.asarray(np.concatenate((scalars, marks)) % self.prime, dtype=np.uint32, order="C")
        output = np.empty(50, dtype=np.uint32)
        if self.library.f4_albert_run(self.context, inputs.ctypes.data, output.ctypes.data):
            raise ArithmeticError(self.library.f4_albert_error().decode())
        return output

    def reference(self, point, closed=False):
        boundary = self.boundaries.evaluate(point, self.prime, closed)
        product = self.cycles.point(point, self.prime, self.columns, self.coefficients, closed)
        return np.concatenate((boundary, product))

    def close(self):
        if self.context:
            self.library.f4_albert_destroy(self.context)
            self.context = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--points", type=int, nargs="+", default=[12, 32, 219])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    runner = AlbertRunner(args.inputs, args.library, args.prime)
    checks = []
    for point in args.points:
        started, cpu = perf_counter(), process_time()
        value = runner.evaluate(point)
        elapsed, cpu = perf_counter() - started, process_time() - cpu
        expected = runner.reference(point)
        if not np.array_equal(value, expected):
            bad = np.flatnonzero(value != expected)
            raise ArithmeticError(str([(int(j), int(value[j]), int(expected[j])) for j in bad]))
        closed = runner.evaluate(point, True)
        if closed.any():
            raise ArithmeticError("classical split does not vanish")
        check = {"point": point, "verified_entries": 50, "classical_split_nonzeros": 0,
                 "integer_seconds": elapsed, "integer_cpu_seconds": cpu,
                 "residues": list(map(int, value)),
                 "residue_sha256": sha256(value.astype("<u2").tobytes()).hexdigest()}
        checks.append(check)
        print(json.dumps(check), flush=True)
    receipt = {"prime": args.prime, "tensor_sha256": runner.constants_sha256,
               "product_terms_sha256": runner.product_terms_sha256,
               "checks": checks,
               "scope": "Complete 49-column boundary bank plus the rationally closed product, at the stated control points."}
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    runner.close()


if __name__ == "__main__":
    main()
