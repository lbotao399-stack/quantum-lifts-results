"""Check both primitive choices directly on Dot's signed moment tables."""

import argparse
import base64
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import tarfile


ARCHIVES = {
    "so11-j36": "d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b",
    "so12-j42": "fd0f6d814f0d601b25df7b80715e7e67b9d9da674b7bfa8a5e938bc636ca152d",
    "so13-j42": "66e1080b0f98d83688e0545730aa1103fd0ea435e8a8b8c81077f4185c61bfeb",
}


def read_bundle(folder, name):
    encoded = (folder / (name + "-replay.tar.gz.b64")).read_bytes().strip()
    raw = base64.b64decode(encoded, validate=True)
    assert sha256(raw).hexdigest() == ARCHIVES[name]
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        rows = {}
        manifest = archive.extractfile(name + "/MANIFEST.sha256").read().decode()
        for line in manifest.splitlines():
            digest, file_name = line.split("  ", 1)
            file_name = file_name.removeprefix("./")
            data = archive.extractfile(name + "/" + file_name).read()
            assert sha256(data).hexdigest() == digest
            rows[file_name] = data
    return rows


def moments(record, scale):
    row = {key: Fraction(int(value), scale) for key, value in record["gamma"]}
    assert len(row) == len(record["gamma"])
    return row


def check_cancellation(terms, coordinates):
    keys = set().union(*(row for _, row in terms))
    assert len(keys) == coordinates
    for key in keys:
        assert sum(weight * row.get(key, 0) for weight, row in terms) == 0, key


def check_eleven(files, alternate):
    suffix = "-alt-h-native-projected.json" if alternate else "-native-projected.json"
    combination = (("cycle", 1), ("path3", -5), ("path3-alt", Fraction(500, 3)))
    terms = []
    value = Fraction()
    for stem, weight in combination:
        record = json.loads(files[stem + suffix])
        scale = int(record["scale"])
        terms.append((weight, moments(record, scale * 8)))
        value += weight * Fraction(int(record["score"]), scale * 512)
    check_cancellation(terms, 122)
    assert value == 3072000
    return value


def source(files, name):
    record = json.loads(files[name].splitlines()[-1])
    scale = int(record["projectorScale"])
    assert int(record["nativeCombinedScale"]) == 512 * scale
    return moments(record, 8 * scale), Fraction(int(record["zeta"]), 512 * scale)


def check_twelve(files, alternate):
    name = "first-native-full-alternate-h.log" if alternate else "first-native-full.log"
    row, value = source(files, name)
    assert len(row) == 89
    terms = [(1, row)]
    for index, weight in ((14, -374400), (146, 187200)):
        tail = json.loads(files[f"tail-full/{index}.json"])
        terms.append((weight, moments(tail, int(tail["scale"]))))
    check_cancellation(terms, 121)
    assert value == -33177600
    return value


def check_thirteen(files, alternate):
    suffix = "-alt-h-native.log" if alternate else "-batch-native.log"
    path, path_value = source(files, "promoted" + suffix)
    cycle, cycle_value = source(files, "promoted-cycle" + suffix)
    tail = moments(json.loads(files["actual-p7-tail-score.json"]), 1)
    weight = Fraction(5, 2)
    check_cancellation(((weight, path), (1, cycle), (-2400, tail)), 122)
    value = weight * path_value + cycle_value
    assert value == -921600
    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_directory", type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("run with Python assertions enabled")
    for name, check in (("so11-j36", check_eleven), ("so12-j42", check_twelve),
                        ("so13-j42", check_thirteen)):
        files = read_bundle(args.bundle_directory, name)
        first, second = (check(files, alternate) for alternate in (False, True))
        print(f"{name}   {first}   {second}")
