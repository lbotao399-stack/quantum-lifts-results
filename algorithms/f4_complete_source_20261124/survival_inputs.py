"""Materialize and authenticate only the explicitly listed public proof inputs."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import quote
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent


def pinned_files():
    manifest = json.loads((HERE / "survival_input_manifest.json").read_text())
    dependencies = json.loads((HERE / "proof_dependencies.json").read_text())
    records = {}
    for record in manifest["files"] + dependencies["inputs"]:
        if not re.fullmatch(r"[0-9a-f]{40}", record["commit"]):
            raise ValueError("the public input commit must be immutable")
        path = Path(record["path"])
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("the public input path must stay below its commit")
        key = record["commit"], record["path"]
        if key in records:
            assert (records[key]["sha256"], records[key]["bytes"]) == (record["sha256"], record["bytes"])
        else:
            records[key] = record
    return manifest, list(records.values())


def verify(root):
    manifest, records = pinned_files()
    for record in records:
        path = root / record["commit"] / record["path"]
        raw = path.read_bytes()
        if len(raw) != record["bytes"] or sha256(raw).hexdigest() != record["sha256"]:
            raise ValueError("public proof input differs from its pin: " + str(path))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--repository", type=Path,
                        help="read the pinned blobs from an existing public Git checkout")
    source.add_argument("--download", action="store_true",
                        help="download only the 100 listed blobs from their immutable public GitHub URLs")
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    manifest, records = pinned_files()
    if args.repository or args.download:
        for record in records:
            path = args.destination / record["commit"] / record["path"]
            if path.is_file():
                existing = path.read_bytes()
                if len(existing) == record["bytes"] and sha256(existing).hexdigest() == record["sha256"]:
                    continue
            if args.repository:
                raw = subprocess.check_output(["git", "show", record["commit"] + ":" + record["path"]], cwd=args.repository)
            else:
                url = ("https://raw.githubusercontent.com/" + manifest["repository"] + "/"
                       + record["commit"] + "/" + quote(record["path"], safe="/"))
                request = Request(url, headers={"User-Agent": "f4-public-proof-replay/1"})
                with urlopen(request, timeout=30) as response:
                    raw = response.read(record["bytes"] + 1)
            if len(raw) != record["bytes"] or sha256(raw).hexdigest() != record["sha256"]:
                raise ValueError("public object differs from its pin: " + record["path"])
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
    manifest = verify(args.destination)
    print(json.dumps({"verified_public_files": len(records),
                      "numerical_input_files": len(manifest["files"]),
                      "additional_theorem_and_provenance_files": len(records) - len(manifest["files"]),
                      "total_bytes": sum(r["bytes"] for r in records)}))


if __name__ == "__main__":
    main()
