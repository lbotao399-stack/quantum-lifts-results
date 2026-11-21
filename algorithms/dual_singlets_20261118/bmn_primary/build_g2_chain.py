"""Build the G2 classical differential against the pinned Taylor trace engine."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


HERE = Path(__file__).resolve().parent


def build(native, sampler, output):
    output.mkdir(parents=True, exist_ok=True)
    pins = json.loads((native/"sha256.json").read_bytes())
    for name in ("core.h", "trace.h", "levels.h", "linalg.h", "batch.h"):
        data = (native/name).read_bytes()
        if hashlib.sha256(data).hexdigest() != pins[name]:
            raise ValueError(f"Native header digest mismatch for {name}")
        if name == "core.h":
            marker = b"  build_weights();"
            if data.count(marker) != 1:
                raise ValueError("The native character initialization has changed")
            data = data.replace(marker, b"  // The G2 driver supplied ADJW and POSR.")
        (output/name).write_bytes(data)
    shutil.copyfile(sampler, output/"bmn_sampler.h")
    shutil.copyfile(HERE/"g2_chain.cpp", output/"g2_chain.cpp")
    binary = output/"g2_chain"
    subprocess.run(["c++", "-std=c++17", "-O3", "-pthread", str(output/"g2_chain.cpp"),
                    "-o", str(binary)], check=True)
    return binary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--native", type=Path, default=HERE.parent/"taylor_hooks/native")
    parser.add_argument("--sampler", type=Path, default=HERE.parent/"b2_level36/bmn_sampler.h")
    args = parser.parse_args()
    build(args.native, args.sampler, args.output)
