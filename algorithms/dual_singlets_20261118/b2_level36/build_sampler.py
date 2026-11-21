"""Build an isolated BMN engine from the archived so7_sp3_act2 source."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess


def build(source, output):
    here = Path(__file__).resolve().parent
    output.mkdir(parents=True, exist_ok=True)
    filenames = ("batch.h", "core.h", "d1.h", "q1.h", "trace.h",
                 "linalg.h", "engine.cpp", "levels.h")
    for name in filenames:
        shutil.copyfile(source / name, output / name)
    trace = output / "trace.h"
    levels = output / "levels.h"
    assert hashlib.sha256(trace.read_bytes()).hexdigest() == "57eaf1989bf2dd2fdd7861a896c0095bae5effd17f20aad940922655ba95d900"
    assert hashlib.sha256(levels.read_bytes()).hexdigest() == "ea2d3498ff69cdd3ff8d0ae3abd39db797caab2746555b51b34cb1c9158b4161"
    shutil.copyfile(here / "bmn_sampler.h", output / "bmn_sampler.h")
    with trace.open("a") as stream:
        stream.write('\n#include "bmn_sampler.h"\n')
    old = levels.read_text()
    assert old.count("!random_gen(sec, S, G)") == 1
    levels.write_text(old.replace("!random_gen(sec, S, G)", "!bmn_gen(sec, S, G)"))
    executable = output / "bmn_engine"
    subprocess.run(["c++", "-O3", "-std=c++17", "-pthread", str(output / "engine.cpp"),
                    "-o", str(executable)], check=True)
    return executable


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.source.resolve(), args.output.resolve())
