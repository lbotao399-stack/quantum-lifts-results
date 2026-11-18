"""Build the archived orthogonal trace, extremal-projector and Q1 kernels."""
import argparse
import fcntl
import os
import subprocess
import tempfile
from pathlib import Path


def build(folder):
    builds = [("full_projection_native.cpp", "full_projection_native.dylib", []),
              ("so7_q1_native_v1.cpp", "so7_q1_native_v1.dylib", [])]
    for prime in (1000003, 1000033):
        options = ["-DTRACE_FRAME_PRIME=%d" % prime, "-DTRACE_FRAME_MATRIX_SIZE=7"]
        builds.append(("trace_frame_terminal_contract_v1.cpp",
                       "trace_frame_terminal_contract_v1_d7_p%d.dylib" % prime, options))
    with (folder / ".compile.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for source, binary, options in builds:
            target = folder / binary
            dependencies = [folder / source]
            if source == "full_projection_native.cpp":
                dependencies.append(folder / "q1_native.cpp")
            if target.exists() and target.stat().st_mtime_ns >= max(p.stat().st_mtime_ns for p in dependencies):
                continue
            with tempfile.TemporaryDirectory(dir=folder) as local:
                temporary = Path(local) / binary
                subprocess.run([os.environ.get("CXX", "c++"), "-std=c++17", "-O3", "-shared", "-fPIC",
                                *options, str(folder / source), "-o", str(temporary)], check=True)
                temporary.replace(target)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--native", type=Path, default=Path(__file__).with_name("native"))
    build(p.parse_args().native)
