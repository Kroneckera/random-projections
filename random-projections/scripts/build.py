#!/usr/bin/env python3
"""
Simple one‑shot build helper; identical to what `pip install .` will do,
but faster for day‑to‑day hacking.

Run from the repo root:
    $ python scripts/build.py
"""

from pathlib import Path
import multiprocessing as mp
import subprocess
import sys

ROOT   = Path(__file__).resolve().parent.parent
BUILD  = ROOT / "build"

def run(cmd, **kwargs):
    print(" ·", *cmd)
    subprocess.check_call(cmd, **kwargs)

def main() -> None:
    BUILD.mkdir(exist_ok=True)
    run(["cmake", "-S", str(ROOT), "-B", str(BUILD),
         "-DCMAKE_BUILD_TYPE=Release",
         f"-DPYTHON_EXECUTABLE={sys.executable}"])
    run(["cmake", "--build", str(BUILD), "--config", "Release",
         "--", f"-j{mp.cpu_count()}"])

if __name__ == "__main__":
    main()
