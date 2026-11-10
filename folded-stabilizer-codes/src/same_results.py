#!/usr/bin/env python3
"""Compare the results of a rerun with one of our logs.

    python3 src/same_results.py RERUN LOG

Result lines: for src/mitm.c the statistics line (without "cpu=...") and the verdict line ("NONE of weight <= w"
or "FOUND weight w: ..."); for Magma the lines "code: ...", "Final Results: ..." and "MinimumDistance: d" (without
the CPU time). The check passes if RERUN and LOG have the same result lines, in the same order.
"""
import re
import sys

CPU = re.compile(r" cpu=[0-9.]+| cputime: [0-9.]+")
KEEP = ("n=", "NONE of weight", "FOUND weight", "code:", "Final Results:", "MinimumDistance:")


def results(path):
    with open(path, encoding="utf-8") as f:
        return [CPU.sub("", line.rstrip()) for line in f if line.startswith(KEEP)]


def main():
    a, b = results(sys.argv[1]), results(sys.argv[2])
    ok = bool(a) and a == b
    print(("same results: " if ok else "DIFFERENT results: ") + f"{sys.argv[1]} ({len(a)} lines) vs {sys.argv[2]}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
