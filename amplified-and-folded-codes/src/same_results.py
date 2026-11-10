#!/usr/bin/env python3
"""Compare the results of a rerun with one of our logs.

    python3 src/same_results.py RERUN LOG

Result lines are the lines starting with "RESULT", "t=", "witness", "logical:", "weight ", "n " or "code ", with
timings ("cpu 1.2 s") and file names removed. The check passes if the result lines of RERUN occur, in order and
contiguously, among the result lines of LOG (a log may hold several runs).
"""
import re
import sys

KEEP = ("RESULT", "t=", "  witness", "witness", "logical:", "weight ", "n ", "code ")
CPU = re.compile(r",? ?cpu [0-9.]+ s|cumulative [0-9.]+ s|\([0-9.]+ ms/trial\)")
PATH = re.compile(r"\S+\.code\b")


def results(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(KEEP):
                out.append(PATH.sub("CODE", CPU.sub("", line.rstrip())))
    return out


def main():
    a, b = results(sys.argv[1]), results(sys.argv[2])
    ok = bool(a) and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
    print(("same results: " if ok else "DIFFERENT results: ") + f"{sys.argv[1]} ({len(a)} lines) vs {sys.argv[2]}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
