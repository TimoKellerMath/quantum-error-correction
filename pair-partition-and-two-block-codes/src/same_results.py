#!/usr/bin/env python3
"""Compare a rerun of exactA or exactB with our log of the same run.

    python3 src/same_results.py RERUN.log OUR.log

The RESULT line without its timing fields (wall, cpu, search, setup), every WORD line and the CONCLUSION line must
agree (both programs are deterministic for a fixed instance, generator file and number of threads). Exit status 1
otherwise.
"""
import re
import sys

TIMING = re.compile(r"^(wall|cpu|search|setup)=")


def key(fn):
    res, words, concl = None, [], None
    with open(fn, encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith("RESULT"):
                res = " ".join(t for t in ln.split() if not TIMING.match(t))
            elif ln.startswith("WORD"):
                words.append(ln)
            elif ln.startswith("CONCLUSION"):
                concl = ln
    return res, words, concl


def main():
    new, old = key(sys.argv[1]), key(sys.argv[2])
    if new[0] is None or new != old:
        print(f"DIFFERENT: {sys.argv[1]} against {sys.argv[2]}\n  rerun: {new}\n  log:   {old}")
        return 1
    print(f"same: {sys.argv[1]} ({new[0]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
