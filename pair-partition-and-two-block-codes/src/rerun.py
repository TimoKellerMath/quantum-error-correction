#!/usr/bin/env python3
"""Rerun the exact searches recorded in <slug>/records/ and compare them with our logs.

    python3 src/rerun.py [SLUG ...]          (after make bin; all six codes if no SLUG is given)

For each run of <slug>/records/searches.json, the instance and generator files are copied to rerun/<slug>/ and
bin/exactA or bin/exactB is started there with the recorded arguments (one thread). The new log is compared with
ours by src/same_results.py: the RESULT line without its timing fields, every WORD line and the CONCLUSION line
must agree, which they do because both programs are deterministic for a fixed instance, generator file and number
of threads. Exit status 1 if a run differs.
"""
import json
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from same_results import key  # noqa: E402

SLUGS = ["294-42-11", "322-46-12", "434-62-13", "518-74-14", "742-106-15", "672-28-24"]


def main():
    slugs = sys.argv[1:] or SLUGS
    bad = 0
    for slug in slugs:
        rec = os.path.join(ROOT, slug, "records")
        with open(os.path.join(rec, "searches.json"), encoding="utf-8") as f:
            runs = json.load(f)["runs"]
        out = os.path.join(ROOT, "rerun", slug)
        os.makedirs(out, exist_ok=True)
        for r in runs:
            for a in r["args"]:
                if os.path.exists(os.path.join(rec, a)):
                    shutil.copyfile(os.path.join(rec, a), os.path.join(out, a))
            prog = os.path.join(ROOT, "bin", "exactA" if r["impl"] == "A" else "exactB")
            log = os.path.join(out, r["log"])
            with open(log, "w", encoding="utf-8") as fh:
                p = subprocess.run([prog] + r["args"], cwd=out, stdout=fh, stderr=subprocess.STDOUT)
            new, old = key(log), key(os.path.join(rec, r["log"]))
            same = p.returncode == 0 and new[0] is not None and new == old
            print(("same" if same else "DIFFERENT") + f": {slug} exact{r['impl']} side {r['side']}: {new[0]}")
            if not same:
                print(f"  ours: {old}")
                bad += 1
    print("rerun: all runs agree with our logs" if not bad else f"rerun: {bad} runs differ")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
