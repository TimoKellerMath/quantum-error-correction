#!/bin/sh
# Runs every search of src/jobs.txt (name, then the command), P at a time; logs go to rerun/<name>.log.
# Usage: sh src/run_all.sh [P]
cd "$(dirname "$0")/.." || exit 1
P=${1:-4}
mkdir -p rerun
grep -v '^ *$' src/jobs.txt | while read -r name cmd; do
  printf '%s\t%s\n' "$name" "$cmd"
done | xargs -P "$P" -L 1 sh -c 'name=$0; "$@" > rerun/$name.log 2>&1; echo "$name exit=$?"'
# node counts are deterministic: apart from the times and the instance directory (codes/ in our runs),
# every log must agree with the one in runs/
fail=0
for f in runs/A_*.log runs/B_*.log runs/W_*.log; do
  norm='s/ time=[0-9.]*//; s/ [0-9.]* s$//; s#codes/#instances/#'
  sed "$norm" "$f" > rerun/expected.tmp
  sed "$norm" "rerun/${f#runs/}" > rerun/got.tmp 2>/dev/null
  cmp -s rerun/expected.tmp rerun/got.tmp || { echo "DIFFERS: rerun/${f#runs/}"; fail=1; }
done
rm -f rerun/expected.tmp rerun/got.tmp
if [ "$fail" = 0 ]; then echo "SEARCH OK: every log agrees with runs/"; else exit 1; fi
