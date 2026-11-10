#!/bin/sh
# Validation on the gross code (d = 12) and the two-gross code (d = 18): nothing at W = d - 2, a logical
# operator of weight d at W = d, for both programs. Logs go to rerun/validation.log.
cd "$(dirname "$0")/.." || exit 1
mkdir -p rerun
# job files for vsearch_orb -L as used for C600: Lemma 2.4 for the translations, whose orbits are the two blocks
echo "0 0" > rerun/jobs_V144.txt; echo "72 72 $(seq -s ' ' 0 71)" >> rerun/jobs_V144.txt
echo "0 0" > rerun/jobs_V288.txt; echo "144 144 $(seq -s ' ' 0 143)" >> rerun/jobs_V288.txt
for c in "V144 10 12" "V288 16 18"; do
  set -- $c
  for W in $2 $3; do
    echo "== $1 W=$W"
    bin/vsearch_orb instances/vinst_$1_Z.txt -K num -W $W 2>&1 | grep VRESULT
    bin/vsearch_orb instances/vinst_$1_Z.txt -L rerun/jobs_$1.txt -W $W 2>&1 | grep VRESULT
    bin/cluster2 instances/cinst_$1_Z.txt -W $W -a 0 -m 2 -b 2 | grep RESULT
    bin/cluster2 instances/cinst_$1_Z.txt -W $W -R -m 2 -b 2 | grep RESULT
  done
done | tee rerun/validation.log
[ "$(grep -c ' found=0 ' rerun/validation.log)" -eq 8 ] && [ "$(grep -c ' found=1 ' rerun/validation.log)" -eq 8 ] || { echo 'VALIDATION FAILED' >&2; exit 1; }
# node counts are deterministic: everything but the times must agree with runs/validation.log
sed 's/ time=[0-9.]*//' runs/validation.log > rerun/validation.expected
sed 's/ time=[0-9.]*//' rerun/validation.log > rerun/validation.got
if cmp -s rerun/validation.expected rerun/validation.got; then
  echo "VALIDATION OK"
else
  echo "VALIDATION FAILED: rerun/validation.log differs from runs/validation.log"
  exit 1
fi
