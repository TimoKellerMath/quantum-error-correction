#!/usr/bin/env python3
"""Check that perm.json relabels the board entry codes/480-16-26.json onto the tropidist inputs.

    python3 check_perm.py PATH/codes/480-16-26.json

perm.json holds pi with pi[i] = qubit of the inputs (inputs/480-16-26.hx, .hz) that corresponds to qubit i of the
entry. The check: pi is a permutation of 0..479, and it maps the rows of H_X and of H_Z of the entry onto the rows of
inputs/480-16-26.hx and .hz (as sets of rows). A qubit permutation preserves the minimum distance, so the tropidist
verdicts in records/ on the inputs hold for the entry. Standard library only."""
import json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
def rows(path):
    lines = open(path).read().split('\n'); n, m = map(int, lines[0].split())
    return n, {tuple(sorted(map(int, l.split()))) for l in lines[1:1 + m]}
entry = json.load(open(sys.argv[1]))
pi = json.load(open(os.path.join(here, 'perm.json')))['pi']
n, hx = rows(os.path.join(here, 'inputs', '480-16-26.hx')); _, hz = rows(os.path.join(here, 'inputs', '480-16-26.hz'))
ok = sorted(pi) == list(range(n)) == list(range(entry['n']))
for key, target in (('X', hx), ('Z', hz)):
    mapped = {tuple(sorted(pi[q] for q in r)) for r in entry['checks'][key]}
    ok &= len(mapped) == len(entry['checks'][key]) and mapped == target
print('perm maps the entry onto the inputs:', ok)
sys.exit(0 if ok else 1)
