"""
fullpos/fp_contain.py: rigour patch check (FULLPOS_CERT.md §6).  The first moment runs called acb.integral (for K_0 and the two
log-piece integrals) without the holomorphy check on log(1+z^2).  The patched runs (suffix _an) use log(..., analytic=True).
This script verifies that every patched (rigorous) ball is CONTAINED in the corresponding ball used by the certificates, so every
certificate built on the old balls remains valid.
Run: python3 fp_contain.py old.json new.json
"""
import sys, json
sys.dont_write_bytecode = True
from flint import arb, ctx

ctx.prec = 7000
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
n = min(len(a['Z']), len(b['Z']))
bad = 0; worst = 0.0
for key in ('Z', 'Y'):
    for k in range(n):
        old = arb(a[key][k]); new = arb(b[key][k])
        if not old.contains(new):
            bad += 1
            if bad < 5:
                print('NOT contained', key, k)
print(f'{sys.argv[1]} vs {sys.argv[2]}: k <= {n-1}: {2*n - bad}/{2*n} patched balls contained in the old balls')
