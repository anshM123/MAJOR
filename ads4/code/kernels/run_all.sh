#!/bin/sh
# Reproduces every number quoted in ads4/KERNELS_D3.md (about 10 minutes on one core; python3 with sympy, mpmath,
# python-flint, scipy).  Intermediate results are cached in ./cache (override with KERNELS_CACHE).
set -e
cd "$(dirname "$0")"
python3 fock_free.py                 # Fock space: free scalar and free Majorana fermion
python3 free_checks.py 6             # position space -> light transform -> Fourier transform, free fields
python3 ttt.py                       # <TTT>: Ward identity, c, parity-odd structure
for J in 0 1 2 3 4 5 6 7 8 9 10 11 12; do python3 compute_kernels.py $J > /dev/null; done
python3 verify_closed_forms.py       # universal m-dependence (R0), (R4), (Z)
python3 canonical.py                 # canonical structures, closed-form kernels, data/structures_J*.json
python3 free_canonical.py            # closed forms at the unitarity bound against free currents
python3 family4.py                   # exact four-operator solutions of the sum rules
python3 lp_balance.py 10 36          # linear programme (gap, Jmax)
python3 lp_balance.py 15 40
python3 lp_balance.py 20 60
