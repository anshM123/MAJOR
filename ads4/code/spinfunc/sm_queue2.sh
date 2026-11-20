#!/bin/bash
# SMOOTH_BOUND.md §2.2: Delta_g = 300 near-gap exact scalar with small-N0 n-sum configurations (h13's default NCONF diverges at
# nu/Lam >~ 0.85 for J >~ 40 there, sm_nsum.py).  Key tag 'n0s'.  Then a p-node convergence check at (1040, 40) (tag 'npf2').
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
SM_NCONF=80:960:1,160:960:1,120:960:1 SM_KEYTAG=n0s nice -n 5 python3 sm_exact.py 300 24,40,64,100
SM_NCONF=80:960:1,160:960:1,120:960:1 SM_KEYTAG=n0s nice -n 5 python3 sm_exact.py 1000 64
SM_NPF=2 SM_KEYTAG=npf2 nice -n 5 python3 sm_exact.py 1000 40
echo "# queue2 done"
