#!/bin/sh
# Regenerates every number of ads4/SPINNING_MFT.md.  Exact Fock-space data are cached in cache/ (git-ignored).
# Run time (4 cores, empty cache): about 15 min (mft_exact.py, l = 2 dominates); the rest about 2 min.
# Requires python-flint, sympy, mpmath, and the code (not the caches) of ../e2e.
set -e
cd "$(dirname "$0")"
python3 mft_exact.py > out_mft_exact.txt      # exact MFT OPE coefficients vs closed forms (145 primaries)
python3 odd_checks.py > out_odd_checks.txt    # parity-odd structure: D~ identity, infinity and heavy-frame limits
python3 mft_heavy.py > /dev/null              # heavy limit, MFT normalisation, end-to-end with e2e tables -> out_mft_heavy.txt
python3 goal2.py > out_goal2.txt                # spin dependence of the 1/m^2 coefficient after MFT normalisation
echo done
