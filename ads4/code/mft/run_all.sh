#!/bin/sh
# Regenerates every number of ads4/SPINNING_MFT.md.  Exact Fock-space data are cached in cache/ (git-ignored).
# Run time (4 cores, empty cache): mft_exact.py about 2 min; lightcone_check.py and ope_odd_check.py a few minutes each.
# Requires python-flint, sympy, mpmath, and the code (not the caches) of ../e2e.
set -e
cd "$(dirname "$0")"
python3 mft_exact.py > out_mft_exact.txt      # exact MFT OPE coefficients vs closed forms (145 primaries)
python3 odd_checks.py > out_odd_checks.txt    # parity-odd structure: D~ identity, infinity and heavy-frame limits
python3 mft_heavy.py > /dev/null              # heavy limit, MFT normalisation, end-to-end with e2e tables -> out_mft_heavy.txt
python3 goal2.py > out_goal2.txt                # spin dependence of the 1/m^2 coefficient after MFT normalisation
python3 lightcone_check.py > out_lightcone_check.txt    # sect. 5.1: leading-twist block decomposition of the MFT correlator
python3 ope_odd_check.py > out_ope_odd_check.txt      # sect. 5.2: odd normalisation and sign
python3 st_compare.py > out_st_compare.txt            # sect. 3.4: comparison with Sleight-Taronna
echo done
