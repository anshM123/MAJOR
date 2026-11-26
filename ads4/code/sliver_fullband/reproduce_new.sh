#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Requires Python 3, GMP C++ headers, MPFR development library or ABI fallback,
# and SymPy for the retained real raw kernel bound.
g++ -O3 -std=c++17 gaussian_seed.cpp -lgmpxx -lgmp -lmpfr -o gaussian_seed
g++ -O3 -std=c++17 gaussian_fixed.cpp -lgmpxx -lgmp -lmpfr -o gaussian_fixed
g++ -O3 -std=c++17 long_transform.cpp -lgmpxx -lgmp -o long_transform
g++ -O3 -std=c++17 direct_sign_triangle.cpp -lgmpxx -lgmp -o direct_sign_triangle
python3 verify_new_seeds.py
python3 spectral_cut.py
python3 full_grid.py
python3 raw_seed_budget_full.py
python3 inventory_cut.py weights.json
# Regenerate full data only if explicitly requested; the published full
# vector contains the original full-scale enclosures.
if [[ "${1:-}" == "--regenerate" ]]; then
 python3 run_fixed_inventory.py 6000 3000 12
 python3 prepare_long.py
 ./long_transform < raw_moments.in > raw_long.txt
 ./long_transform < ideal4.in > ideal_long_corrected.txt
 python3 combine_long.py 17200
 python3 verify_long_transport.py coeff_long_17200.txt
fi
# Reproduce the exact completed scalar checker inputs/outputs, using the
# archived inputs (not a rounded reconstruction of those inputs).
for f in next_129_500_67_250_{0,4}.in; do
 ./direct_sign_triangle < "$f" > "${f%.in}.reproduced.out"
 cmp "${f%.in}.reproduced.out" "${f%.in}.out"
done
python3 verify_new_signs.py
python3 verify_long_transport.py coeff_long_17200.txt
python3 termination_split.py
