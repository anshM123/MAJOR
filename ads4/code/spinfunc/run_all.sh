#!/bin/bash
# spinfunc: route-A milestone M1 (H0_SPINNING_FUNCTIONAL.md). Single core, ~10 min total. Outputs in cache/ (git-ignored).
set -e
cd "$(dirname "$0")"
mkdir -p cache
python3 f1_bulkpoint.py            # symbolic: CLS-D crystal, R_k densities, scalar cross-check   (~1 min)
python3 f2_flat.py                 # flat limit of the D-safe functional with graviton densities   (seconds)
python3 f2b_check.py               # robustness of the flat failure; twist-gap reachability        (~1 min)
python3 f3_ads.py 100 300 1000     # near-gap model, objects A / PsiG / PsiG2 / PsiS; tolerances     (~3 min)
python3 f4_etainv.py 100 300 1000  # eta-inversion cost (large-nu symbol)                           (~3 min)
python3 f5_coherent.py 100 300 1000  # coherent deviation tolerances                               (~3 min)
