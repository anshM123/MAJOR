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
python3 g1_crystal4d.py            # M2: D = 4 crystal, helicity content of the CLS-D kernels (seconds)
python3 f7_flat_stripped.py        # M2: flat limits under the stripped dictionary (~1 min)
python3 f6_routeA2.py 100 300 1000 # M2/M3: near-gap model of the route-A' functional psi2 C2 + psi3 C3 (~15 min)
python3 g2_structures.py           # M2b: CLS-D (2,4) structures at the D = 4 crystal: flip-free family, Q^2, O*Q (seconds)
python3 h1_etainv_exact.py 100 300 1000   # M3b: exact eta-inversion of psi_2, checks V1-V3, dk0 (~40 min at 1000)
python3 h2_neargap_exact.py 100 300 1000  # M3b: exact near-gap shift of route A' (needs h1 pickles)
python3 g3_homogeneity.py          # M2c: building-block normalisation check at the crystal (seconds)
python3 g4_images.py               # M2c: independent crystal images of the CLS-D building blocks (seconds)
python3 g5_covariant_Q.py          # M2c: covariant flip-free structure q, dual-coefficient consistency (seconds)
python3 h3_H17_model.py            # M5'b: H17 model estimate for route A' (~1 min)
python3 m1_scalar_scaling.py       # M6: scaling of the exact scalar remainder with m at fixed x (seconds)
python3 m2_tolerance_c2.py 100 300 1000   # M6: H4 tolerance for route A' in c2/m^2 units (~20 min)
python3 l2_scalar_validate.py 100 300 1000   # M6 step 1: Laplace expansion of the CLS-D crystal vs exact scalar skeleton (~10 min)
python3 -c "import l2_scalar_validate as V; V.summary(print)" > cache/l2_richardson.txt   # M6 step 1: Richardson over the gaps
python3 l3_shadow_laplace.py         # M6 step 2 (piece): shadow coefficient from the x5 saddle to O(1/Delta^2), scalar check (~2 min)
