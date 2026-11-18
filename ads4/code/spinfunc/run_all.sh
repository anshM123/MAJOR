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
python3 l4_tto_even.py              # M6 step 2 (piece): conserved parity-even <TTO> structures, d = 3 (~20 min)
python3 l5_spin.py                  # M6 step 2 round 3: contraction / structure machinery self-tests (seconds)
python3 l5_shadow_matrix.py         # item 1: shadow MATRIX beta_cb (parity even), Delta^-3 validation (~10 min)
python3 l7_dc2_scan.py 4 6 8 10 12 0 2   # items 2-3/5 data: graviton C^_2 vs scalar k = 2 crystal actions on (x, Delta) grids (~2 h; resumable)
python3 l7_dc2_analysis.py 0.05     # LO check (2D span, helicity shapes) + per-helicity dc2 fits (UNSTABLE, see doc 21.4)
python3 l7_hel_relative.py 0.05     # helicity-relative correction and mixing (UNSTABLE fits, see doc 21.4)
python3 g6_odd_structure.py > cache/g6_odd_structure.txt   # C^_3 ingredient: covariant odd (2,4) structure, crystal image check (seconds)
