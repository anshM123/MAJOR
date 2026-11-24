"""gaplarge/gl4b_c2conv.py (GAP_LARGESPIN.md §2.2): Mehler-node convergence check of gl4_tol at the J* image (beta = 0.37, 1.203 Lam):
NPM = node multiplier.  Converged at 1000 and 2000; at 5000 the C2 / vector entries move with NPM (likely the finite-difference d/dbeta with
h = 1e-4 beta loses accuracy at nu ~ Lam: error ~ (h nu)^2) -> those entries are not reported at 5000.  -> cache/gl4_c2conv.log"""
import os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import gl4_tol as T
for Dg, J, mu, npms in ((1000, 223, None, (1, 2, 4)), (2000, 447, 1.203, (1, 2, 4)), (5000, 1118, 1.203, (1, 2, 2.5))):
    g = T.Gap(Dg)
    m = mu*g.L if mu else np.sqrt((Dg-2.0)*(Dg - 1.0 + 2*J))
    beta = T.C.beta_of(J, m); F = T.fU(beta)/g.L**2
    for npm in npms:
        T.NPM = npm
        t = g.tol(J, m, F)
        print(f'Dg={Dg} NPM={npm} J={J} beta={beta:.4f}: vec {t["Lam_vec"]:.1f} KC2 {t["KC2"]:.1f} Ksh1 {t["Ksh1"]:.4g} mono {t["Lam_mono"]:.1f}', flush=True)
