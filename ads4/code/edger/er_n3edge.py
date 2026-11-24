"""edger/er_n3edge.py (EDGE_R.md §2.4): Gamma (kappa subtracted, c2) at the n3 pilot nodes of the J = 29 twist-edge row, through the same
hs_core injection as HIGHSPIN (merged n1 + n3 point files, er_row.merge), with the failed-node test restricted to ONE kind at a time
(n3's C3H fails at the state, EDGE_R §2.4, so a node is kept for the kind whose records pass rel. err <= 0.2).  -> stdout"""
import er_row as ER, er_core as R, hs_core as HC, hs_odd as HO, numpy as np
J, D = 29, 1029.0
info = ER.setup(J, D)
_f0 = HC.failed_x
for kind in ('EE', 'C3H'):
    def fx(J_, D_, rtol=0.2, kind=kind):
        m = HC.HN.mphys(D_, J_); bad = set()
        for (DD, nu), r in HC.HN.pts(kind, J_, D_).items():
            if ER.relerr(r) > rtol or not np.all(np.isfinite(r['E'])):
                bad.add(round(nu/m, 6))
        return bad
    HC.failed_x = fx
    rK, rI, isdir = HO.rows(J, D)
    kd = 'C2' if kind == 'EE' else 'C3H'
    print(f"{kind} ({kd}) kept direct nodes x >= 0.9: " + ", ".join(
        f"{x:.6f}: {rI[kd][i]:+.4g}{' (n3)' if info.get(round(x, 6)) else ''}" for i, x in enumerate(rI['xs']) if x >= 0.9 and isdir[i]))
HC.failed_x = _f0
