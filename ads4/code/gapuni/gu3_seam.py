"""
gu3_seam (GAP_UNIFORM.md §3.3): tolerance of the twist-edge ('seam') states that have no Delta_g = 1000 counterpart.
Gamma(x; J, m) is a property of the heavy block alone (it does not know Lam), so a seam state (J, mu_s, Lam') with
mu_s in [mu_edge(Dg'), mu_edge(1000)) needs Gamma on x <= 1/mu_s at m = mu_s Lam'; its TOLERANCE is computed here directly at its own
(mu_s, Lam') with the same chain definitions as gu1_tol (read-only large_m_map.state):
   G_cert = rig/W,  G_hyb = lam/W  (c2),  and their ratio to the Delta_g = 1000 twist-edge tolerance x (Lam'/Lam0)^2.
Also the weight of the extra edge strip x in (1/mu_edge(1000), 1/mu_s] (cost of |Gamma| = 1 c2 there, N_s units, x Lam'^2).
Usage: python3 gu3_seam.py -> cache/gu3_seam.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gu_guard  # noqa: F401
sys.path.insert(0, os.path.join(HERE, '..', 'poserr'))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import large_m_map as LM
E = LM.E
Q = E.Q
C = os.path.join(HERE, 'cache')
OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def strip_weight(Dg, J, mu, x_lo, h):
    """cost of |Gamma| = 1 c2 on x > x_lo (N_s units), same node weights as large_m_map (W)."""
    Lm = E.lam_of(Dg)
    st = E.hyb(Dg, J, mu, h, LM.npn_of(J, mu, Lm))
    x, m, u = st['x'], st['m'], st['u']
    w = (np.abs(st['p2'])*Q.A.env(h, x, 2) + u*np.abs(st['p3'])*Q.A.env(h, x, 3))/m**2/st['Ns']
    return float(w[x > x_lo].sum()), float(w.sum())


def main():
    L0 = E.lam_of(1000.0)
    pr("# seam states: mu_s = twist edge at Dg'; tolerances vs the 1000 twist-edge tolerance scaled by (Lam'/Lam0)^2")
    pr(f"{'J':>3} {'h':>3} {'Dg_':>7} {'mu_s':>9} {'mu_1000':>8} | {'Gcert_s':>9} {'ratio':>7} | {'Ghyb_s':>9} {'ratio':>7} | "
       f"{'W_strip*Lam^2':>13} {'W_tot*Lam^2':>11}")
    for J in [0, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16]:
        mu0 = E.mgap(1000.0, J)*1.000001
        hs = ['++'] if J < 4 else (['+-'] if J % 2 else ['++', '+-'])
        s0 = LM.state(1000.0, J, mu0)
        for Dg in (2000.0, 5000.0, 1.0e5):
            mus = E.mgap(Dg, J)*1.000001
            Lm = E.lam_of(Dg)
            s = LM.state(Dg, J, mus)
            for h in hs:
                W = s[h]['W_lo'] + s[h]['W_hi']; W0 = s0[h]['W_lo'] + s0[h]['W_hi']
                Gc, Gh = s[h]['rig']/W, s[h]['lam']/W
                Gc0, Gh0 = s0[h]['rig']/W0, s0[h]['lam']/W0
                r = (Lm/L0)**2
                ws, wt = strip_weight(Dg, J, mus, 1/mu0, h)
                pr(f"{J:>3} {h:>3} {Dg:7.0f} {mus:9.6f} {mu0:8.5f} | {Gc:9.3e} {Gc/(Gc0*r) if Gc0 > 0 else float('nan'):7.3f} | "
                   f"{Gh:9.3e} {Gh/(Gh0*r):7.3f} | {ws*Lm**2:13.4e} {wt*Lm**2:11.4e}")
    open(os.path.join(C, 'gu3_seam.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
