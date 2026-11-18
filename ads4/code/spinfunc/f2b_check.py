"""
f2b: robustness of the f2 finding (flat limit of the current Omega^{NF,SC} with graviton densities negative at m = M, +-, J ~ 20):
quadrature 600 / 2000 / 6000 nodes, J <= 120, fine m/M grid in [1, 1.003]; per-k decomposition at the worst state.
Usage: python3 f2b_check.py -> cache/f2b_check.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from f2_flat import action
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
pr("# f2b: +- flat action of the current functional (PsiG) near m = M; min F/N_s over J <= 120 (J >= 4)")
for npq in (600, 2000, 6000):
    for mfl in (1.0, 1.0001, 1.0003, 1.0005, 1.001, 1.002, 1.003):
        u = 1/mfl**2
        cells = []
        for obj in ('A', 'PsiG'):
            Js, parts, Ns, _ = action('+-', u, 120, obj, npq=npq)
            rat = parts.sum(1)/Ns; ok = Js >= 4
            i = int(np.argmin(rat[ok])); cells.append(f"{obj}: {rat[ok][i]:+.5f} (J={int(Js[ok][i])})")
            if obj == 'PsiG':
                neg = Js[ok][rat[ok] < 0]
                cells.append(f"negative J: {list(map(int, neg))[:12]}{'...' if len(neg) > 12 else ''}")
        pr(f"  npq={npq} m/M={mfl}: " + " | ".join(cells))
pr("# per-k parts / N_s at m = M, +-, J = 16..24 (npq = 6000): A | PsiG")
Js, pA, Ns, _ = action('+-', 1.0, 40, 'A', npq=6000)
Js, pG, Ns, _ = action('+-', 1.0, 40, 'PsiG', npq=6000)
for J in range(14, 27, 2):
    i = int(np.searchsorted(Js, J))
    pr(f"   J={J}: A " + " ".join(f"{v:+.4f}" for v in pA[i]/Ns[i]) + f" = {pA[i].sum()/Ns[i]:+.4f} | PsiG " +
       " ".join(f"{v:+.4f}" for v in pG[i]/Ns[i]) + f" = {pG[i].sum()/Ns[i]:+.4f}")
pr("# same at J = 15..25 odd")
for J in range(15, 26, 2):
    i = int(np.searchsorted(Js, J))
    pr(f"   J={J}: A {pA[i].sum()/Ns[i]:+.4f} | PsiG {pG[i].sum()/Ns[i]:+.4f}  (k=5 part A {pA[i][3]/Ns[i]:+.4f}, PsiG {pG[i][3]/Ns[i]:+.4f})")
open(os.path.join(HERE, 'cache', 'f2b_check.txt'), 'w').write("\n".join(out) + "\n")

# --- where is the flat failure kinematically reachable in AdS with a twist gap?  m_min(J)^2 = (Dg-2)(Dg+2J-1), Lam^2 = (Dg-2)(Dg-1)
pr("# failure window per J: largest m/M with F_PsiG < 0 (bisection), and the Delta_g above which AdS states with that J reach it")
def Fpsig(mfl, J):
    Js, parts, Ns, _ = action('+-', 1/mfl**2, J, 'PsiG', npq=2000)
    return parts[-1].sum()/Ns[-1]
for J in range(10, 34, 2):
    if Fpsig(1.0, J) >= 0:
        pr(f"   J={J}: F(m=M) = {Fpsig(1.0, J):+.5f} >= 0, no window"); continue
    lo, hi = 1.0, 1.01
    for _ in range(40):
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if Fpsig(mid, J) < 0 else (lo, mid)
    mstar = lo
    Dgc = 1 + 2*J/(mstar**2 - 1)
    pr(f"   J={J}: F(m=M) = {Fpsig(1.0, J):+.5f}; window m/M < {mstar:.6f}; reachable only for Delta_g > {Dgc:.3g}")
open(os.path.join(HERE, 'cache', 'f2b_check.txt'), 'w').write("\n".join(out) + "\n")
