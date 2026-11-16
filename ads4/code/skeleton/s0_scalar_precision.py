"""
s0: how far in nu (at fixed m) are double-precision exact Regge moments reliable?  Identical external scalars (Df), exact blocks
(e2e.moments.ScalarMoments), on the standard (u, t) grid and on refined grids; compared with the CHMRS order-8 series.
Output R = exact / leading Regge formula (SPINNING_E2E.md 4.1 with D1 = D2 = Df) and the series ratio truncated at 1/m^8.
Usage: python3 s0_scalar_precision.py Delta J Df [N]
"""
import sys, os, time
from common import CACHE, grid_kw, info
import numpy as np, mpmath as mp
from moments import ScalarMoments
from formulas import gamma_a, log_b, Pnu
import skel
mp.mp.dps = 20

def lead(Delta, J, Df, k, nu):
    m, beta, eta = info(Delta, J)
    a = 2*Df + k - 1
    return mp.exp(mp.log(mp.pi/4) + 2*Df*mp.log(4) - log_b(mp.mpf(Delta), J, 0) - a*mp.log(m*m))*gamma_a(a, nu)**2*Pnu(nu, eta)

def series_ratio(Delta, J, Df, k, nu, jmax=4):
    """CHMRS order-2*jmax series of the moment divided by the leading formula (seriesPi normalisation of uniform/skel)."""
    m, beta, eta = info(Delta, J); m = float(m); eta = float(eta)
    lq, la = skel.log_q(Delta, J, Df), mp.log(skel.a_Df(Df))
    g2 = skel.gamma_a(2*Df + k - 1, nu)**2
    pref = mp.exp(lq - la - k*mp.log(m*m))*g2
    Pd = skel.Pders(nu, eta)
    t = skel.series_terms('seriesPi', k, Df, nu, eta, Pd=Pd)
    L = lead(Delta, J, Df, k, nu)
    return [float(pref*sum(t[i]/(m*m)**i for i in range(j + 1))/mp.re(L)) for j in range(jmax + 1)]

if __name__ == "__main__":
    Delta, J, Df = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    N = int(sys.argv[4]) if len(sys.argv) > 4 else int(min(12000, max(600, 45*Delta)))
    m, beta, eta = info(Delta, J)
    out = [f"# s0: Delta={Delta} J={J} Df={Df} m={float(m):.3f} beta={float(beta):.3f} N={N}"]
    grids = {'std': grid_kw(), 'u2': grid_kw(fu=2), 'u3s2': grid_kw(fu=3, usplit=2), 'u2t2': grid_kw(fu=2, ft=2)}
    moms = {}
    for name, kw in grids.items():
        t0 = time.time(); moms[name] = ScalarMoments(Delta, J, Df, Df, N, **kw); out.append(f"#  grid {name}: {time.time()-t0:.1f} s")
    print("\n".join(out), flush=True)
    for k in (4, 6):
        hdr = f"#  k={k}  nu   x=nu/m   " + "  ".join(f"R[{g}]" .ljust(13) for g in grids) + "   series8      series6"
        print(hdr); out.append(hdr)
        for nu in np.arange(0, 10.01, 1.0):
            L = lead(Delta, J, Df, k, nu)
            Rs = []
            for g in grids:
                s, lg = moms[g].moment(k, nu)
                Rs.append(float(s*mp.exp(lg)/mp.re(L)))
            ser = series_ratio(Delta, J, Df, k, nu)
            line = f"   {nu:5.1f}  {nu/float(m):.3f}  " + "  ".join(f"{r:13.8f}" for r in Rs) + f"  {ser[4]:11.6f}  {ser[3]:11.6f}"
            print(line, flush=True); out.append(line)
    open(os.path.join(CACHE, f"out_s0_D{Delta:.3f}_J{J}_Df{Df}.txt"), 'w').write("\n".join(out) + "\n")
