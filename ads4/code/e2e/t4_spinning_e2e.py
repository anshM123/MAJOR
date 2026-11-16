"""
End-to-end test (tasks 2, 3): tensor Regge moments of the exact spinning partial waves versus the derived leading
heavy-block formula (predict.py), for all Fourier modes (q1,q3) of the stripped frame tensor.
   exact:      (-1)^J Pihat^{(q1,q3)}_{k,nu}[ s13^{l+1} s24^{Dphi} W_diff ]     (moments.TensorMoments; no large-Delta approximation)
   predicted:  predict.prediction   (derived; no free constant)
The Fourier modes of the exact stripped tensor on the integration grid are stored in data_modes/ (one file per block), so that
moments for any (k, nu) can be recomputed without re-evaluating the blocks (t7_corrections.py, t8_uniformity.py).
Usage: python3 t4_spinning_e2e.py l [quick | small]      ('small': the small-impact-parameter cases J = 2, 3, 5 used in t8)
"""
import sys, os, time, numpy as np, mpmath as mp
from moments import TensorMoments
from predict import prediction
from cases import cases_fixed_beta, case, info
mp.mp.dps = 25
KS = (2, 4, 6)
NUS = (0.0, 1.0, 3.0)
DPHI = 1.75
HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data_modes')

def default_N(l, Delta, c=400):
    """radial truncation.  The region of small u = sqrt(rhobar/rho) that can be covered is u >~ 34 alpha/(N t) (moments.series_converged);
    the cost grows like N^3 (more points become usable, each costing N^2), so N is capped.  The neglected region contributes
    ~ u_min^{k+1/2}; it is estimated for every moment by tightening the mask (TensorMoments.tail_estimate) and by the
    N-dependence study t4b_convergence.py."""
    return int(min(10000 if l == 1 else 8000, max(4000, c*Delta)))

def fname(l, Delta, J):
    return os.path.join(DATADIR, f"modes_l{l}_D{Delta:.6f}_J{J}.npz")

def get_moments(l, Delta, J, Dphi=DPHI, N=None, cache=True, **gkw):
    fn = fname(l, Delta, J)
    if cache and os.path.exists(fn):
        return TensorMoments.load(fn)
    N = default_N(l, Delta) if N is None else N
    mom = TensorMoments(l, Delta, J, Dphi, N, **gkw)
    if cache:
        os.makedirs(DATADIR, exist_ok=True)
        mom.save(fn)
    return mom

def moments_of(mom, ks, nus):
    l, J = mom.l, mom.J
    res = {}
    for k in ks:
        for nu in nus:
            for q1 in range(-l, l+1):
                for q3 in range(-l, l+1):
                    val, sh = mom.moment(k, nu, q1, q3)
                    res[(k, nu, q1, q3)] = (-1)**J*mp.mpc(complex(val))*mp.exp(mp.mpf(sh))
    return res

def run_case(l, Delta, J, Dphi, ks, nus, N=None, **gkw):
    """uncached evaluation (used by the convergence study)"""
    N = default_N(l, Delta) if N is None else N
    return moments_of(TensorMoments(l, Delta, J, Dphi, N, **gkw), ks, nus)

def main_cases(quick=False):
    out = []
    for beta in ((1.0,) if quick else (0.4, 1.0, 2.0)):
        mlist = [10, 20, 40] if quick else ([10, 20, 40, 70, 100] if beta < 1.5 else [10, 20, 40, 70])
        out += cases_fixed_beta(beta, mlist)
    return out

def small_cases(l):
    """small impact parameter: J = 2, 3, 5 (also J = 1 for currents) at m ~ 20, 40, 80;  Delta = m + J... chosen as m_target + 3/2"""
    out = []
    for mt in (20.0, 40.0, 80.0):
        for J in ((1, 2, 4) if l == 1 else (2, 3, 5)):
            # alpha alphabar = m^2 with alphabar - alpha = 2J+1
            al = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*mt*mt))/2
            out.append((float(al + J + 2), J))
    return out

if __name__ == "__main__":
    l = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else 'full'
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    modes = [(1, -1), (1, 1), (0, 0), (1, 0)] if l == 1 else [(2, -2), (2, 2), (0, 0), (1, -1), (1, 1), (2, 0), (2, -1), (1, 0)]
    pr(f"# <O phi O phi>, O = conserved spin l = {l} (Delta_O = {l+1}), Dphi = {DPHI}.  Entries: exact/predicted for the mode (q1,q3)")
    pr("# (complex ratio; real part shown, imaginary part printed when > 1e-6).  The exact moments are symmetric under (q1,q3) -> (q3,q1) and (-q1,-q3).")
    caselist = small_cases(l) if mode == 'small' else main_cases(mode == 'quick')
    for (Delta, J) in caselist:
        t0 = time.time()
        m, b_, eta = info(Delta, J)
        mom = get_moments(l, Delta, J, cache=(mode != 'quick'))
        res = moments_of(mom, KS, NUS)
        pr(f"# Delta = {Delta:.5f}, J = {J}, m = {float(m):.2f}, beta = {float(b_):.3f}   ({time.time()-t0:.0f} s)")
        for k in KS:
            for nu in NUS:
                row = []
                for (q1, q3) in modes:
                    ex = res[(k, nu, q1, q3)]
                    pd = prediction(l, Delta, J, DPHI, k, nu, q1, q3)
                    rat = ex/pd if abs(pd) > 0 else mp.nan
                    row.append(f"({q1:+d},{q3:+d}) {float(mp.re(rat)):9.6f}" + ("" if abs(mp.im(rat)) < 1e-6 else f"{float(mp.im(rat)):+.1e}i"))
                pr(f"   k={k} nu={nu:3.1f}: " + "  ".join(row))
    if mode != 'quick':
        open(os.path.join(HERE, f'out_t4_spinning_e2e_l{l}' + ('_small' if mode == 'small' else '') + '.txt'), 'w').write("\n".join(out) + "\n")
