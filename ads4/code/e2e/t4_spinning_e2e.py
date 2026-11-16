"""
End-to-end test (tasks 2, 3): tensor Regge moments of the exact spinning partial waves versus the derived leading
heavy-block formula (predict.py), for all Fourier modes (q1,q3) of the stripped frame tensor.
   exact:      (-1)^J Pihat^{(q1,q3)}_{k,nu}[ s13^{l+1} s24^{Dphi} W_diff ]     (moments.TensorMoments; no large-Delta approximation)
   predicted:  predict.prediction   (derived; no free constant)
Usage: python3 t4_spinning_e2e.py l [quick]         (writes out_t4_spinning_e2e_l{l}.txt and data_t4_l{l}.pkl)
"""
import sys, time, pickle, numpy as np, mpmath as mp
from moments import TensorMoments
from predict import prediction
from cases import cases_fixed_beta, info
mp.mp.dps = 25
KS = (2, 4, 6)
NUS = (0.0, 1.0, 3.0)

def default_N(l, Delta, c=160):
    cap = 20000 if l == 1 else 16000
    return int(min(cap, max(2000, c*Delta)))

def run_case(l, Delta, J, Dphi, ks, nus, N=None, **gkw):
    N = default_N(l, Delta) if N is None else N
    mom = TensorMoments(l, Delta, J, Dphi, N, **gkw)
    res = {}
    for k in ks:
        for nu in nus:
            for q1 in range(-l, l+1):
                for q3 in range(-l, l+1):
                    val, sh = mom.moment(k, nu, q1, q3)
                    res[(k, nu, q1, q3)] = (-1)**J*mp.mpc(complex(val))*mp.exp(mp.mpf(sh))
    res['N'] = N
    return res

if __name__ == "__main__":
    l = int(sys.argv[1])
    quick = len(sys.argv) > 2 and sys.argv[2] == 'quick'
    Dphi = 1.75
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    betas = (1.0,) if quick else (0.4, 1.0, 2.0)
    store = {}
    modes = [(1, -1), (1, 1), (0, 0), (1, 0)] if l == 1 else [(2, -2), (2, 2), (0, 0), (1, -1), (1, 1), (2, 0), (2, -1), (1, 0)]
    pr(f"# <O phi O phi>, O = conserved spin l = {l} (Delta_O = {l+1}), Dphi = {Dphi}.  Entries: exact/predicted for the mode (q1,q3)")
    pr("# (complex ratio; real part shown, imaginary part printed when > 1e-6).  The exact moments are symmetric under (q1,q3) -> (q3,q1) and (-q1,-q3).")
    for beta in betas:
        mlist = [10, 20, 40] if quick else ([10, 20, 40, 80, 120] if beta < 1.5 else [10, 20, 40, 80])
        for (Delta, J) in cases_fixed_beta(beta, mlist):
            t0 = time.time()
            m, b_, eta = info(Delta, J)
            res = run_case(l, Delta, J, Dphi, KS, NUS)
            store[(Delta, J)] = res
            pr(f"# Delta = {Delta:.5f}, J = {J}, m = {float(m):.2f}, beta = {float(b_):.3f}   (N = {res['N']}, {time.time()-t0:.0f} s)")
            for k in KS:
                for nu in NUS:
                    row = []
                    for (q1, q3) in modes:
                        ex = res[(k, nu, q1, q3)]
                        pd = prediction(l, Delta, J, Dphi, k, nu, q1, q3)
                        rat = ex/pd if abs(pd) > 0 else mp.nan
                        row.append(f"({q1:+d},{q3:+d}) {float(mp.re(rat)):9.6f}" + ("" if abs(mp.im(rat)) < 1e-6 else f"{float(mp.im(rat)):+.1e}i"))
                    pr(f"   k={k} nu={nu:3.1f}: " + "  ".join(row))
            if not quick:
                pickle.dump(store, open(f'data_t4_l{l}.pkl', 'wb'))
    if not quick:
        open(f'out_t4_spinning_e2e_l{l}.txt', 'w').write("\n".join(out) + "\n")
