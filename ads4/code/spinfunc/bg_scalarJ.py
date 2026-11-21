"""
bg_scalarJ (BELOW_FOLD2.md §3.2): the scalar analytic-radial expansion RESUMMED in J a = J/Delta (large spin).
P_J(t), t = 1 + delta, delta = a^2 (tau2 + a tau3 + a^2 tau4 + ...):  P_J(1 + delta) = sum_n c_n(J) a^{2n} Dhat^n,  c_n = (J+n)!/((J-n)! 2^n n!^2).
At fixed J a the n-series is not truncated by the a-order: keep n <= N (until converged) and the residual order r <= R in
  I = sum_{n<=N} sum_{r<=R} c_n(J) Delta^{-(2n+r)} sum_j (i nu/2)^j M_{n,r,j}(nu),  M = int lead e^{-S} [Dhat^n x corr_j x B-corr]_{a^r}.
Usage: python3 bg_scalarJ.py N R   (builds cache/bg_scalarJ_N{N}_R{R}.pkl on the validation nu grid)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
import bg_ser as S, bg_geom as GG, bg_scalar as BS

C = os.path.join(HERE, 'cache')


def cJ(J, n):
    if n == 0:
        return 1.0
    v = 1.0
    for i in range(1, n + 1):
        v *= (J*(J + 1) - i*(i - 1))/(2.0*i*i)
    return v


def build(N, R, nus, verbose=True):
    t = time.time()
    Kg = R + 3
    logs, polys, G, _ = GG.scalar_pieces(Kg + 2)
    A = GG.assemble(logs, polys, R)
    t0, dl = GG.spin_t(G, R + 2)
    Dh = [dl.c[k + 2] for k in range(R + 1)]                 # delta = a^2 Dhat
    Bf = GG.bfac(R)
    base = {j: S.lmul(cj, Bf.c, R) for j, cj in A['corr'].items()}
    out = dict(N=N, R=R, nus=np.asarray(nus), E0=A['E0'], E1=A['E1'], M={}, nterms={})
    pw = [S.Poly.const(1.0)] + [S.Poly() for _ in range(R)]
    for n in range(0, N + 1):
        if n:
            pw = S.lmul(pw, Dh, R)
        for j, b in base.items():
            P = S.lmul(pw, b, R)
            for r in range(R + 1):
                out['M'][(n, r, j)] = S.integrate_fast(P[r].E, P[r].C, A['E0'], A['E1'], out['nus']) if len(P[r]) else 0*out['nus']
                out['nterms'][(n, r, j)] = len(P[r])
        if verbose:
            print(f'  n={n}: terms {[len(p) for p in pw]}  {time.time() - t:.0f}s', flush=True)
    return out


def ratio(r, Dl, J, Dref=2e4, Nt=None, Rt=None):
    Nt = r['N'] if Nt is None else Nt; Rt = r['R'] if Rt is None else Rt
    def val(D):
        v = 0
        for (n, rr, j), M in r['M'].items():
            if n <= Nt and rr <= Rt:
                v = v + cJ(J, n)*D**(-(2*n + rr))*(0.5j*r['nus'])**j*M
        return v
    return val(Dl)/val(Dref)


if __name__ == '__main__':
    N, R = int(sys.argv[1]), int(sys.argv[2])
    import bg_validate as BV
    nus = BV.NUSV
    r = build(N, R, nus)
    pickle.dump(r, open(os.path.join(C, f'bg_scalarJ_N{N}_R{R}.pkl'), 'wb'))
