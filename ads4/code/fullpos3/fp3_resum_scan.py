"""
fullpos3/fp3_resum_scan.py (FULLPOS3.md §2) [DV]: coordinator's resummed-ladder test, integrated.
  R_M(nu) = sum_{k>=0} g_k(nu) c^k exp[-sum_{l<=M} h^l L_l(k)],  L_l(k) = ((c^{-l}-1)/l) sum_{j<k} v_j^l,  v_j = j(j+2b+1),
  h = 1/(T - b(b+1)), c = (C - b(b+1))/(T - b(b+1)); no termination at k = n.   M = 0 is the Legendre kernel at c.
  exact bare kernel K3(nu) = sum_{k<=n} U_k g_k(nu),  g_k = [z^k](1-z)^{-3} P_{-1/2+i nu}((1+z)/(1-z)).
Measured: I_M = 2 int_0^1 |psi2(p)| w^3 |K3 - R_M|(Lam p) dp / N_s  (N_s: T0 hyb_fast, largespin/ls_core.py read-only),
and the signed 2 int psi2 w^3 K3 vs 2 int psi2 w^3 R_M.   Midpoint rule in p (npts), arb with radii discarded (DV).
Usage: fp3_resum_scan.py J,Delta,b [npts] [Mmax]
"""
import sys, os, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import numpy as np
from flint import arb, ctx
import ls_core as LSC

Dg = 1000; lam = (Dg - 2)*(Dg - 1); Lam = math.sqrt(lam)
Fw = LSC.NFw(Dg)

def run(J, Delta, b, npts=600, Mmax=6):
    t0 = time.time()
    n = J - b; C = J*(J + 1); T = (Delta - 1)*(Delta - 2); u0 = b*(b + 1)
    numax = Lam*(1 - 0.5/npts)
    ctx.prec = int((0.75*numax + 80)*3.33) + 64
    X = arb(T - u0); h = 1/X; c = arb(C - u0)/X
    w = 1 - arb(C)/T
    lc = -math.log(float(c.mid()))
    Kmax = max(n + 60, int((math.pi/2*numax + 120)/lc) + 50)
    # ladder weights
    U = []; Uk = arb(1)
    for k in range(n + 1):
        U.append(Uk); u = (b + k)*(b + k + 1); Uk = (Uk*(C - u)/(T - u)).mid()
    WM = [[None]*(Kmax + 1) for _ in range(Mmax + 1)]
    S = [arb(0)]*(Mmax + 1); ck = arb(1)
    coef = [None] + [(c**(-l) - 1)/l*h**l for l in range(1, Mmax + 1)]
    for k in range(Kmax + 1):
        expo = arb(0)
        for M in range(Mmax + 1):
            if M >= 1:
                expo = expo + coef[M]*S[M]
            WM[M][k] = (ck*(-expo).exp()).mid()
        v = arb(k*(k + 2*b + 1)); vl = arb(1)
        for l in range(1, Mmax + 1):
            vl = vl*v; S[l] = S[l] + vl
        ck = (ck*c).mid()
    ps = (np.arange(npts) + 0.5)/npts
    psi = Fw.psi2(ps)
    K3 = np.zeros(npts); RM = np.zeros((Mmax + 1, npts))
    for i, p in enumerate(ps):
        nu = Lam*p
        old = ctx.prec; ctx.prec = int((0.75*nu + 80)*3.33) + 64
        s = arb(nu)**2 + arb(1)/4
        g0, gm = arb(1), arb(0)
        ex = arb(0); acc = [arb(0)]*(Mmax + 1)
        for k in range(Kmax + 1):
            if k <= n:
                ex = ex + U[k]*g0
            for M in range(Mmax + 1):
                acc[M] = acc[M] + WM[M][k]*g0
            g1 = (((2*k*k + 6*k + 3 - s)*g0 - (k + 2)**2*gm)/((k + 1)**2)).mid()
            gm, g0 = g0, g1
        K3[i] = float(ex.mid()); RM[:, i] = [float(a.mid()) for a in acc]
        ctx.prec = old
    m = math.sqrt(T - C); hy = LSC.hyb_fast(Dg, J, m/Lam, '++' if b == 0 else '+-')
    Ns = hy['Ns']; w3 = float(w.mid())**3
    L = -math.log(1 - float(c.mid()))
    I_ex = 2*np.mean(psi*K3)*w3
    out = [f'J={J} Delta={Delta} b={b} a={Delta-J-2} c={float(c.mid()):.5f} L={L:.3f} bhat={(2*J+1)*Lam/m:.1f} h={float(h.mid()):.3e} '
           f'n={n} Kmax={Kmax} npts={npts}: N_s={Ns:.4e}  2int psi2 w3 K3={I_ex:.6e}  [{time.time()-t0:.0f}s]']
    for M in range(Mmax + 1):
        D = np.abs(K3 - RM[M]); I_abs = 2*np.mean(np.abs(psi)*D)*w3; I_s = 2*np.mean(psi*RM[M])*w3
        jmax = int(np.argmax(np.abs(psi)*D))
        out.append(f'   M={M}: int|psi||K3-R_M|/N_s={I_abs/Ns:.3e}  (signed int psi R_M - int psi K3)/N_s={(I_s-I_ex)/Ns:+.3e}  '
                   f'max rel err={np.max(D/np.maximum(np.abs(K3),1e-300)):.2e} worst p={ps[jmax]:.3f}')
    return '\n'.join(out)

if __name__ == '__main__':
    st = sys.argv[1]; npts = int(sys.argv[2]) if len(sys.argv) > 2 else 600; Mm = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    J, D, b = (int(x) for x in st.split(','))
    s = run(J, D, b, npts, Mm)
    print(s, flush=True)
    open(os.path.join(HERE, 'cache', 'fp3_resum_scan.txt'), 'a').write(s + '\n')
