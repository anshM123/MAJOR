"""deshift/ab_check.py -- certification of an analytic-basis de-shift weight (cache/ab_<piece>.json), DESHIFT_CLOSE.md §2.

(1) residue part (must equal Phi exactly): constant and t^2 coefficients, a2;
(2) eps_W = int_0^oo |E_even(mu)| dmu on an INDEPENDENT fine grid (step 0.01, to mu = 110), with the tail beyond 80 reported separately;
    eps_W bounds |W - Phi| for all t >= 0 and |image - target| / P_{-1/2}(cosh beta) for all beta;
(3) moment bounds  int_0^oo |E_even| mu^{2n} dmu, n = 0..8 (light-side (P) moments of degree 2n are exact up to these);
(4) direct check of the shift formula: W(t) by real-axis quadrature (trapezoid on nu in [-110, 110], step 0.01; no shifting) vs
    Phi(t) + 1/2 int E cos(mu t) dmu, t = 0 ... 10;
(5) comparison with the stored Stage-7 XP weight on its own grid nu_k = (k+1/2)/25: max |w - w_XP| / max |w_XP|, and the change of the
    stored functional's per-state image bounded by  dA = int |w - w_XP| sum_j |R_j| dnu  (C4)  or  int |c - c_XP| sum_j |hhat_j| (C5);
(6) amplification A0 = int |w| sum_j |R_j| (C4) / int |c| sum_j |hhat_j| (C5), and the weight's decay;
(7) the HK cutoff: c(nu) = 1/2[erf((nu+150)/6) - erf((nu-150)/6)] changes w by <= max_{nu} |w (c - 1)| (reported).
Usage: python3 ab_check.py <piece>  -> cache/ab_check_<piece>.txt
"""
import os, sys, time, json
import numpy as np
from flint import acb, arb, ctx
import kern
import ab_solve as AB

ctx.prec = int(os.environ.get('AB_PREC', '384'))
I = acb(0, 1)


def load(piece):
    d = json.load(open(os.path.join(kern.CACHE, f'ab_{piece}.json')))
    P = AB.Piece(piece, d['h'], d['sigma'], d['K'])
    x = [arb(s) for s in d['coeffs']]
    return P, x, d


class Sol:
    def __init__(self, P, x):
        self.P, self.x = P, x
        self.LA = {j: P.laurent(j) for j in (-P.jp, P.jp)}
        self.a = {j: [sum(self.LA[j][i][k]*x[i] for i in range(P.n)) for k in range(P.p)] for j in self.LA}

    def u(self, z):
        return self.P.Q(z)*sum(b*xi for b, xi in zip(self.P.basis(z), self.x))

    def f(self, j, mu):
        z = mu + I*j
        return self.u(z)*self.P.R[j](z)

    def S(self, j, mu):
        a = self.a[j]; e = (-(mu*mu)/2).exp()
        if self.P.p == 1:
            return a[0]*e/mu
        return e*(a[2]/mu**3 + a[1]/mu**2 + (a[0] + a[2]/2)/mu)

    def E(self, mu):
        tot = acb(0)
        for j in self.P.js:
            v = self.f(j, mu)
            if abs(j) == self.P.jp:
                v -= self.S(j, mu)
            tot += v
        return tot

    def Eeven(self, mu):
        return (self.E(acb(mu)) + self.E(acb(-mu)))/2

    def weight(self, nu):
        z = acb(nu); w = self.u(z)
        return w if self.P.sym == 'even' else w*self.P.Nhat(z)

    def W_direct(self, t, L=None, d=0.01):
        """1/2 sum_j int_R u R_j cos((nu - ij) t) dnu, trapezoid (integrand analytic in a strip, decays like a Gaussian)."""
        t = arb(t); tot = acb(0)
        L = L or min(110.0, float(self.P.h.mid())*self.P.K + 40.0)
        n = int(round(L/d))
        for k in range(-n, n + 1):
            nu = arb(k)*arb(d)
            if k == 0:
                continue            # integrand is O(nu) at 0 for both pieces (u = O(nu^2) / O(nu), R_j = O(1/nu))
            z = acb(nu); uz = self.u(z)
            for j in self.P.js:
                tot += uz*self.P.R[j](z)*((z - I*j)*t).cos()
        return (tot*arb(d)/2)


def main(piece):
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    P, x, d = load(piece)
    S = Sol(P, x)
    pr(f"# ab_check {piece}: h = {d['h']}, sigma = {d['sigma']}, K = {d['K']} (unknowns {P.n}), lam = {d['lam']}, prec {ctx.prec}")
    # (1)
    const = sum(I*arb.pi()*(-1 if j < 0 else 1)/2*S.a[j][0] for j in S.a)
    t2 = sum(I*arb.pi()*(-1 if j < 0 else 1)/2*(-S.a[j][2]/2) for j in S.a) if P.p == 3 else acb(0)
    pr(f"(1) residue part: const = {const.real.str(15, radius=False)} + {const.imag.str(3, radius=False)}i (target {P.Phi[0]}), "
       f"t^2-coeff = {t2.real.str(15, radius=False)} + {t2.imag.str(3, radius=False)}i (target {P.Phi[1]})"
       + (f"; a2 = " + ", ".join(S.a[j][1].str(3, radius=False) for j in S.a) if P.p == 3 else ""))
    dconst = abs(float((const - P.Phi[0]).real.mid())) + abs(float(const.imag.mid())) + abs(float((t2 - P.Phi[1]).real.mid())) \
        + abs(float(t2.imag.mid()))
    # (2), (3)
    dm = 0.01
    mus = (np.arange(int(110/dm)) + 0.5)*dm
    Ev = np.empty(len(mus)); Eim = np.empty(len(mus)); Ec = np.empty(len(mus), complex)
    for i, mu in enumerate(mus):
        e = S.Eeven(arb(repr(float(mu))))
        Ec[i] = complex(e.mid()); Ev[i] = abs(Ec[i]); Eim[i] = float(e.imag.mid())
    epsW = float(np.sum(Ev)*dm); tail = float(np.sum(Ev[mus > 80])*dm)
    pr(f"(2) eps_W = int_0^110 |E_even| = {epsW:.3e}  (part beyond mu = 80: {tail:.1e}; |E_even| at 100, 110: {Ev[int(100/dm)]:.1e}, "
       f"{Ev[-1]:.1e});  max|E_even| = {Ev.max():.2e} at mu = {mus[np.argmax(Ev)]:.2f};  max |Im E_even| = {np.max(np.abs(Eim)):.1e}  "
       f"[{time.time() - t0:.0f}s]")
    pr("    |E_even| profile: " + ", ".join(f"mu={m}: {Ev[int(m/dm)]:.1e}" for m in (0.005, 1, 2, 5, 10, 20, 30, 40, 60, 80)))
    pr(f"    => |W(t) - Phi(t)| <= {epsW + dconst:.2e} for ALL t >= 0;  |image - target| <= {epsW + dconst:.2e} P_(-1/2)(cosh beta) for "
       f"ALL beta (value/C5); for the derivative the target is -D(beta) and the bound is the same multiple of P_(-1/2)")
    mom = [float(np.sum(Ev*mus**(2*n))*dm) for n in range(9)]
    pr("(3) moment bounds int|E_even| mu^(2n), n = 0..8: " + ", ".join(f"{m:.1e}" for m in mom))
    # (4)
    pr("(4) shift formula vs direct real-axis W(t):")
    Eeval = Ev*0
    for t in [float(v) for v in os.environ.get('AB_TLIST', '0,0.5,1,2,4,6,8,10').split(',') if v]:
        Wd = S.W_direct(t)
        phi = P.Phi[0] + P.Phi[1]*t*t
        # 1/2 int_R E_even cos = int_0^oo E_even cos  (fine grid, complex E_even)
        corr = float(np.sum(Ec.real*np.cos(mus*t))*dm)
        pr(f"    t = {t:5.1f}: W_direct = {Wd.real.str(20, radius=False)} (Im {Wd.imag.str(2, radius=False)});  W - Phi = "
           f"{float((Wd.real - phi).mid()):+.6e};  shift formula int_0^oo E_even cos(mu t) = {corr:+.6e}   [{time.time() - t0:.0f}s]")
    # (5), (6), (7)
    X = kern.load_xp()
    key = {'c4v': 'C4_value', 'c4d': 'C4_deriv', 'c5': 'C5'}[piece]
    nug, wg = X[key]['nu'], X[key]['wf']
    wn = np.array([complex(S.weight(arb(repr(float(v)))).mid()) for v in nug])
    if P.sym == 'even':
        wn = wn.real
        dens = sum(np.abs(P.R[j].np(nug)) for j in P.js)
    else:
        nh = P.Nhat.np(nug)
        dens = sum(np.abs(P.R[j].np(nug)/nh) for j in P.js)
    dw = np.abs(wn - wg)
    sc = np.max(np.abs(wg))
    dnu = nug[1] - nug[0]
    pr(f"(5) vs stored XP weight ({key}, grid 1/25, nu <= {nug[-1]:.2f}): max|w - w_XP|/max|w_XP| = {dw.max()/sc:.2e} (at nu = "
       f"{nug[np.argmax(dw)]:.2f}); on nu <= 30: {dw[nug <= 30].max()/sc:.2e};  dA = int |w - w_XP| sum|R_j| = {np.sum(dw*dens)*dnu:.2e}")
    A0n = np.sum(np.abs(wn)*dens)*dnu; A0x = np.sum(np.abs(wg)*dens)*dnu
    pr(f"(6) amplification A0: new {A0n:.6f}, XP {A0x:.6f};  |w|/max|w| at nu = 20, 40, 60, 80, 100, 120: " +
       ", ".join(f"{abs(complex(S.weight(arb(v)).mid()))/sc:.1e}" for v in (20, 40, 60, 80, 100, 120)))
    cmax = 0.0
    for v in np.arange(100, 200, 2.0):
        c = ((arb(v) + 150)/6).erf()/2 - ((arb(v) - 150)/6).erf()/2
        cmax = max(cmax, abs(complex(S.weight(arb(v)).mid()))*abs(float((c - 1).mid())))
    pr(f"(7) HK cutoff (nu_c = 150, a = 6): max_nu |w (c - 1)| / max|w| <= {cmax/sc:.1e} (nu in [100, 200]; below 100 |c - 1| < e^-300)")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(kern.CACHE, f'ab_check_{piece}.txt'), 'w').write("\n".join(out) + "\n")


if __name__ == '__main__':
    main(sys.argv[1])
