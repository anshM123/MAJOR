"""
cert_quad: certified error bound for the 64-node Gauss-Legendre p-rule in the D-safe action, so that the certified statement of
cert_action.py (F^{(64)} > 0) becomes a statement about the exact p-integral  F = sum_k u^{k-2} int_0^1 phi_k(p) Mod_k(p) dp.

Bound [Trefethen, ATAP Thm 19.3, mapped to [0,1]]: if f is analytic in the Bernstein ellipse E_rho of [0,1] and |f| <= M there, then
   |int_0^1 f - sum_i w_i f(p_i)| <= (1/2)(64/15) M rho^{-2n} / (rho^2 - 1),   n = 64.
Analyticity: f(p) = sum_k u^{k-2} phi_k(p) Mod_k(p) is analytic except at y = u p^2 = 1 (p = +-m/Lam) and the branch cut y in [1, oo) of
(1-y)^{-(k+1)/2} (odd k), i.e. real |p| >= m/Lam; E_rho is chosen with right vertex < m/Lam.  The moments ex_k(nu) = sum_j sum_m h^(j)_m(nu)
S_{j,m} are entire in nu (g_m are polynomials in nu^2, absolutely convergent series).  M is bounded on the boundary of E_rho (maximum
modulus principle) by covering it with K complex balls and using:
  |ex_k| <= |pre| [ sum_{j, m<M} hhat^(j)_m (|S_{j,m}| + U_{j,m}) + sum_j K_j sum_{M<=m<M2} hhat^(j)_m bx_m + far ],  hhat = majorant of h
     built from ghat_m = g_m with nu^2 -> |nu|^2 (|(A/2+i)^2 + nu^2/4| <= (A/2+i)^2 + |nu|^2/4), M2 re-certified for Lam_eff = max|nu|;
  odd k: |ex_k| <= |E_k| sum_l |L_l(k)| |ex_l|/|E_l|;
  |P_{-1/2+i nu}(cosh beta)| <= cosh(|Im nu| beta) P_{-1/2}(cosh beta)  (Mehler-Dirichlet);  |J_0(z)| <= cosh(|Im z|);
  E_k, (1-y)^{-(k+1)/2}, phi_k, P_J, d^J_44 evaluated in acb on the balls.
Usage: python3 cert_quad.py Dg Df [mfl list] [J list]  -> cache/certq_<Dg>_Df<Df>.pkl, .log   (reads cache/cert_<Dg>_Df<Df>.pkl)
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
Dg = int(sys.argv[1]); DF = float(sys.argv[2])
os.environ['OMEGAU_DF'] = '%g' % DF
import numpy as np
from fractions import Fraction as Fr
from flint import arb, acb, ctx
from nf_common import NF, KS, SIGMA
from common import state
from cert_moments import CertMoments, af, poly_nonneg_shift, pmul, ppow, psub

KE = (2, 4, 6, 8, 10)
NQ = 64
F = NF()


def cpow(z, s):
    """z^s (principal branch) as exp(s log z): much tighter balls than acb.__pow__ for non-integer s."""
    if float(s) == int(float(s)):
        n = int(float(s)); r = acb(1); b = z if n >= 0 else 1/z
        for _ in range(abs(n)):
            r = r*b
        return r
    return (af(s)*z.log()).exp()


def jacobi_acb(nmax, a, b, x):
    P = [acb(1)]
    if nmax >= 1:
        P.append((a + 1) + (a + b + 2)*(x - 1)/2)
    for n in range(2, nmax + 1):
        c = 2*n + a + b
        P.append((((c - 1)*(a*a - b*b) + (c - 2)*(c - 1)*c*x)*P[n - 1] - 2*(n + a - 1)*(n + b - 1)*c*P[n - 2])/(2*n*(n + a + b)*(c - 2)))
    return P


_PHC = {}
def phi_acb(k, p):
    """Horner evaluation of the exact phi_k (less ball overestimation than monomial powers)."""
    if k not in _PHC:
        _PHC[k] = [af(c) for c in F.PHx[k]]
    c = _PHC[k]
    r = acb(0)
    for a in reversed(c):
        r = r*p + a
    return r


def ex_abs_bound(cm, k, numax):
    """upper bound (arb) on |m^{2k} C_{k,nu}| for all complex nu with |nu| <= numax (cm: CertMoments with S, U, K for k)."""
    J = cm.J
    A = 2*cm.Dfa + k - arb(1.5)
    gA = A.gamma()
    M = cm.M; x0 = cm.x0; kap = k//2; gk = arb(kap).gamma()
    S, U, K = cm.S[k], cm.U[k], cm.K[k]
    Le = arb(numax)
    # M2: g_m(Le) decreasing for m >= M2 and the certified ratio check
    M2 = int(float((A*A/4 + Le*Le/4 - A).mid())) + 10
    M2 = max(M2, M + 10)
    x0q = (cm.Dq - J)/2
    while True:
        lhs = pmul(pmul(pmul([Fr(0), Fr(0), Fr(1)], ppow([Fr(J), Fr(1)], J)), [x0q + kap, Fr(1)]), [2*x0q + kap - 1, Fr(1)])
        rhs = pmul(pmul(pmul(ppow([Fr(1), Fr(1)], 2), ppow([Fr(1 + J), Fr(1)], J)), [Fr(kap), Fr(1)]), [2*x0q + kap, Fr(1)])
        if poly_nonneg_shift(psub(lhs, rhs), Fr(M2)):
            break
        M2 = int(M2*1.5) + 10
    L = M2 + J + 2
    g = []; gm = 1/gA
    for m in range(L):
        g.append(gm); gm = gm*((A/2 + m)**2 + Le*Le/4)/((A + m)*(m + 1))
    bx = []; b = (x0.lgamma() + arb(kap).lgamma() - (x0 + kap).lgamma()).exp()
    for m in range(M2 + 1):
        bx.append(b*(1 + (x0 + m + kap)/(x0 - 1))/gk); b = b*(m + kap)/(x0 + m + kap)
    tot = arb(0)
    h = g
    for j in range(J + 1):
        if j > 0:
            h = [2*(m + 1)*h[m + 1] + 2*m*h[m] for m in range(len(h) - 1)]      # majorant of |D h|
        for m in range(M):
            tot += h[m]*(abs(S[j][m]) + U[j][m])
        for m in range(M, M2):
            tot += h[m]*K[j]*bx[m]
    far = sum((K[j] for j in range(J + 1)), arb(0))*(4*(arb(M2) + J))**J*g[M2]*bx[M2]*(1 + M2)
    logpre = (arb(2).log()*(1 - 2*k) + arb.pi().log()*arb(1.5) + cm.aDf.log() + k*cm.m2.log() - cm.logq)
    return logpre.exp()*(tot + far)


def run():
    ctx.prec = 192
    Lamf = float(np.sqrt((Dg - 2.0)*(Dg - 1.0)))
    Lam = arb((Dg - 2)*(Dg - 1)).sqrt()
    src = pickle.load(open(os.path.join(HERE, 'cache', 'cert_%d_Df%g.pkl' % (Dg, DF)), 'rb'))
    pk = os.path.join(HERE, 'cache', 'certq_%d_Df%g.pkl' % (Dg, DF))
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    keys = sorted(src)
    if len(sys.argv) > 3 and sys.argv[3]:
        ms = [float(a) for a in sys.argv[3].split(',')]; keys = [k for k in keys if k[0] in ms]
    if len(sys.argv) > 4 and sys.argv[4]:
        js = [int(a) for a in sys.argv[4].split(',')]; keys = [k for k in keys if k[1] in js]
    for key in keys:
        if key in db:
            continue
        r0 = src[key]
        mfl, J = key
        t0 = time.time()
        cm = None
        for Ntry in ([None] if mfl < 1.1 else [2000, int(0.06*r0['Delta']**2), int(0.15*r0['Delta']**2), None]):
            try:
                cm = CertMoments(r0['Delta'], J, DF, Lamf, ks=KE, N=Ntry)   # crude M suffices away from the gap
                break
            except AssertionError:
                continue
        D = cm.D
        al, alb = D - J - 2, D + J - 1
        m2 = al*alb
        u = Lam*Lam/m2
        eta = (al*al + alb*alb)/(2*m2)
        beta = (eta + (eta*eta - 1).sqrt()).log()
        P0 = arb(eta).legendre_p(arb(-0.5), 0, 3) if hasattr(arb, 'legendre_p') else None
        P0 = acb(eta).legendre_p(acb(-0.5), 0, 3).real
        pmax = 1/u.sqrt()                                   # singularity at p = m/Lam
        best = None
        exc = {}
        for frac in ((0.4, 0.6, 0.8, 0.9) if mfl < 1.3 else (0.3,)):
            R = 1 + frac*(pmax - 1)                          # right vertex of the ellipse
            s = 4*R - 2                                       # rho + 1/rho
            rho = (s + (s*s - 4).sqrt())/2
            Kc = 720
            Mb = {}
            for c in range(Kc):
                th = 2*np.pi*(c + 0.5)/Kc
                dth = np.pi/Kc
                ctr = acb(0.5) + (acb(rho)*acb(np.cos(th), np.sin(th)) + acb(np.cos(th), -np.sin(th))/rho)/4
                rad = (rho + 1/rho)/4*arb(dth)*1.01
                P = ctr + acb(arb(0, rad.upper()), arb(0, rad.upper()))
                y = u*P*P
                nu = Lam*P
                numax = (Lam*(abs(P).upper()))
                imnu = abs(nu.imag).upper()
                nm = float(np.ceil(float(numax.upper())*4)/4)
                exb = {}
                for k in KE:
                    if (k, nm) not in exc:
                        exc[(k, nm)] = ex_abs_bound(cm, k, nm)
                    exb[k] = exc[(k, nm)]
                x = P*u.sqrt()
                Pc_b = (arb(imnu)*beta).cosh()*P0 + (abs(((2*J + 1)*x).imag).upper()*arb(1)).cosh()
                PJ = jacobi_acb(J, 0, 0, 1 - 2*y)[J]
                d44 = (1 - y)**4*jacobi_acb(J - 4, 0, 8, 1 - 2*y)[J - 4] if J >= 4 else acb(0)
                z1 = abs(1 - y); z2 = abs(2 - y)
                zlo, zhi = arb(z1.lower()), arb(z1.upper()); wlo, whi = arb(z2.lower()), arb(z2.upper())
                if not (zlo > 0):
                    zlo = arb(1e-300)
                Eup = {k: whi/zlo**(arb(2 + k)/2) for k in range(2, 11)}     # real-modulus bounds (robust for wide balls)
                Elo = {k: wlo/zhi**(arb(2 + k)/2) for k in range(2, 11)}
                for h in SIGMA:
                    if (h == '++' and J % 2) or (h == '+-' and J < 4):
                        continue
                    Dh = PJ if h == '++' else d44
                    tot = arb(0); dbg = {}
                    for k in KS:
                        if k % 2 == 0:
                            Mod = exb[k] + Eup[k]*abs(Dh - PJ)
                        else:
                            Lc = {l: abs(arb(float(np.prod([(k - lb)/(l - lb) for lb in KE if lb != l])))) for l in KE}
                            exk = Eup[k]*sum((Lc[l]*exb[l]/Elo[l] for l in KE), arb(0))
                            wk = zlo**(-arb(k + 1)/2)
                            Mod = wk*abs(Dh - PJ) + Pc_b + wk*(exk + 2*Pc_b)/Elo[k]
                        dbg[k] = (Mod, abs(phi_acb(k, P)))
                        tot += u**(k - 2)*abs(phi_acb(k, P))*Mod
                    v = float(tot.upper())
                    if not np.isfinite(v):
                        if os.environ.get('CQDBG'):
                            print('nonfinite', h, c, P, dbg, flush=True)
                        v = float('inf')
                    Mb[h] = max(Mb.get(h, 0.0), v)
            errs = {h: (arb(0.5)*arb(64)/15*arb(Mb[h])*rho**(-2*NQ)/(rho*rho - 1) if np.isfinite(Mb[h]) else arb('inf')) for h in Mb}
            cand = (max(float(e.upper()) for e in errs.values()), frac, float(rho.mid()), {h: float(e.upper()) for h, e in errs.items()},
                    dict(Mb))
            if best is None or cand[0] < best[0]:
                best = cand
        res = {}
        for h, v in r0['res'].items():
            e = best[3][h]
            res[h] = dict(F_lo=v['F_lo'], err=e, Fexact_lo=v['F_lo'] - e, Ns=v['Ns'], cert=(v['F_lo'] - e) > 0)
        db[key] = dict(res=res, rho=best[2], frac=best[1], M=best[4], time=time.time() - t0)
        msg = " ".join(f"{h}: F_lo/Ns {v['F_lo']/v['Ns']:+.4f} quad.err/Ns {v['err']/v['Ns']:.1e} -> {'CERT' if v['cert'] else 'FAIL'}" for h, v in res.items())
        print(f"Dg={Dg} Df={DF:g} m/Lam={mfl:.2f} J={J:2d} rho={best[2]:.3f}: {msg}  ({time.time()-t0:.0f}s)", flush=True)
        pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print("# done")


if __name__ == '__main__':
    run()
