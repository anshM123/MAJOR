"""
cert_actionA: certified near-gap positivity of the ROUTE-A' D-safe action (H0_SPINNING_FUNCTIONAL.md s10/R; model of spinfunc/f6_routeA2.py,
stripped envelope, no eta-inversion -- see CERTIFIED_NEARGAP_A.md for the inversion) in the scalar-skeleton model, exact p-integral.

Per state (Delta, J), u = Lam^2/m^2, nu = Lam p, x = p sqrt(u), y = x^2, sigma = +-1:
  F = int_0^1 [ psi2(p) Mod2(p) + u sigma psi3(p) Mod3(p) ] dp,
  Mod2 = (2-y) Dt + 2 Pc + fac (ex2 - E2 P_J - 2 Pc),
  Mod3 = Dt + Pc + (1-y)^{-2} fac ( sum_l L_l(3) ex_l/E_l - P_J - 2 Pc/E3 ),
  Dt = P_J(1-2y) (++), P^{(0,8)}_{J-4}(1-2y) (+-) [= d^J_44/(1-y)^4];  fac = (1-y)^2 (++), (1-y)^{-2} (+-);
  E_k = (2-y)/(1-y)^{1+k/2};  Pc = P_{-1/2+i nu}(cosh beta) - J0((2J+1)x);  ex_l = m^{2l} C_{l,nu}, l = 2,4,...,10 (cert_moments.py);
  L_l(3) the quartic Lagrange weights in k (odd-k definition [M], as cert_action.py).
This is exactly f6's  B2 + B3 + psi2 2(P-J0) + psi3 u sigma (P-J0) + psi2 E2 fac rho2 + u sigma psi3 env3 fac rho3  written per node.
Step 1: F^(64) = 64-node Gauss sum, every quantity an Arb ball (192 bits).  Step 2: certified Gauss error (Trefethen ATAP Thm 19.3) from
|f| <= M on a Bernstein ellipse with right vertex below p = m/Lam, M bounded on 720 complex balls exactly as cert_quad.py (ex_abs_bound,
real-modulus bounds for (1-y), (2-y) powers, Mehler-Dirichlet bound on P, cosh bound on J0).  Certified: F_exact >= F64_lo - err > 0.
Usage: python3 cert_actionA.py Dg Df [mfl list] [J list]  -> cache/certA_<Dg>_Df<Df>.pkl (incremental, resumable), stdout log
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
Dg = int(sys.argv[1]); DF = float(sys.argv[2])
os.environ['OMEGAU_DF'] = '%g' % DF
import numpy as np
from flint import arb, acb, ctx
from nf_common import NF, SIGMA
from common import state, gl01
from cert_moments import CertMoments, af
import cert_quad as CQ                       # cpow, jacobi_acb, ex_abs_bound (read-only reuse)
from h13_common import flat_rows
from h13a_common import hybridA

MGRID = [1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.12, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0]
JS = [0, 2, 4, 5, 6, 7, 8, 9, 10]
if len(sys.argv) > 3 and sys.argv[3]:
    MGRID = [float(a) for a in sys.argv[3].split(',')]
if len(sys.argv) > 4 and sys.argv[4]:
    JS = [int(a) for a in sys.argv[4].split(',')]
KE = (2, 4, 6, 8, 10)
NQ = 64
pnf, pwf = gl01(NQ)
F = NF()
L3 = {l: float(np.prod([(3 - lb)/(l - lb) for lb in KE if lb != l])) for l in KE}
C2 = [af(c) for c in F.s2]; C3 = [af(c) for c in F.s3]


def horner(c, p, z):
    r = z
    for a in reversed(c):
        r = r*p + a
    return r


def jac(nmax, a, b, x):
    P = [x*0 + 1]
    if nmax >= 1:
        P.append((a + 1) + (a + b + 2)*(x - 1)/2)
    for n in range(2, nmax + 1):
        c = 2*n + a + b
        P.append((((c - 1)*(a*a - b*b) + (c - 2)*(c - 1)*c*x)*P[n - 1] - 2*(n + a - 1)*(n + b - 1)*c*P[n - 2])/(2*n*(n + a + b)*(c - 2)))
    return P


def action64(cm, J, Lam, u, eta, Lamf):
    nus = [arb(float(v)) for v in Lamf*pnf]
    ws = [arb(float(v)) for v in pwf]
    ex = {k: cm.values(k, [float(v) for v in Lamf*pnf]) for k in KE}
    out = {}
    for h, s in SIGMA.items():
        if (h == '++' and J % 2) or (h == '+-' and J < 4):
            continue
        tot = arb(0)
        for i in range(NQ):
            p = nus[i]/Lam; x = p*u.sqrt(); y = x*x
            PJ = jac(J, 0, 0, 1 - 2*y)[J]
            Dt = PJ if h == '++' else jac(J - 4, 0, 8, 1 - 2*y)[J - 4]
            fac = (1 - y)**2 if h == '++' else 1/(1 - y)**2
            E = {k: (2 - y)/(1 - y)**(1 + arb(k)/2) for k in (2, 3, 4, 6, 8, 10)}
            Pnu = acb(eta).legendre_p(acb(arb(-0.5), nus[i]), 0, 3).real
            Pc = Pnu - ((2*J + 1)*x).bessel_j(0)
            r3 = sum((arb(L3[l])*ex[l][i]/E[l] for l in KE), arb(0))
            Mod2 = (2 - y)*Dt + 2*Pc + fac*(ex[2][i] - E[2]*PJ - 2*Pc)
            Mod3 = Dt + Pc + fac/(1 - y)**2*(r3 - PJ - 2*Pc/E[3])
            tot += ws[i]*(horner(C2, p, arb(0))*Mod2 + u*s*horner(C3, p, arb(0))*Mod3)
        out[h] = tot
    return out


def quad_error(cm, J, Lam, u, eta, mfl):
    beta = (eta + (eta*eta - 1).sqrt()).log()
    P0 = acb(eta).legendre_p(acb(-0.5), 0, 3).real
    pmax = 1/u.sqrt()
    best = None; exc = {}
    for frac in ((0.4, 0.6, 0.8, 0.9) if mfl < 1.3 else (0.3,)):
        R = 1 + frac*(pmax - 1); sE = 4*R - 2; rho = (sE + (sE*sE - 4).sqrt())/2
        Kc = 720; Mb = {}
        for c in range(Kc):
            th = 2*np.pi*(c + 0.5)/Kc; dth = np.pi/Kc
            ctr = acb(0.5) + (acb(rho)*acb(np.cos(th), np.sin(th)) + acb(np.cos(th), -np.sin(th))/rho)/4
            rad = (rho + 1/rho)/4*arb(dth)*1.01
            P = ctr + acb(arb(0, rad.upper()), arb(0, rad.upper()))
            y = u*P*P; nu = Lam*P
            numax = Lam*(abs(P).upper()); imnu = abs(nu.imag).upper()
            nm = float(np.ceil(float(numax.upper())*4)/4)
            exb = {}
            for k in KE:
                if (k, nm) not in exc:
                    exc[(k, nm)] = CQ.ex_abs_bound(cm, k, nm)
                exb[k] = exc[(k, nm)]
            x = P*u.sqrt()
            Pc_b = (arb(imnu)*beta).cosh()*P0 + (abs(((2*J + 1)*x).imag).upper()*arb(1)).cosh()
            PJ = abs(CQ.jacobi_acb(J, 0, 0, 1 - 2*y)[J])
            z1 = abs(1 - y); z2 = abs(2 - y)
            zlo, zhi = arb(z1.lower()), arb(z1.upper()); wlo, whi = arb(z2.lower()), arb(z2.upper())
            if not (zlo > 0):
                zlo = arb(1e-300)
            Eup2 = whi/zlo**2
            Elo = {k: wlo/zhi**(1 + arb(k)/2) for k in (3, 4, 6, 8, 10, 2)}
            a2 = abs(horner(C2, P, acb(0))); a3 = abs(horner(C3, P, acb(0)))
            for h in SIGMA:
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                Dt = PJ if h == '++' else abs(CQ.jacobi_acb(J - 4, 0, 8, 1 - 2*y)[J - 4])
                fac = zhi**2 if h == '++' else 1/zlo**2
                M2 = z2.upper()*Dt + 2*Pc_b + fac*(exb[2] + Eup2*PJ + 2*Pc_b)
                r3 = sum((abs(arb(L3[l]))*exb[l]/Elo[l] for l in KE), arb(0))
                M3 = Dt + Pc_b + fac/zlo**2*(r3 + PJ + 2*Pc_b/Elo[3])
                v = float((a2*M2 + u*a3*M3).upper())
                if not np.isfinite(v):
                    v = float('inf')
                Mb[h] = max(Mb.get(h, 0.0), v)
        errs = {h: (float((arb(0.5)*arb(64)/15*arb(Mb[h])*rho**(-2*NQ)/(rho*rho - 1)).upper()) if np.isfinite(Mb[h]) else float('inf'))
                for h in Mb}
        cand = (max(errs.values()), frac, float(rho.mid()), errs)
        if best is None or cand[0] < best[0]:
            best = cand
    return best


def run():
    ctx.prec = 192
    Lamf = float(np.sqrt((Dg - 2.0)*(Dg - 1.0)))
    Lam = arb((Dg - 2)*(Dg - 1)).sqrt()
    pk = os.path.join(HERE, 'cache', 'certA_%d_Df%g.pkl' % (Dg, DF))
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    for mfl in MGRID:
        for J in JS:
            key = (mfl, J)
            if key in db:
                continue
            st = state(Dg, Lamf, mfl, J)
            if st is None:
                continue
            t0 = time.time()
            cm = CertMoments(st['Delta'], J, DF, Lamf, ks=KE)
            D = cm.D; al, alb = D - J - 2, D + J - 1; m2 = al*alb
            u = Lam*Lam/m2; eta = (al*al + alb*alb)/(2*m2)
            F64 = action64(cm, J, Lam, u, eta, Lamf)
            t1 = time.time() - t0
            qb = quad_error(cm, J, Lam, u, eta, mfl)
            res = {}
            for h, v in F64.items():
                B2, B3 = flat_rows(h, 1/mfl**2, [J])[J]
                hy = hybridA(h, Dg, mfl, J, B2, B3, beta=st['beta'])
                lo = float(v.lower()); e = qb[3][h]
                res[h] = dict(F_lo=lo, F_hi=float(v.upper()), F_mid=float(v.mid()), err=e, Fex_lo=lo - e, Ns=hy['Ns'],
                              cert64=lo > 0, cert=(lo - e) > 0)
            db[key] = dict(Delta=st['Delta'], J=J, mfl=mfl, res=res, rho=qb[2], frac=qb[1], N=cm.N, M=cm.M, time=time.time() - t0, t64=t1)
            msg = " ".join(f"{h}: F/Ns in [{v['F_lo']/v['Ns']:+.5f}, {v['F_hi']/v['Ns']:+.5f}] quad.err/Ns {v['err']/v['Ns']:.1e} "
                           f"{'CERT' if v['cert'] else ('CERT64' if v['cert64'] else 'FAIL')}" for h, v in res.items())
            print(f"Dg={Dg} Df={DF:g} m/Lam={mfl:.2f} J={J:2d} rho={qb[2]:.3f}: {msg}  ({t1:.0f}+{time.time()-t0-t1:.0f}s)", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print("# done", flush=True)


if __name__ == '__main__':
    run()
