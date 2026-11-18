"""
cert_action: rigorous (Arb) enclosure of the D-safe heavy action F (NF02 uplift Omega^{NF,SC}, OMEGA_NF.md; no forward terms) on heavy states
(m, J) of the scalar-skeleton model, both channels, with the 64-node Gauss-Legendre p-rule of OMEGA_NF as the *definition* of the p-integral:

  F_h(m, J) = sum_{k=2..7} u^{k-2} sum_i w_i phi_k(p_i) Mod_k^h(p_i),   u = Lam^2/m^2,  nu_i = Lam p_i,  x_i = p_i sqrt(u),  y = x^2,
  even k:  Mod = ex_k + E_k (D_h - P_J)
  odd k:   Mod = sigma_h [ w_k (D_h - P_J) + Pc + w_k (ex_k - 2 Pc)/E_k ],   w_k = (1-y)^{-(k+1)/2},  E_k = (2-y)/(1-y)^{1+k/2},
  Pc = P_{-1/2+i nu}(cosh beta) - J_0((2J+1) x),  D_++ = P_J(1-2y),  D_+- = d^J_44 = (1-y)^4 P^{(0,8)}_{J-4}(1-2y),  sigma_++ = 1, sigma_+- = -1
  (this is exactly OMEGA_NF's model: e_k D_h + 2 eps_k sigma^k Pc + sigma^k env_k rho_k with rho_k = (ex_k - E_k P_J - 2 Pc)/E_k),
  ex_k = m^{2k} C_{k, nu}[Delta, J] for k = 2,4,6,8,10 enclosed by cert_moments.py (rigorous n- and m-tails),
  odd k: ex_k := E_k * (quartic Lagrange interpolant in k of ex_l/E_l, l = 2,4,6,8,10)   [a definition of the odd moments, M].
Everything (nodes as given doubles, exact NF02 rationals, Legendre/Jacobi recurrences, P_{-1/2+i nu} via Arb legendre_p, J_0 via Arb) in Arb.
Output: lower bound of F (certified > 0 or not), F/Ns (Ns in floating point, for reporting only).
Usage: python3 cert_action.py Dg Df [mfl list comma] [J list comma]  -> cache/cert_<Dg>_Df<Df>.pkl (incremental, resumable), .log
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
Dg = int(sys.argv[1]); DF = float(sys.argv[2])
os.environ['OMEGAU_DF'] = '%g' % DF
import numpy as np
from fractions import Fraction as Fr
from flint import arb, acb, ctx
from nf_common import NF, TransformsNF, flat_parts_nf, action_nf, KS, SIGMA
from common import state, gl01
from cert_moments import CertMoments, af

MGRID = [1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.12, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0]
JS = [0, 2, 4, 5, 6, 7, 8, 9, 10]
if len(sys.argv) > 3 and sys.argv[3]:
    MGRID = [float(a) for a in sys.argv[3].split(',')]
if len(sys.argv) > 4 and sys.argv[4]:
    JS = [int(a) for a in sys.argv[4].split(',')]
KE = (2, 4, 6, 8, 10)
NP = 64
pnf, pwf = gl01(NP)
F = NF(); TRf = TransformsNF(F)


def jacobi_arb(nmax, a, b, x):
    P = [arb(1)]
    if nmax >= 1:
        P.append((a + 1) + (a + b + 2)*(x - 1)/2)
    for n in range(2, nmax + 1):
        c = 2*n + a + b
        P.append((((c - 1)*(a*a - b*b) + (c - 2)*(c - 1)*c*x)*P[n - 1] - 2*(n + a - 1)*(n + b - 1)*c*P[n - 2])/(2*n*(n + a + b)*(c - 2)))
    return P


def phi_arb(k, p):
    return sum((af(c)*p**i for i, c in enumerate(F.PHx[k])), arb(0))


def run():
    ctx.prec = 192
    Lam = arb(98*99).sqrt() if Dg == 100 else arb((Dg - 2)*(Dg - 1)).sqrt()
    Lamf = float(np.sqrt((Dg - 2.0)*(Dg - 1.0)))
    pk = os.path.join(HERE, 'cache', 'cert_%d_Df%g.pkl' % (Dg, DF))
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    for mfl in MGRID:
        for J in JS:
            key = (mfl, J)
            if key in db and all(v['cert'] for v in db[key]['res'].values()):
                continue
            st = state(Dg, Lamf, mfl, J)
            if st is None:
                continue
            t0 = time.time()
            Dq = st['Delta']
            cm = CertMoments(Dq, J, DF, Lamf, ks=KE)
            D = cm.D
            al, alb = D - J - 2, D + J - 1
            m2 = al*alb
            u = Lam*Lam/m2
            eta = (al*al + alb*alb)/(2*m2)
            ps = [arb(float(v)) for v in pnf]; ws = [arb(float(v)) for v in pwf]
            nus = [Lam*p for p in ps]
            ex = {k: cm.values(k, [float(v) for v in Lamf*pnf]) for k in KE}
            # nodes nu_i = Lam p_i used in cert_moments as floats Lamf*p_i: re-define nu_i consistently as those floats
            nus = [arb(float(v)) for v in Lamf*pnf]
            xs = [nu*u.sqrt()/Lam for nu in nus]
            res = {}
            for h in SIGMA:
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                s = SIGMA[h]
                Ftot = arb(0)
                parts = {}
                for k in KS:
                    pk_ = arb(0)
                    for i in range(NP):
                        x = xs[i]; y = x*x; p = nus[i]/Lam
                        Ek = lambda kk: (2 - y)/(1 - y)**(1 + arb(kk)/2)
                        PJ = jacobi_arb(J, 0, 0, 1 - 2*y)[J]
                        Dh = PJ if h == '++' else (1 - y)**4*jacobi_arb(J - 4, 0, 8, 1 - 2*y)[J - 4]
                        if k % 2 == 0:
                            Mod = ex[k][i] + Ek(k)*(Dh - PJ)
                        else:
                            r = sum((ex[l][i]/Ek(l)*np.prod([(k - lb)/(l - lb) for lb in KE if lb != l]) for l in KE), arb(0))
                            exk = Ek(k)*r
                            Pnu = acb(eta).legendre_p(acb(arb(-0.5), nus[i]), 0, 3).real
                            Pc = Pnu - ((2*J + 1)*x).bessel_j(0)
                            wk = (1 - y)**(-(arb(k) + 1)/2)
                            Mod = s*(wk*(Dh - PJ) + Pc + wk*(exk - 2*Pc)/Ek(k))
                        pk_ += ws[i]*phi_arb(k, p)*Mod
                    parts[k] = u**(k - 2)*pk_
                    Ftot += parts[k]
                # Ns (floating point, reporting only) from the OMEGA_NF action
                uf = 1/mfl**2
                pf, wf = gl01(max(400, 2*J + 100))
                Js, B2, B3, fparts = flat_parts_nf(F, h, uf, J, pf, wf)
                ii = int(np.searchsorted(Js, J))
                skz = {k: (np.zeros(NP), np.zeros(NP), np.ones(NP, bool)) for k in KS}
                stf = dict(m=st['m'], m2=st['m2'], eta=st['eta'], beta=st['beta'], Delta=st['Delta'])
                rr = action_nf(F, TRf, h, stf, J, mfl, Lamf, pnf, pwf, skz, fparts[ii], B2[ii], B3[ii])
                lo = float(Ftot.lower()); hi = float(Ftot.upper())
                res[h] = dict(F_lo=lo, F_hi=hi, F_mid=float(Ftot.mid()), rad=float(Ftot.rad()), Ns=rr['Ns'], cert=lo > 0,
                              parts={k: (float(v.mid()), float(v.rad())) for k, v in parts.items()})
            db[key] = dict(Delta=Dq, J=J, mfl=mfl, res=res, N=cm.N, M=cm.M, P={k: cm.info[k]['P'] for k in KE}, time=time.time() - t0)
            msg = " ".join(f"{h}: F/Ns in [{v['F_lo']/v['Ns']:+.5f}, {v['F_hi']/v['Ns']:+.5f}] {'CERT' if v['cert'] else 'FAIL'}"
                           for h, v in res.items())
            print(f"Dg={Dg} Df={DF:g} m/Lam={mfl:.2f} J={J:2d} Delta={Dq:.4f} N={cm.N} M={cm.M}: {msg}  ({time.time()-t0:.0f}s)", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print("# done")


if __name__ == '__main__':
    run()
