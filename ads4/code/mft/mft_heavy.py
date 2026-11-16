"""
Task B: large-Delta (heavy) limit of the MFT OPE coefficients and the MFT-normalised heavy formula.

Inputs:  closed forms of mft_exact.py (verified exactly against the Fock-space computation);
         SPINNING_E2E.md sect. 4.5 heavy formula (predict.py), Lorentzian shadow coefficient b^{(D34)} (formulas.py).

Checks performed (all numbers printed):
 (1) exact identity  P_FK(D1,D2)/b^{(D12)} = (2/pi) A+ A- m^2 G(alpha) G(alphabar) / [Gamma(D1-1/2)Gamma(D2-1/2)Gamma(D1)Gamma(D2)],
     G(x) = Gamma((x+D1+D2-2)/2)/Gamma((x+4-D1-D2)/2)                                              (random points)
 (2) scalar MFT normalisation:  P_FK (pi/4) 4^{D1+D2}/b * pi m^2/(2 n_J (2Delta-3)) * m^{-2(D1+D2-1)} * a_{D1,D2} -> 1,
     a_{D1,D2} = 4 Gamma(D1)Gamma(D1-1/2)Gamma(D2)Gamma(D2-1/2)        [= CHMRS a_{Dphi} for D1 = D2]
 (3) spinning, parity even (l = 1, 2): with q^{(l)} = (1/(pi p^MFT_E)) (2 n_J/m^2)(2Delta-3) [CHMRS (2.32) with p^MFT_E]:
       (a_{Delta_O,Dphi}/(gamma_ax gamma_ay)) m^{2k} * [prediction per lambda^2 of SPINNING_E2E 4.5] / q^{(l)}
          =  Q_E(Delta,J) * (1/4) [M_pres P + M_flip F^nu_{2l}],
       Q_E = Gamma-ratio(1,Dphi) * (A+A-)^{2l} / prod_{i<=l} (A+^2-(i-1/2)^2)(A-^2-(i-1/2)^2)  -> 1 + O(1/m^2)
     (evaluated with predict.prediction itself, i.e. with the code that was verified end-to-end)
 (4) parity odd, l = 1:  Q_O = Gamma-ratio(2,Dphi);  vertex constants of the odd structure (heavy frame and infinity frame)
 (5) S-part symmetry of the M table: S^{q1q3}_{-h1,-h3} = -S^{q1q3}_{h1h3}  (so S cancels for both parities)
 (6) end-to-end: MFT-normalised measured moments (out_t9_tables.txt of e2e) / (1/4)[M_pres P + M_flip F]
"""
import os, sys, re
import mpmath as mp, sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
E2E = os.path.join(here, '..', 'e2e')
sys.path.insert(0, here); sys.path.insert(0, E2E)
from formulas import log_b, m2_eta, gamma_a, Fnu_mp
from predict import prediction, components
mp.mp.dps = 30
half = mp.mpf(1)/2

def logPFK(D1, D2, Delta, J):
    n = (Delta - D1 - D2 - J)/2; h = mp.mpf(3)/2; g = mp.loggamma
    lrf = lambda a, k: g(a + k) - g(a)
    return (lrf(D1-h+1, n) + lrf(D2-h+1, n) + lrf(D1, J+n) + lrf(D2, J+n)
            - (g(n+1) + g(J+1) + lrf(J+h, n) + lrf(D1+D2+n-2*h+1, n) + lrf(D1+D2+2*n+J-1, J) + lrf(D1+D2+n+J-h, n)))

def a12(D1, D2):
    return 4*mp.gamma(D1)*mp.gamma(D1-half)*mp.gamma(D2)*mp.gamma(D2-half)

def Aplusminus(Delta, J):
    return mp.mpf(Delta) - mp.mpf(3)/2, mp.mpf(J) + half

def gamma_ratio(D1, D2, Delta, J):
    """(4/m^2)^{S-3} G(alpha)G(alphabar), S = D1 + D2: the exact finite-m factor of the scalar MFT normalisation"""
    al, alb = mp.mpf(Delta) - J - 2, mp.mpf(Delta) + J - 1
    S = mp.mpf(D1) + D2
    G = lambda x: mp.exp(mp.loggamma((x + S - 2)/2) - mp.loggamma((x + 4 - S)/2))
    return (4/(al*alb))**(S - 3)*G(al)*G(alb)

def pol_ratio(l, Delta, J):
    Ap, Am = Aplusminus(Delta, J)
    num = (Ap*Am)**(2*l); den = 1
    for i in range(1, l+1):
        den *= (Ap**2 - (i-half)**2)*(Am**2 - (i-half)**2)
    return num/den

C_L = {1: 4, 2: 24}
def log_pMFT_L(l, Dphi, Delta, J):
    """log of p^MFT for the structure L = sum L_n S_n  (= L_0 S_cons)"""
    return (J*mp.log(2) + logPFK(1, Dphi, Delta, J) - mp.log(C_L[l]) - mp.log(mp.rf(J-l+1, 2*l)) - mp.log(mp.rf(Delta-l-1, 2*l)))

def L0(l, Dphi, Delta, J):
    if l == 1: return -J*(Delta-1)*(Delta-Dphi-J-1)
    return Delta*J*(Delta-1)*(J-1)*(Delta-Dphi-J-3)*(Delta-Dphi-J-1)/2

def q_factor_inv(Delta, J, logp):
    """1/q = pi p^MFT m^2/(2 n_J (2 Delta - 3)),  n_J = 16 pi (2J+1)  (d = 3)"""
    m2, _ = m2_eta(mp.mpf(Delta), J)
    nJ = 16*mp.pi*(2*J + 1)
    return mp.exp(logp)*mp.pi*m2/(2*nJ*(2*mp.mpf(Delta) - 3))

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def check1():
    pr("## (1) exact identity P_FK/b  (ratio lhs/rhs at random points)")
    import random
    random.seed(3)
    for _ in range(5):
        D1 = mp.mpf(random.randint(5, 40))/10; D2 = mp.mpf(random.randint(5, 40))/10
        n = random.randint(0, 30); J = random.randint(0, 40)
        Delta = D1 + D2 + 2*n + J
        al, alb = Delta - J - 2, Delta + J - 1
        Ap, Am = Aplusminus(Delta, J)
        S = D1 + D2
        G = lambda x: mp.gamma((x + S - 2)/2)/mp.gamma((x + 4 - S)/2)
        rhs = (2/mp.pi)*Ap*Am*al*alb*G(al)*G(alb)/(mp.gamma(D1-half)*mp.gamma(D2-half)*mp.gamma(D1)*mp.gamma(D2))
        lhs = mp.exp(logPFK(D1, D2, Delta, J) - log_b(Delta, J, D1 - D2))
        pr("  D1=%s D2=%s n=%d J=%d   lhs/rhs - 1 = %s" % (D1, D2, n, J, mp.nstr(lhs/rhs - 1, 3)))

def check2():
    pr("## (2) scalar MFT normalisation  R_s a_{D1,D2} (-> 1, = gamma_ratio exactly)")
    for (D1, D2) in [(1, mp.mpf(7)/4), (2, mp.mpf(7)/4), (3, mp.mpf(7)/4), (mp.mpf(7)/4, mp.mpf(7)/4), (2, mp.mpf(13)/10)]:
        for (n, J) in [(10, 5), (100, 50), (1000, 500), (1000, 4000)]:
            Delta = D1 + D2 + 2*n + J
            m2, _ = m2_eta(Delta, J)
            lg = (logPFK(D1, D2, Delta, J) + mp.log(mp.pi/4) + (D1 + D2)*mp.log(4) - log_b(Delta, J, D1 - D2)
                  - (D1 + D2 - 1)*mp.log(m2))
            val = mp.exp(lg)*q_factor_inv(Delta, J, 0)*a12(D1, D2)
            pr("  (D1,D2)=(%s,%s) n=%d J=%d  m=%.1f:  R_s a = %s   gamma_ratio = %s   m^2(R_s a - 1) = %s" % (
                mp.nstr(D1, 4), mp.nstr(D2, 4), n, J, float(mp.sqrt(m2)), mp.nstr(val, 12), mp.nstr(gamma_ratio(D1, D2, Delta, J), 12),
                mp.nstr(m2*(val - 1), 6)))

def check3(Dphi=mp.mpf(7)/4):
    pr("## (3) parity even: MFT-normalised leading formula / [(1/4)(M_pres P + M_flip F)] = Q_E  (predict.py evaluated)")
    pr("| l | Delta | J | m | beta | k | nu | mode | normalised/((1/4)(MP+MF)) | Q_E (analytic) | ratio | m^2(Q_E-1) | kappa_l(1/sh^2+1/ch^2) - (8/3)B3((Dphi-1)/2) cosh b |")
    pr("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    import numpy as np
    for l in (1, 2):
        DO = l + 1
        kap = {1: mp.mpf(1)/4, 2: mp.mpf(5)/2}[l]
        for beta in (0.4, 1.0, 2.0):
            for m in (20, 100, 1000):
                J = int(round((2*m*np.sinh(beta/2) - 1)/2)); mm = (2*J + 1)/(2*np.sinh(beta/2))
                Delta = mp.mpf((mm*np.exp(-beta/2) + mm*np.exp(beta/2) + 3)/2)
                m2, eta = m2_eta(Delta, J); bb = mp.acosh(eta)
                for k in (4, 6):
                    for nu in (0.0, 3.0):
                        q1, q3 = l, -l
                        pd = (-1)**l*prediction(l, Delta, J, Dphi, k, nu, q1, q3)        # (-1)^{J+l} Pihat per W_diff ...
                        # per lambda^2 of S_cons:  divide by 2^J L_0^2 ; times 1/q with p^MFT_E = p_L L_0^2
                        lp = log_pMFT_L(l, Dphi, Delta, J)
                        val = pd/(2**J)*q_factor_inv(Delta, J, lp)              # L_0^2 cancels exactly
                        ax, ay = 2*Dphi + k - 1, 2*DO + k - 1
                        val *= a12(DO, Dphi)/(gamma_a(ax, nu)*gamma_a(ay, nu))*m2**k
                        Mp, Mf = components(l, k, nu, q1, q3)
                        ref = (Mp*Fnu_mp(nu, 0, bb) + Mf*Fnu_mp(nu, 2*l, bb))/4
                        QE = gamma_ratio(1, Dphi, Delta, J)*pol_ratio(l, Delta, J)
                        s = (Dphi - 1)/2; B3 = s*(s - half)*(s - 1)
                        pred2 = kap*(1/mp.sinh(bb/2)**2 + 1/mp.cosh(bb/2)**2) - mp.mpf(8)/3*B3*mp.cosh(bb)
                        r = mp.re(val/ref)
                        pr("| %d | %.4f | %d | %.1f | %.3f | %d | %.0f | (%+d,%+d) | %s | %s | %s | %s | %s |" % (
                            l, float(Delta), J, float(mp.sqrt(m2)), float(bb), k, nu, q1, q3, mp.nstr(r, 12), mp.nstr(QE, 12),
                            mp.nstr(r/QE, 15), mp.nstr(m2*(QE - 1), 6), mp.nstr(pred2, 6)))

def check4(Dphi=mp.mpf(7)/4):
    pr("## (4) parity odd, l = 1: exact finite-m factor Q_O = Gamma-ratio(2, Dphi) (-> 1 + O(1/m^2)); no normalisation polynomial")
    import numpy as np
    for beta in (0.4, 1.0, 2.0):
        for m in (20, 100, 1000):
            J = int(round((2*m*np.sinh(beta/2) - 1)/2)); mm = (2*J + 1)/(2*np.sinh(beta/2))
            Delta = mp.mpf((mm*np.exp(-beta/2) + mm*np.exp(beta/2) + 3)/2)
            m2, eta = m2_eta(Delta, J)
            QO = gamma_ratio(2, Dphi, Delta, J)
            s = Dphi/2; B3 = s*(s - half)*(s - 1)
            pr("  beta=%.1f m=%.1f: Q_O = %s   m^2(Q_O-1) = %s   predicted -(8/3)B3(Dphi/2) cosh beta = %s" % (
                beta, float(mp.sqrt(m2)), mp.nstr(QO, 12), mp.nstr(m2*(QO - 1), 6), mp.nstr(-mp.mpf(8)/3*B3*eta, 6)))

def check5():
    pr("## (5) symmetry of the S-part of the M table (laplace.py):  max |S^{q1q3}_{h1h3} + S^{q1q3}_{-h1,-h3}|")
    from laplace import M_table
    for l in (1, 2):
        T = M_table(l)
        worst = 0
        for (q1, q3, h1, h3), v in T.items():
            d = sp.simplify(v[1] + T[(q1, q3, -h1, -h3)][1])
            if d != 0: worst = d
        pr("  l = %d: %s  (0 means: S cancels in sum_h and in sum_h sgn(h1)sgn(h3))" % (l, worst))

def check6(Dphi=mp.mpf(7)/4):
    pr("## (6) end-to-end: exact moments of SPINNING_E2E (out_t9_tables.txt), MFT-normalised, / (1/4)[M_pres P + M_flip F]")
    pr("##     = (meas/pred of e2e) x Q_E;  helicity components: R^MFT = R_e2e x Q_E")
    pr("| l | Delta | J | m | beta | k | nu | R_pres(e2e) | R_pres^MFT | R_flip(e2e) | R_flip^MFT | m^2(R_pres-1) | m^2(R_pres^MFT-1) |")
    pr("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    l = None
    for line in open(os.path.join(E2E, 'out_t9_tables.txt')).read().splitlines():
        if line.startswith('## l ='): l = int(line.split('=')[1].split(':')[0]); continue
        if not line.startswith('| ') or line.startswith('| Delta'): continue
        f = [x.strip() for x in line.strip('|').split('|')]
        Delta, J, m, beta, k, nu, mode = float(f[0]), int(f[1]), float(f[2]), float(f[3]), int(f[4]), float(f[5]), f[6]
        if mode != '(+1,-1)' and mode != '(+2,-2)': continue
        if k != 6: continue
        Rp, Rf = mp.mpf(f[11]), mp.mpf(f[12])
        Q = gamma_ratio(1, Dphi, mp.mpf(Delta), J)*pol_ratio(l, mp.mpf(Delta), J)
        mm = mp.mpf(m)
        pr("| %d | %.4f | %d | %.1f | %.2f | %d | %.0f | %s | %s | %s | %s | %s | %s |" % (
            l, Delta, J, m, beta, k, nu, mp.nstr(Rp, 7), mp.nstr(Rp*Q, 7), mp.nstr(Rf, 7), mp.nstr(Rf*Q, 7),
            mp.nstr(mm**2*(Rp - 1), 5), mp.nstr(mm**2*(Rp*Q - 1), 5)))

if __name__ == "__main__":
    check1(); check2(); check3(); check4(); check5(); check6()
    open(os.path.join(here, 'out_mft_heavy.txt'), 'w').write("\n".join(out) + "\n")
