"""u1: flat (bulk-point) limit of Omega^U in the B' dictionary.

(a) [D, exact sympy] pointwise density identity  sum_k phi_k(p) eps_k m^{2k} C_k^flat = psi2 (B^(1)_2 density) + psi3 (B^(1)_3 density)
    for both helicities (stripped +- density with dtilde44 = d44/(1-y)^4), and the forward map  a(2JJ - c_h) + 5 a sigma -> stripped.
(b) [NE] numerically: sum_k part_k + FW  (high-order quadrature) = stripped flat action of rcommon.flat_parts, for M02, M01, F_g3.
(c) [NE] cancellation amplification N_U/N_s, N_U = sum_k |part_k| + |FW|, N_s = |B2| + |B3| + |FW| (stripped component norm), vs m/Lam.
(d) [NE] odd-k prefactor sensitivity: replacing the odd B' density sigma/(1-y)^{(k+1)/2} by the scalar continuation (eps_k E_k) changes
    the flat action by dOdd; reported relative to N_s and to the action.
Output: cache/u1.log
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
from common import load, flat_unstripped_parts, flat_parts, SIGMA, KS, EPS, ek, Ek_scalar, Dh_all, peval, FUNS

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

# ---------------------------------------------------------------- (a) exact identity
p, y, u, P, d, JJ, sg = sp.symbols('p y u P d JJ sigma')
s2, s3 = sp.symbols('psi2 psi3')
q = sp.Symbol('q', positive=True)
PHI = {2: s2, 3: s3 - 4*p**2*s2, 4: 2*p**4*s2 - p**2*s3, 5: -3*p**6*s2 + 2*p**4*s3, 6: (p**8*s2 - p**6*s3)/2, 7: (-p**10*s2 + p**8*s3)/2}
okall = True
for h, sgv in (('++', 1), ('+-', -1)):
    D = P if h == '++' else d          # d = full d44; stripped +- uses dtilde44 = d/(1-y)^4
    yy = u*p**2
    lhs = sum(u**(k - 2)*PHI[k]*(((2 - yy)/(1 - yy)**sp.Rational(k + 2, 2)) if k % 2 == 0 else sgv/(1 - yy)**sp.Rational(k + 1, 2))*D for k in KS)
    Dt = P if h == '++' else d/(1 - yy)**4
    rhs = s2*(2 - yy)*Dt + sgv*u*s3*Dt
    ok = sp.simplify(sp.together(lhs - rhs)) == 0
    okall &= ok
    pr(f"(a) {h}: sum_k u^(k-2) phi_k e_k D_h  ==  psi2 (2-y) Dtilde + sigma u psi3 Dtilde  (pointwise, exact): {ok}")
# forward: -a d_y[E_4 D](0) + (5a/2) * 2 sigma D(1);  P_J'(1) = JJ/2, d44'(1) = JJ/2 - 8 (as functions of cos theta)
for h, dprime, sgv, target in (('++', JJ/2, 1, 2*JJ), ('+-', JJ/2 - 8, -1, 2*JJ - 42)):
    E4 = (2 - y)/(1 - y)**3
    val = -(sp.diff(E4, y).subs(y, 0)*1 + E4.subs(y, 0)*(-2)*dprime) + 5*sgv
    ok = sp.simplify(val - target) == 0
    okall &= ok
    pr(f"(a) forward {h}: -d_y[E_4 D](0) + 5 sigma = {sp.expand(val)}  (stripped: {target})  {ok}")
pr(f"(a) VERDICT exact identities: {'PASS' if okall else 'FAIL'}")

# ---------------------------------------------------------------- (b)-(d) numerics
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(1500); pn = (xg + 1)/2; pw = wg/2
worst = 0.0
for name in FUNS:
    F, PH = load(name)
    pr(f"== {name}: a = {F.a:.6g}, k(0) = {F.k():.6g}, phi_k degrees = {[len(PH[k]) - 1 for k in KS]}, phi_k(0) = {[round(float(PH[k][0]), 4) for k in KS]}")
    for h in SIGMA:
        for mfl in (1.0005, 1.01, 1.05, 1.2, 1.5, 2.0, 3.0, 6.0, 20.0):
            u = 1/mfl**2
            Jmax = 60
            Js, parts, FW = flat_unstripped_parts(F, PH, h, u, Jmax, pn, pw)
            Jsr, B2, B3, FWs = flat_parts(F, h, u, Jmax)
            tot = parts.sum(axis=1) + FW
            ref = B2 + B3 + FWs
            Ns = np.abs(B2) + np.abs(B3) + np.abs(FWs)
            NU = np.abs(parts).sum(axis=1) + np.abs(FW)
            sel = (Js % 2 == 0) if h == '++' else np.ones_like(Js, bool)
            err = np.max(np.abs(tot - ref)[sel]/Ns[sel])
            worst = max(worst, err) if mfl >= 1.01 else worst
            # odd-prefactor sensitivity
            yv = u*pn*pn
            _, Dm = Dh_all(h, Jmax, yv)
            dodd = sum(u**(k - 2)*(Dm @ (pw*peval(PH[k], pn)*SIGMA[h]*(EPS[k]*Ek_scalar(k, yv) - 1/(1 - yv)**((k + 1)/2)))) for k in KS if k % 2)
            amp = NU[sel]/Ns[sel]
            pr(f"  {h} m/Lam={mfl:8.4f}: max|sum - stripped|/N_s = {err:.2e};  N_U/N_s: max {amp.max():9.3g} (J={Js[sel][np.argmax(amp)]}), J<=6 {amp[:4].max():9.3g};"
               f"  odd-prefactor shift max|dOdd|/N_s = {np.max(np.abs(dodd[sel])/Ns[sel]):.3g}, min (F+dOdd)/N_s = {np.min((ref + dodd)[sel]/Ns[sel]):+.4f} vs min F/N_s = {np.min(ref[sel]/Ns[sel]):+.4f}")
pr(f"(b) VERDICT flat-limit identity (numeric, m >= 1.01 Lam, J <= 60): max relative deviation {worst:.2e} -> {'PASS' if worst < 1e-8 else 'FAIL'}")
open(os.path.join(HERE, 'cache', 'u1.log'), 'w').write("\n".join(out) + "\n")
