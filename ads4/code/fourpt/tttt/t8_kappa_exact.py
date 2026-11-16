"""
Exact shadow eigenvalue of CHL's helicity structures, kappa_h:  T T^{E,h}(Delta,J) = kappa_h T^{E,h}(3-Delta,-1-J)  (T = word-wise shadow map
c_w -> c_w b^(0)/b_w of fourpt), and the identity  kappa_h A_h = 2^{3J+3/2}  with A_h the MFT helicity factors (now derived from our own Fock
data, t7_fock_l2).  Method (no conservation solve needed): CHL's structure in position space (chl_xspace) is T^{E,h} = C_l(Delta,J) X_h with
X_h rational and
    C_l(Delta,J) = (-1)^l 2^{3l-1/2} 2^{J/2} ((2-Delta+J)/2)_{2l} ((Delta+J+1-4l)/2)_{2l}      [all pi's cancel; derived in SPINNING_TTTT.md 3.5];
X_h(Delta) and X_h(Delta~) are decomposed exactly into the differential-basis words (modulo the d = 3 relation for l = 2), T is applied, and
lambda_h is read off from  T X_h(Delta) = lambda_h X_h(Delta~)  (exact; the residual must vanish identically);  kappa_h = lambda_h C(Delta)/C(Delta~).
Frame: CHL (O_1(0), O_2(x)); our word structures (-zeta.x/2)^J |x|^{-2al} R_w map to it by x -> -x.
Usage: python3 t8_kappa_exact.py l Dv Jv [Dv Jv ...]       Output appended to out_t8_kappa_exact.txt
"""
import os, sys, time
import sympy as sp
from chl_xspace import xspace_structure
from ft_helicity import R_inv, AX, BX, WX, XX, AB, AW, BW
from infframe import words
from embed2 import shifts, Delta as sD, J as sJ
from t1_structures_l1 import b_ratio
INV = (AX, BX, WX, XX, AB, AW, BW)

def C_l(l, D, J):
    return (-1)**l*2**(3*l - sp.Rational(1, 2))*2**(J/2)*sp.rf((2 - D + J)/2, 2*l)*sp.rf((D + J + 1 - 4*l)/2, 2*l)

def decompose(l, Dv, Jv, kind):
    """word coefficients of X_h at (Dv, Jv): particular solution and null directions"""
    al = (4*l + 2 - Dv + Jv)/sp.Integer(2)
    Xn = sp.expand(sp.powsimp(sp.expand(xspace_structure(l, Dv, Jv, kind)*WX**(-Jv)*XX**al)))
    ws = words(l)
    Rw = {w: sp.together(R_inv(l, {w: 1}, Dv, Jv).xreplace({AX: -AX, BX: -BX, WX: -WX})*sp.Integer(2)**(-Jv)) for w in ws}
    cs = sp.symbols('c0:%d' % len(ws))
    rel = 0
    if l == 2:     # d = 3: the Gram determinant of (x, a, b, zeta) vanishes; it is the source of the word relation
        G4 = sp.Matrix([[XX, AX, BX, WX], [AX, 0, AB, AW], [BX, AB, 0, BW], [WX, AW, BW, 0]]).det()
        mus = sp.symbols('mu0:7'); cs = cs + mus; rel = sum(mus[k]*G4*XX**(k - 3)/WX**2 for k in range(7))   # homogeneity picks XX^2
    expr = sp.together(Xn - sum(cs[i]*Rw[w] for i, w in enumerate(ws)) - rel)
    num = sp.expand(sp.fraction(expr)[0])
    eqs = sp.Poly(num, *INV).coeffs()
    A, b = sp.linear_eq_to_matrix(eqs, cs)
    sol, params = A.gauss_jordan_solve(b)
    return ws, sol[:len(ws), :], params

def kappa(l, Dv, Jv, kind):
    Dt, Jt = 3 - Dv, -1 - Jv
    ws, s1, p1 = decompose(l, Dv, Jv, kind)
    s1 = s1.subs({p: 0 for p in p1})                     # any representative (relation ambiguity)
    ws2, s2, p2 = decompose(l, Dt, Jt, kind)
    Tc = sp.Matrix([s1[i]/b_ratio(shifts(w)[0] - shifts(w)[1]).subs({sD: Dv, sJ: Jv}) for i, w in enumerate(ws)])
    lam = sp.Symbol('lam')
    nu = sp.symbols('nu0:%d' % max(1, len(p2)))
    # T c = lam * (s2 particular + sum tau_i r_i)  ->  unknowns lam, nu_i = lam tau_i
    base = s2.subs({p: 0 for p in p2})
    dirs = [s2.diff(p) for p in p2]
    eqs = list(Tc - lam*base - sum((nu[i]*dirs[i] for i in range(len(p2))), sp.zeros(len(ws), 1)))
    sol = sp.solve(eqs, [lam] + list(nu[:len(p2)]), dict=True)
    assert len(sol) == 1, "T X_h(Delta) is not proportional to X_h(Delta~)"
    lam_v = sp.nsimplify(sp.simplify(sol[0][lam]))
    kap = sp.simplify(lam_v*C_l(l, Dv, Jv)/C_l(l, Dt, Jt))
    A = sp.rf(-Jv, 2*l)/sp.rf(Jv + 1, 2*l) if kind == 'opp' else sp.rf(Dv - 1, 2*l)/sp.rf(2 - Dv, 2*l)
    test = sp.nsimplify(sp.simplify(sp.powsimp(kap*A/2**(3*Jv + sp.Rational(3, 2)), force=True)))
    return lam_v, kap, test, len(p1), len(p2)

if __name__ == "__main__":
    l = int(sys.argv[1]); pts = sys.argv[2:]
    out = []
    for i in range(0, len(pts), 2):
        Dv, Jv = sp.Rational(pts[i]), sp.Rational(pts[i + 1])
        for kind in ('same', 'opp'):
            t = time.time()
            lam, kap, test, n1, n2 = kappa(l, Dv, Jv, kind)
            s = (f"  l={l} (Delta,J)=({Dv},{Jv}) {kind}: word-relation dims {n1},{n2}; lambda_h = {lam};  kappa_h = {sp.N(kap, 15)};  "
                 f"kappa_h A_h / 2^(3J+3/2) = {test}   ({time.time() - t:.0f} s)")
            print(s, flush=True); out.append(s)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t8_kappa_exact.txt'), 'a') as f:
        f.write("\n".join(out) + "\n")
