"""
f1: bulk-point (flat-limit) densities of the graviton functionals, from the CLS-D factorisation (arXiv:2311.04271 sec. 3.3-3.5).

CLS-D: in the bulk-point limit (nu, m >> 1, x = nu/m < 1 fixed) the action of a light-ray-commutator functional Psi^{(a),(c)}_{k,nu,j'}[f G]
on a (spinning) block localises on the "scattering crystal"; the fast factor is the same for scalars and spinning operators, and the subtraction
factor f(u',v') (nu- and m-independent) enters only through its value at the crystal (their eq. for Psi_{k,nu,f^scalar} and sec. 3.5.2).
Hence two functionals with the SAME (a),(c) (same light-ray operators, dual and celestial structures) and different weights f1, f2 have
bulk-point actions in the ratio f1/f2 |_crystal.  CLS-D's own weights f^{+-}_{k;J_e} give exactly the flat sum rules (shock amplitudes),
i.e. the B' densities.  So the bulk-point density of the project's Psi_k (weights A_k, OMEGA_NF/SUPERCONVERGENT_UPLIFT) relative to B' is
   R_k(x) = [A_k/f_{k;2}](crystal) / [A_k/f_{k;2}](x -> 0)      (normalised as the leading Regge kernel, OMEGA_PRIME N_k convention).

Checks done here (sympy, exact):
 (1) the crystal is the CHMRS Regge-frame point eta = 1, r = r*(x) = (1-s)/(1+s), s = sqrt(1-x^2): all 9 printed CLS-D saddle values
     f^{+-}_{k;J_e}|crystal (J_e = 0: k = 2, 4, 6; J_e = 2: k = 2..7) are reproduced (up to the u'<->v' orientation sign);
 (2) scalar cross-check: with the scalar CLS-D weights f^scalar_k (= CHMRS C~_k) and the scalar-analogue weights A = r^{k+1} of
     OBJECT_IDENTITY.md sec. 2.3, the method gives (1-r*^2)^{k+2}/(1+r*^2), i.e. OBJECT_IDENTITY's independently derived (and refereed)
     Psi/C density ratio, for k = 2, 4, 6;
 (3) leading Regge (r -> 0) form of f^{+-}_{k;2}: c_k eta^j r^{k-3}(1 + O(r^2)) with j = 1,0,3,2,1,0 for k = 2..7 (SUPERCONVERGENT_UPLIFT
     item 1 quotes the same j); the project's A_k are eta-free at leading order;
 (4) R_k(x) for the graviton Psi_k: k = 3, 7: 1 (A_3, A_7 are CLS-D's weights); k = 4: (1-r^2)^2/(1+r^2); k = 6: (1-r^2)^4/(1+r^2);
     k = 5: (1+6r^2+r^4)/(1-r^2)^2, i.e. in s: 8s^2/((1+s)^2(1+s^2)), 128 s^4/((1+s)^6(1+s^2)), (1+s^4)/(2s^2).
Usage: python3 f1_bulkpoint.py -> cache/f1_bulkpoint.txt
"""
import os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))

r, eta, x, s = sp.symbols('r eta x s', positive=True)
# CHMRS u-channel Regge frame (OBJECT_IDENTITY 2.3, SPINNING_MELLIN 5.2): u' = (1-2 eta r + r^2)^2/(16 r^2), v' = (1+2 eta r + r^2)^2/(16 r^2)
U = (1 - 2*eta*r + r**2)**2/(16*r**2)
V = (1 + 2*eta*r + r**2)**2/(16*r**2)
rstar = (1 - s)/(1 + s)          # s = sqrt(1 - x^2)
def at_crystal(e):
    return sp.simplify(sp.factor(e.subs(eta, 1).subs(r, rstar)))

def f_cls(k, Je, uu, vv):
    """CLS-D eqs. (subtraction_factor_evenjprime), (subtraction_factor_oddjprime); Je = 0 reproduces (subtracction_scalars)."""
    if k % 2 == 0:
        if (k - 2*Je - 2) % 4 == 0:
            return (vv - uu)/(uu*vv)**sp.Rational(k - 2*Je + 2, 4)
        return ((vv - uu)**3 - (vv - uu))/(4*(uu*vv)**sp.Rational(k - 2*Je + 4, 4))
    if (k - 2*Je - 3) % 4 == 0:
        return -1/(uu*vv)**sp.Rational(k - 2*Je + 1, 4)
    return -((uu - vv)**2 - 1)/(4*(uu*vv)**sp.Rational(k - 2*Je + 3, 4))

def f_cls_saddle(k, Je):
    """CLS-D's printed values at the scattering crystal (m = 1, nu = x)."""
    m, nu = 1, x
    if k % 2 == 0:
        return -nu**(2*(k - 2*Je + 1))*(2*m**2 - nu**2)/(m**2*(m**2 - nu**2))/(m**2*(m**2 - nu**2))**sp.Rational(k - 2*Je, 2)
    return -nu**(2*(k - 2*Je))*nu**2/(m**2*(m**2 - nu**2))/(m**2*(m**2 - nu**2))**sp.Rational(k - 2*Je - 1, 2)

pr("# f1: bulk-point densities of graviton functionals via the CLS-D crystal (arXiv:2311.04271)")
pr("## (1) crystal = CHMRS frame point (eta=1, r=r*(x)); CLS-D orientation is u'<->v' of the CHMRS frame")
ok_all = True
for Je in (0, 2):
    for k in range(2, 8):
        if Je == 0 and k % 2:
            continue                      # no scalar odd functional
        mine = at_crystal(f_cls(k, Je, V, U))          # CLS-D orientation: their u' = our v'
        mine_x = sp.simplify(mine.subs(s, sp.sqrt(1 - x**2)))
        theirs = sp.simplify(f_cls_saddle(k, Je))
        d = sp.simplify(sp.factor(mine_x - theirs))
        rel = [float(((mine_x - theirs)/theirs).subs(x, xv)) for xv in (0.13, 0.5, 0.91)]
        ok = all(abs(v) < 1e-12 for v in rel)
        ok_all &= ok
        pr(f"   J_e={Je} k={k}: f|crystal (ours) = {sp.factor(mine_x)};  CLS-D = {sp.factor(theirs)};  rel.diff at x=.13,.5,.91: "
           + ", ".join(f"{v:.1e}" for v in rel) + ("  OK" if ok else "  MISMATCH"))
pr(f"   all match: {ok_all}")

pr("## (3) leading Regge form of the CLS-D graviton weights: f_{k;2} = c_k eta^j r^(k-3) (1 + O(r^2))  [CHMRS orientation]")
lead = {}
for k in range(2, 8):
    f = sp.simplify(f_cls(k, 2, U, V))
    ser = sp.series(f, r, 0, k - 3 + 1).removeO()
    c = sp.factor(sp.simplify(ser.coeff(r, k - 3)))
    lead[k] = c
    pr(f"   k={k}: r^{k-3} coefficient = {c}")

pr("## project weights A_k (OMEGA_NF / SUPERCONVERGENT_UPLIFT), CHMRS orientation (v'-u' > 0 on the moment domain)")
A = {3: sp.Integer(1), 5: (U + V)/(32*U*V), 7: 1/(256*U*V), 4: r, 6: r**3}   # even k: sgn(v'-u') r^(k-3), sgn = +1 on the domain
for k in (4, 5, 6):
    ser = sp.series(A[k], r, 0, k - 3 + 1).removeO()
    pr(f"   A_{k}: leading Regge term {sp.factor(ser.coeff(r, k - 3))} r^{k-3}  (eta-free)")

pr("## (2) scalar cross-check: A = r^(k+1) against f^scalar_k (= CHMRS C~_k) must give (1-r*^2)^(k+2)/(1+r*^2) (OBJECT_IDENTITY 2.5)")
for k in (2, 4, 6):
    fs = f_cls(k, 0, U, V)
    ratio = sp.simplify((r**(k + 1)/fs).subs(eta, 1))
    norm = sp.limit(ratio, r, 0)
    Rk = sp.factor(sp.simplify(ratio/norm))
    target = (1 - r**2)**(k + 2)/(1 + r**2)
    pr(f"   k={k}: R = {Rk};  equals OBJECT_IDENTITY: {sp.simplify(Rk - target) == 0}")

pr("## (4) graviton: R_k = [A_k/f_{k;2}](eta=1, r) / (r -> 0), density of the project's Psi_k relative to B' (= CLS-D's flat sum rule)")
Rg = {}
for k in range(3, 8):
    f = f_cls(k, 2, U, V)
    ratio = sp.simplify((A[k]/f).subs(eta, 1))
    norm = sp.limit(ratio, r, 0)
    Rk = sp.factor(sp.simplify(ratio/norm))
    Rs = sp.factor(sp.simplify(Rk.subs(r, rstar)))
    Rg[k] = (Rk, Rs)
    pr(f"   k={k}: R_k(r) = {Rk};   in s = sqrt(1-x^2): {Rs}")
pr("   table R_k(x) at x = 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99:")
for k in range(3, 8):
    vals = [float(Rg[k][1].subs(s, sp.sqrt(1 - xv**2))) for xv in (0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99)]
    pr(f"     k={k}: " + "  ".join(f"{v:.6g}" for v in vals))
pr("   small-x: R_4 - 1, R_5 - 1, R_6 - 1 at x = 0.1: " +
   ", ".join(f"{float(Rg[k][1].subs(s, sp.sqrt(1 - sp.Rational(1, 100)))) - 1:.2e}" for k in (4, 5, 6)))
for k in (4, 5, 6):
    pr(f"   series in x: R_{k} = {sp.series(Rg[k][1].subs(s, sp.sqrt(1 - x**2)), x, 0, 9)}")

open(os.path.join(HERE, 'cache', 'f1_bulkpoint.txt'), 'w').write("\n".join(out) + "\n")
