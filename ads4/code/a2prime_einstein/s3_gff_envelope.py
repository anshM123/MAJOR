"""S3 [computed exponents; overall constants not tracked]. Large-J size of the GFF [TT] OPE coefficients, hence of the
stripped detector amplitudes, on every light family -- and two defects found on the way.

(a) 1d GFF coefficients.  In the mirror ordering of Sec. sec:angcross the crossed identity (1-z)^{-2b} decomposes as
    sum_n c_n z^n 2F1(H+d, H+d; 2H; z), H = a+b+n, d = b-a, with c_n = (2a)_n (2b)_n / (n! (2a+2b+n-1)_n) > 0
    (checked against a forward solve), and c_n 4^H ~ H^{2a+2b-3/2}.
    A collinear component T_m is an sl(2)xsl(2) primary of weights (h_m, hb_m) = ((3-m)/2, (3+m)/2).  So for the pair
    (T_{m1}, T_{m3}) the GFF weight per state at fixed h and large s = (Delta+J)/2 is
        4^Delta 2^{-J} |F_{m1 m3}|^2  ~  s^{9/2 + m1 + m3}          (summed over a shell [S,2S]: S^{11/2+M}, M = m1+m3,
    exactly the power of the two-scale crossing bound, which GFF saturates).  At h = 1 only (m1,m3) = (2,2) is populated
    (h_m = 1/2 only for m = 2): the 2d state with h = 1 is the e_+ component of the twist-2 primaries alone.
(b) Collinear identities for ALL components (angproofs/ap3_collinear.py, symbolic structures, read-only): the m <= -1
    components used in the Supplement vanish identically at h = 1,2,3,4 (even) and h = 3/2,5/2,7/2 (odd), i.e. exactly on the
    light double-twist families below h = 5 -- so the "h < 5, same inequalities, D-sketch" extension of Lemma 4 cannot work
    with them.  The m >= +1 components do not vanish there: F_{2,2} (no h-zero), F_{1,1} (zero only at h = 1), F_{1,2}
    (odd, no h-zero).  Inverting them gives per state
        4^Delta 2^{-J} |lam_pp|^2, |lam_pm|^2 ~ s^{1/2},   4^Delta 2^{-J} |lam_odd|^2 ~ s^{3/2}     (h >= 2 even; all odd),
    the same powers as the Supplement's h >= 5 result.  With S2 (stripped kernels ~ J^{-13/2}, J^{-15/2}):
        stripped |a_{[TT]_{n,J},b,0}|^2 / C_T  ~  J^{1/2-13/2} = J^{3/2-15/2} = J^{-6}   per state, fixed n.
(c) h = 1 (n = 0, parity even): the two canonical even structures become proportional (rank 1), and
    N(1) := lim_{h->1} (S^{+-} - c S^{++})/(h-1) is a second conserved twist-2 structure whose detector amplitude on the b = 4
    ladder is lim M4_pm(h)/(h-1) != 0 (no zero).  In GFF its coefficient B vanishes (<T|E E|T>_GFF = 0), but its O(1/c_T)
    value B^(1)_J is an independent tree-level datum, NOT controlled by gamma: |a_{[TT]_{0,J},4,0}|^2 ~ c_T |B^(1)_J/c_T|^2
    x (envelope).  So the OPE-data form (ii') of A2' misses one datum on the n = 0 even family."""
import os, pickle
import mpmath as mp
import sympy as sp
from common_a2 import Tee
import ap3_collinear as ap
from common import S, J, CACHE
import closed_forms as cf
mp.mp.dps = 50
log = Tee('s3_gff_envelope.txt')

# ---------------- (a)
def c_closed(a, b, n):
    return mp.rf(2*a, n)*mp.rf(2*b, n)/(mp.factorial(n)*mp.rf(2*a + 2*b + n - 1, n))
def c_forward(a, b, N):
    d = b - a; c = []
    for k in range(N):
        s = mp.rf(2*b, k)/mp.factorial(k)
        for n in range(k):
            H = a + b + n; m = k - n
            s -= c[n]*mp.rf(H + d, m)**2/(mp.rf(2*H, m)*mp.factorial(m))
        c.append(s)
    return c
log("(a) 1d GFF coefficients")
for a, b in ((0.5, 0.5), (0.5, 2.5), (2.5, 2.5), (1, 2)):
    a, b = mp.mpf(a), mp.mpf(b)
    cf_ = c_forward(a, b, 25)
    dev = max(abs(cf_[n]/c_closed(a, b, n) - 1) for n in range(25))
    f = lambda n: c_closed(a, b, n)*mp.mpf(4)**(a + b + n)
    ex = [mp.log(f(2*n)/f(n))/mp.log((a + b + 2*n)/(a + b + n)) for n in (10**3, 10**4, 10**5)]
    log(f"   (a,b)=({float(a)},{float(b)}): closed form vs forward solve {mp.nstr(dev,3)};  exponent of c_n 4^H:"
        f" {', '.join(f'{float(x):.5f}' for x in ex)}  (2a+2b-3/2 = {float(2*a+2*b-1.5)})")

# ---------------- (b)
structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
hh = sp.symbols('h')
F = {}
for name in ('pp', 'pm', 'odd', 'oddJ'):
    Fn, _ = ap.collinear(structs[name], 'even' if name in ('pp', 'pm') else 'odd')
    F[name] = {k: sp.factor(v.subs(S, hh + J)) for k, v in Fn.items()}
log("\n(b) zeros in h of the collinear components (per canonical structure), light range h < 5:")
for name in ('pp', 'pm', 'odd', 'oddJ'):
    for k in sorted(F[name]):
        roots = sorted(set(r for r in sp.solve(sp.numer(F[name][k]), hh) if r.is_number and 0 < r < 5))
        log(f"   {name:4s} F{k}: deg_s {sp.degree(sp.expand(F[name][k].subs(J, S - hh)), S)}, h-zeros in (0,5): {roots}")
# per-state lambda^2 exponents from (2,2), (1,1), (1,2)
log("\n   per-state exponent of 4^Delta 2^-J |lambda|^2 = (9/2 + m1 + m3) - 2 deg_s F:")
for name, k in (('pp', (2, 2)), ('pm', (2, 2)), ('pp', (1, 1)), ('pm', (1, 1)), ('odd', (1, 2)), ('oddJ', (1, 2)),
                ('pp', (-2, -2)), ('odd', (-2, -1))):
    deg = sp.degree(sp.expand(F[name][k].subs(J, S - hh)), S)
    log(f"   {name:4s} via F{k}: {sp.Rational(9, 2) + k[0] + k[1] - 2*deg}")
log("   F(2,2) = s^4 (lam_pm - 2 lam_pp)/16 + O(s^3);  F(1,1)/(1-h) = s^3 (2 lam_pp + lam_pm) + O(s^2): both combinations,")
log("   hence lam_pp and lam_pm separately, are bounded by s^{1/4} 2^{J/2} 4^{-Delta/2} for h >= 2.")

# ---------------- (c)
keys = sorted(set(structs['pp']) | set(structs['pm']), key=str)
def vec(name, hv, Jv):
    return [sp.sympify(structs[name].get(k, 0)).subs({S: hv + Jv, J: Jv}) for k in keys]
log("\n(c) rank of {S^{++}, S^{+-}} (13 monomial coefficients) on the light even families:")
for hv in (1, 2, 3, 4, sp.Rational(11, 2)):
    log(f"   h = {hv}: rank " + ", ".join(f"J={Jv}: {sp.Matrix([vec('pp', hv, Jv), vec('pm', hv, Jv)]).rank()}"
                                       for Jv in (6, 12, 40)))
k0 = [k for k in keys if sp.sympify(structs['pp'].get(k, 0)).subs({S: J + 1}) != 0][0]
c = sp.simplify(sp.sympify(structs['pm'][k0]).subs(S, J + 1)/sp.sympify(structs['pp'][k0]).subs(S, J + 1))
Nvec = {k: sp.factor(sp.limit((sp.sympify(structs['pm'].get(k, 0)) - c*sp.sympify(structs['pp'].get(k, 0))).subs(S, hh + J)
                               /(hh - 1), hh, 1)) for k in keys}
log(f"   at h = 1: S^{{+-}} = c S^{{++}}, c = {c};  N(1) has {sum(1 for v in Nvec.values() if v != 0)} nonzero monomials"
    f" (so it is a second, independent conserved twist-2 structure)")
log("   detector amplitude of N on the b=4 ladder, lim_{h->1} M4_pm(Delta,J)/(h-1) (normalisation of closed_forms.py):")
for Jv in (6, 10, 40, 200):
    vals = [abs(cf.M4_pm(Jv + 2*(1 + e), Jv))/e for e in (mp.mpf('1e-10'), mp.mpf('1e-20'))]
    log(f"     J={Jv}: |a_4(N)| = {mp.nstr(vals[1], 8)}  (eps 1e-10 -> 1e-20 rel. change {mp.nstr(abs(vals[0]/vals[1]-1), 2)})")
log("   => nonzero: the N-direction is not at a double zero.  Its O(1/c_T) coefficient must be added to (ii').")
