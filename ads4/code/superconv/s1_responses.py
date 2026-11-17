"""s1: leading heavy responses of the superconvergent replacements C^SC_k, k = 3..7 (SUPERCONVERGENT_UPLIFT.md).

Replacement for moment k: Psi_k = int K_k <T [L[T](x0,z3), L[T](x0,z1)] T>, K_k = A_k(r) x [covariant celestial structure of parity j' at (1,3)] x [(2,4) lift],
A_k = r^{k-3} (k odd, symmetric) or sgn(v'-u') r^{k-3} (k even, antisymmetric).  Leading Regge kernel: r^{k-2} x antipodal weights w^{(j')}_q(nu) x F^nu_M x (2,4) kernel,
j' = k mod 2.  Responses A(h1,h3) = sum_q w_q M^q_{h1h3}(nu; a), a = 2l + k + 1 = k + 5 (l = 2), from e2e/laplace.py M-tables (imported, cached).
Project kernels at one pair: even D1^2 (C^+), odd A = D5 D1 (C^-).
Outputs per k: N_k(nu) (ratio of preserved responses, sigma_h structure), flip/pres ratios at (1,3) and of the project kernel; zeros on real nu.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e')); os.chdir(os.path.join(HERE, '..', 'e2e'))
import sympy as sp, mpmath as mp
from laplace import M_table, nu, a, g1

delta = sp.Rational(1, 2) + sp.I*nu
L = 2

def w_celestial(jp):
    w = {}
    for q1 in range(-L, L+1):
        for q3 in range(-L, L+1):
            M = q1 + q3
            if M % 2 != jp: continue
            if jp == 0:
                n = abs(M)//2
                c = sp.rf(delta/2, n)/sp.rf(1 - delta/2, n)
            else:
                mm = (M - 1)//2; s = 1
                if mm < 0: mm = -mm - 1; s = -1
                c = s*sp.rf((delta + 1)/2, mm)/sp.rf((3 - delta)/2, mm)
            w[(q1, q3)] = sp.Integer(-1)**q3*c
    return w

# project kernels in the frame z(th) = yhat + cos th e_par + sin th e_perp:  D1 = cos(th1-th3), D5 = sin(th3-th1)
W_D1sq = {(2, -2): 4, (-2, 2): 4, (0, 0): sp.Rational(2, 9)}   # contraction weights (SPINNING_E2E 3.1; checked in s0)
W_A = {(2, -2): sp.I/4, (-2, 2): -sp.I/4}            # sin(th3-th1)cos(th1-th3) = sin(2(th3-th1))/2

T = M_table(L)
def resp(w, h1, h3, aa):
    tot = 0
    for q, ww in w.items():
        R, S, mu = T[(q[0], q[1], h1, h3)]
        tot += ww*(R + g1*S)
    return sp.simplify(sp.together(sp.expand(tot.subs(a, aa))))

def gam(aa, x): return abs(mp.gamma((aa - mp.mpf(1)/2 + 1j*x)/2))**2

if __name__ == "__main__":
    mp.mp.dps = 20
    out = {}
    for k in range(3, 8):
        aa = k + 5; jp = k % 2
        wc = w_celestial(jp); wp = W_A if jp == 1 else W_D1sq
        R = {}
        for name, w in (('SC', wc), ('proj', wp)):
            for hh in ((1, -1), (-1, 1), (1, 1), (-1, -1)):
                R[(name,) + hh] = resp(w, hh[0], hh[1], aa)
        print(f"=== k = {k} (a = {aa}, j' = {jp}, project kernel {'A=D5D1' if jp else 'D1^2'}) ===")
        for name in ('SC', 'proj'):
            pe = sp.simplify((R[(name, 1, -1)] + R[(name, -1, 1)])/2); po = sp.simplify((R[(name, 1, -1)] - R[(name, -1, 1)])/2)
            fe = sp.simplify((R[(name, 1, 1)] + R[(name, -1, -1)])/2); fo = sp.simplify((R[(name, 1, 1)] - R[(name, -1, -1)])/2)
            print(f"  {name:4s}: pres even = {sp.factor(pe)}")
            print(f"        pres odd  = {sp.factor(po)}")
            print(f"        flip even = {sp.factor(fe)}")
            print(f"        flip odd  = {sp.factor(fo)}")
            out[(k, name)] = (pe, po, fe, fo)
        idx = 1 if jp else 0
        Nk = sp.factor(sp.simplify(out[(k, 'SC')][idx]/out[(k, 'proj')][idx]))
        print(f"  N_{k}(nu) = {Nk}")
        num, den = sp.fraction(sp.together(Nk))
        if not Nk.has(g1):
            rn = sp.Poly(sp.expand(num), nu).nroots() if sp.Poly(sp.expand(num), nu).degree() > 0 else []
            rd = sp.Poly(sp.expand(den), nu).nroots() if sp.Poly(sp.expand(den), nu).degree() > 0 else []
            print(f"  zeros of N_{k} in nu: {[complex(r) for r in rn]}")
            print(f"  poles of N_{k} in nu: {[complex(r) for r in rd]}")
        f = sp.lambdify((nu, g1), Nk, 'mpmath')
        mins = min(abs(complex(f(mp.mpf(x), gam(aa-1, x)/gam(aa, x)))) for x in [i*0.05 for i in range(0, 401)])
        print(f"  min |N_{k}| on nu in [0,20] (step 0.05): {mins:.4e}")
        if jp == 0:
            # flip interference: ratio flip-even / pres-even at (1,3), for SC and for the project kernel D1^2
            for name in ('SC', 'proj'):
                rr = sp.simplify(out[(k, name)][2]/out[(k, name)][0])
                fr = sp.lambdify((nu, g1), rr, 'mpmath')
                vals = [(x, complex(fr(mp.mpf(x), gam(aa-1, x)/gam(aa, x)))) for x in (0, 0.5, 1, 2, 3, 5, 10, 20)]
                print(f"  {name}: flip/pres at (1,3) = {sp.factor(rr)}")
                print("        " + "  ".join(f"nu={x}: {abs(v):.3e}" for x, v in vals))
