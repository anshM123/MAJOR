"""
Item 4: remainder bounds under the hypothesis |gamma_O| <= C / Lambda^(2k-2) = C x^(k-1), x = 1/Lambda^2,
k = (number of derivatives)/2, for monomials in the basis used by the collaborator.               [computed]

(a) Algebraic tower I^m J^n (shear) and I^m J^n P, I^m J^n Q (Hall): exact |response| from the verified formulas
      dV_E = 432 F_II - 648 F_IJ - 1944 F_JJ - 108 F_J,
      dV_H = 108 F_Q + 1296 F_QI + 1944 F_QJ - 432 F_PI - 648 F_PJ,
    summed with weight x^k by brute force, compared with the closed forms E(x), H(x); remainder = C * sum_{k>=5}/x.
(b) Gradient tower F(y, X), y = (I, J, P, Q), X_ij = grad y_i . grad y_j: exact responses of every monomial with
    X-degree <= 2 (X-degree >= 3 gives exactly zero: X_ij = f y_i' y_j' vanishes at the horizon to that order),
    from A_i = F_{y_i} - (1/r^2) d_r( r^2 f sum_j c_ij F_{X_ij} y_j' ) and
      dV_E = -108 A_J - 6 A_I' + 18 A_J',   dV_H = 108 A_Q + 6 A_P' - 18 A_Q'.
"""
import itertools
from fractions import Fraction as Fr
import sympy as sp

x = Fr(1, 997002)            # 1/Lambda^2 at Delta_g = 1000
A, B = 12 * x ** 2, 12 * x ** 3
E = (864 * x ** 4 / ((1 - A) ** 3 * (1 - B)) + 648 * x ** 5 / ((1 - A) ** 2 * (1 - B) ** 2)
     + 3888 * x ** 6 / ((1 - A) * (1 - B) ** 3) + 108 * x ** 3 / ((1 - A) * (1 - B) ** 2))
H = (432 * x ** 4 / ((1 - A) ** 2 * (1 - B)) + 648 * x ** 5 / ((1 - A) * (1 - B) ** 2) + 108 * x ** 3 / ((1 - A) * (1 - B))
     + 1296 * x ** 5 / ((1 - A) ** 2 * (1 - B)) + 1944 * x ** 6 / ((1 - A) * (1 - B) ** 2))
qE = (E - 108 * x ** 3 - 864 * x ** 4) / x ** 5
qH = (H - 108 * x ** 3 - 432 * x ** 4) / x ** 5
print(f"(a) (E - 108x^3 - 864x^4)/x^5 = {float(qE):.9f}   (claim < 1945)")
print(f"    (H - 108x^3 - 432x^4)/x^5 = {float(qH):.9f}   (claim < 3241)")

# brute force check of the generating functions (weights x^k, k = 2m+3n (+3 for Q, +2 for P))
def respE(m, n):
    F_J = n * Fr(12) ** (m + n - 1) if n >= 1 else 0
    F_II = m * (m - 1) * Fr(12) ** (m + n - 2) if m >= 2 else 0
    F_IJ = m * n * Fr(12) ** (m + n - 2) if m >= 1 and n >= 1 else 0
    F_JJ = n * (n - 1) * Fr(12) ** (m + n - 2) if n >= 2 else 0
    return abs(432 * F_II) + abs(648 * F_IJ) + abs(1944 * F_JJ) + abs(108 * F_J), 432 * F_II - 648 * F_IJ - 1944 * F_JJ - 108 * F_J


def respH(m, n, odd):
    if odd == 'Q':
        F_Q = Fr(12) ** (m + n); F_QI = m * Fr(12) ** (m + n - 1) if m else 0; F_QJ = n * Fr(12) ** (m + n - 1) if n else 0
        terms = [108 * F_Q, 1296 * F_QI, 1944 * F_QJ]
    else:
        F_PI = m * Fr(12) ** (m + n - 1) if m else 0; F_PJ = n * Fr(12) ** (m + n - 1) if n else 0
        terms = [-432 * F_PI, -648 * F_PJ]
    return sum(abs(t) for t in terms), sum(terms)


SE, SH, SE5, SH5 = Fr(0), Fr(0), Fr(0), Fr(0)
k5E, k5H = [], []
for m in range(0, 40):
    for n in range(0, 30):
        k = 2 * m + 3 * n
        a, s = respE(m, n)
        SE += a * x ** k
        if k >= 5:
            SE5 += a * x ** k
            if k == 5:
                k5E.append(((m, n), s))
        for odd, kk in (('Q', k + 3), ('P', k + 2)):
            a, s = respH(m, n, odd)
            SH += a * x ** kk
            if kk >= 5:
                SH5 += a * x ** kk
                if kk == 5:
                    k5H.append(((m, n, odd), s))
print(f"    brute-force sums vs closed forms: E rel.diff = {float((SE - E) / E):.2e}, H rel.diff = {float((SH - H) / H):.2e}")
print(f"    k = 5 algebraic monomials: shear {k5E}, Hall {k5H}")
print(f"    algebraic remainder R/C (units 1/Lambda^8): shear {float(SE5 / x / x ** 4):.4f}, Hall {float(SH5 / x / x ** 4):.4f}")

# ---------------- (b) gradient tower, exact responses
r = sp.Symbol('r', positive=True)
f0 = r ** 2 - 1 / r
ysym = sp.symbols('I J P Q')
kk = dict(zip(ysym, (2, 3, 2, 3)))
bg = {ysym[0]: 12 / r ** 6, ysym[1]: 12 / r ** 9, ysym[2]: sp.Integer(0), ysym[3]: sp.Integer(0)}
pairs = [(i, j) for i in range(4) for j in range(i, 4)]
Xs = {pq: sp.Symbol(f"X_{ysym[pq[0]]}{ysym[pq[1]]}") for pq in pairs}
dy = {s: sp.diff(bg[s], r) for s in ysym}
subsbg = {**bg, **{Xs[(i, j)]: f0 * dy[ysym[i]] * dy[ysym[j]] for (i, j) in pairs}}


def A_i(F, i):
    s = ysym[i]
    term = sp.diff(F, s)
    flux = 0
    for (a, b), X in Xs.items():
        if a == i and b == i:
            flux += 2 * sp.diff(F, X) * dy[s]
        elif a == i:
            flux += sp.diff(F, X) * dy[ysym[b]]
        elif b == i:
            flux += sp.diff(F, X) * dy[ysym[a]]
    flux = flux.subs(subsbg)
    return term.subs(subsbg) - sp.diff(r ** 2 * f0 * flux, r) / r ** 2


def responses(F):
    AI, AJ, AP, AQ = [A_i(F, i) for i in range(4)]
    dVE = -108 * AJ.subs(r, 1) - 6 * sp.diff(AI, r).subs(r, 1) + 18 * sp.diff(AJ, r).subs(r, 1)
    dVH = 108 * AQ.subs(r, 1) + 6 * sp.diff(AP, r).subs(r, 1) - 18 * sp.diff(AQ, r).subs(r, 1)
    return sp.nsimplify(dVE), sp.nsimplify(dVH)


print("\n(b) checks: (grad I)^2 ->", responses(Xs[(0, 0)]), " grad I . grad P ->", responses(Xs[(0, 2)]),
      " (grad J)^2 ->", responses(Xs[(1, 1)]))
KMAX = int(__import__('sys').argv[1]) if len(__import__('sys').argv) > 1 else 10
sumE = {}; sumH = {}
cnt = 0
for nX in (1, 2):
    for Xcomb in itertools.combinations_with_replacement(pairs, nX):
        kX = sum(kk[ysym[i]] + kk[ysym[j]] + 1 for (i, j) in Xcomb)
        if kX > KMAX:
            continue
        for a in itertools.product(range(0, 6), range(0, 4), range(0, 3), range(0, 3)):
            k = kX + sum(ai * kk[s] for ai, s in zip(a, ysym))
            if k > KMAX:
                continue
            F = sp.Mul(*[Xs[c] for c in Xcomb]) * sp.Mul(*[s ** ai for s, ai in zip(ysym, a)])
            e, h = responses(F)
            cnt += 1
            sumE[k] = sumE.get(k, 0) + abs(e)
            sumH[k] = sumH.get(k, 0) + abs(h)
print(f"    {cnt} gradient monomials with k <= {KMAX} (X-degree 1, 2)")
for k in sorted(sumE):
    print(f"    k = {k:2d}: sum|dV_E| = {sumE[k]},  sum|dV_H| = {sumH[k]}")
GE6 = sum(float(v) * float(x) ** k for k, v in sumE.items() if k >= 6)
GH6 = sum(float(v) * float(x) ** k for k, v in sumH.items() if k >= 6)
print(f"    sum_(k>=6, k<={KMAX}) |resp| x^k / x^6:  shear {GE6 / float(x) ** 6:.1f}, Hall {GH6 / float(x) ** 6:.1f}   (claim < 129000)")
print(f"    k = 5 gradient monomials: shear {sumE.get(5)}, Hall {sumH.get(5)}")
