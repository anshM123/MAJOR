"""Certificate M (Mehler region): m >= 2 (u <= 1/4) and bhat = (2J+1)/m >= BA, both channels, all J.

Analytic input (Lemma M of FLAT_GRAVITON_BOUND.md): for an integer L >= 0, N = L+1/2, m > 1 and a weight w in
C^4[0,1] with w(1)=w'(1)=w''(1)=0,
   H_w(L) := int_0^1 w(p) P_L(1-2p^2/m^2) dp = w(0)/b - w''(0)/(2 b^3) + R,   b = 2N/m,   |R| <= C_R(w,u)/b^4,
   C_R(w,u) <= (2/pi)[ (2/3)(|w'''(1)| + ||w''''||_1) + 6u ||p w'''||_1 + 4u ||w''||_1 ] + u^{3/2}[ (3pi/4) sup|w''| + sup|w'| ].
 ++ : FF = H_{w+}(J) + (10JJ-5)u^3,   w+ = Psi_+ = 2psi2 + u(psi3 - p^2 psi2)
 +- : FF = sum_{j=-4}^4 c_j(J) H_{w-}(J+j) + (10JJ-210)u^3,  w- = Psi_- (1-u p^2)^{-8},  c_j>=0, sum c_j = 1.
This script bounds the norms of w_h(.,u) uniformly for u in [0,u1] and in [u1,1/4] by interval arithmetic
(Taylor jets in p, arb balls; p-cells of width 1/NP) and checks the final inequalities at b = BA (they are
monotone in b, see the .md).  Usage: python3 cert_M.py BA u1
"""
import sys
import flint
from fractions import Fraction as Fr
from common import PSI2, PSI3, set_prec

set_prec(128)
A = flint.arb
ORD = 5  # jets up to 4th derivative

def fq(x):
    x = Fr(x); return A(flint.fmpq(x.numerator, x.denominator))

class Jet:
    __slots__ = ['c']
    def __init__(self, c): self.c = list(c) + [A(0)] * (ORD - len(c))
    def __add__(self, o):
        o = o if isinstance(o, Jet) else Jet([o]); return Jet([a + b for a, b in zip(self.c, o.c)])
    __radd__ = __add__
    def __sub__(self, o):
        o = o if isinstance(o, Jet) else Jet([o]); return Jet([a - b for a, b in zip(self.c, o.c)])
    def __rsub__(self, o): return Jet([o]) - self
    def __mul__(self, o):
        if not isinstance(o, Jet): return Jet([a * o for a in self.c])
        r = [A(0)] * ORD
        for i in range(ORD):
            for j in range(ORD - i):
                r[i + j] += self.c[i] * o.c[j]
        return Jet(r)
    __rmul__ = __mul__
    def inv(self):
        a0 = self.c[0]; r = [A(0)] * ORD; r[0] = 1 / a0
        for n in range(1, ORD):
            s = A(0)
            for k in range(1, n + 1): s += self.c[k] * r[n - k]
            r[n] = -s / a0
        return Jet(r)
    def der(self, k):
        f = 1
        for i in range(2, k + 1): f *= i
        return self.c[k] * f

def poly_jet(coeffs, pj):
    r = Jet([A(0)])
    for c in reversed(coeffs): r = r * pj + Jet([fq(c)])
    return r

def weight_jet(h, p0, u):
    pj = Jet([p0, A(1)])
    psi2 = poly_jet(PSI2, pj); psi3 = poly_jet(PSI3, pj)
    sig = 1 if h == '++' else -1
    w = psi2 * (2 - u * pj * pj) + psi3 * (sig * u)
    if h == '+-':
        q = (1 - u * pj * pj).inv()
        q2 = q * q; q4 = q2 * q2; w = w * (q4 * q4)
    return w

def norms(h, ulo, uhi, NP):
    u = fq(ulo).union(fq(uhi))
    n1 = {'w4': A(0), 'pw3': A(0), 'w2': A(0)}; s2 = A(0); s1 = A(0)
    for i in range(NP):
        p0 = fq(Fr(i, NP)).union(fq(Fr(i + 1, NP)))
        w = weight_jet(h, p0, u)
        wd = [w.der(k) for k in range(5)]
        width = fq(Fr(1, NP))
        n1['w4'] += A(wd[4].abs_upper()) * width
        n1['pw3'] += A((p0 * wd[3]).abs_upper()) * width
        n1['w2'] += A(wd[2].abs_upper()) * width
        s2 = s2.max(A(wd[2].abs_upper())); s1 = s1.max(A(wd[1].abs_upper()))
    w1 = weight_jet(h, A(1), u); w3at1 = A(w1.der(3).abs_upper())
    uu = fq(uhi)
    CR = 2 / A.pi() * (A(2) / 3 * (w3at1 + n1['w4']) + 6 * uu * n1['pw3'] + 4 * uu * n1['w2']) \
         + uu * uu.sqrt() * (3 * A.pi() / 4 * s2 + s1)
    return A(CR.upper()), {k: A(v.upper()) for k, v in n1.items()}, w3at1

def main():
    BA = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(120)
    u1 = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(1, 100)
    NP = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    ok_all = True
    b = fq(BA)
    for h in ['++', '+-']:
        CR1, n1, w3 = norms(h, 0, u1, NP)
        CRm, _, _ = norms(h, u1, Fr(1, 4), NP)
        print('M %s: u in [0,%s]: |w\'\'\'(1)|<=%s ||w\'\'\'\'||_1<=%s  C_R<=%s ;  u in [%s,1/4]: C_R<=%s' % (
            h, u1, w3.str(6), n1['w4'].str(6), CR1.str(6), u1, CRm.str(6)))
        if h == '++':
            # b >= BA: 162u/b - 21250u/b^3 >= 0 ; FF >= 80/b^3 - C_R/b^4 + 2.5 b^2 u^2 - 7.5 u^3
            c1 = bool(fq(162) - fq(21250) / (b * b) > 0)
            c2 = bool(80 * b - CR1 > 0)                                  # u <= u1
            uu = fq(u1)
            c3 = bool(fq(Fr(5, 2)) * b ** 2 * uu ** 2 * (1 - 3 * fq(Fr(1, 4)) / b ** 2) - CRm / b ** 4 > 0)  # u >= u1
            ok = c1 and c2 and c3
            print('M ++ b>=%s: (162-21250/b^2>0) %s ; 80b - C_R(u<=u1) = %s > 0 %s ; u>=u1: %s' % (
                BA, c1, (80 * b - CR1).str(6), c2, c3))
        else:
            delta = 8 * fq(u1).sqrt() / b
            eta = fq(85) * fq(Fr(1, 4)) / b ** 2
            val = 80 * b / (1 + delta) ** 3 - CR1 / (1 - delta) ** 4 - fq(Fr(26244, 10)) / ((1 - delta) ** 2 * (1 - eta))
            c2 = bool(val > 0)
            # u in [u1,1/4]: b_- >= b-4 ; need 2.5 b^2 u^2 (1-eta) > 162 u/(b-4) + C_R/(b-4)^4  (worst at u=u1)
            bm = b - 4; uu = fq(u1)
            c3 = bool(fq(Fr(5, 2)) * b ** 2 * uu * (1 - eta) > 162 / bm) and \
                 bool(fq(Fr(5, 2)) * b ** 2 * uu ** 2 * (1 - eta) - 162 * uu / bm - CRm / bm ** 4 > 0)
            ok = c2 and c3
            print('M +- b>=%s: b^4*FF lower bound (u<=u1) = %s > 0: %s ; u in [u1,1/4]: %s' % (BA, val.str(6), c2, c3))
        ok_all &= ok
    print('VERDICT cert_M (BA=%s, u1=%s): %s' % (BA, u1, 'PASS' if ok_all else 'FAIL'))

if __name__ == '__main__':
    main()
