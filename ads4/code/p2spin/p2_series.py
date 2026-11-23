"""
p2spin/p2_series.py: truncated Laurent series in eps = 1/m with sympy coefficients (exact), used for the heavy symbols (§2).
A series is (v, [c0, c1, ..., cK]) meaning eps^v (c0 + c1 eps + ... + cK eps^K) + O(eps^{v+K+1}).
"""
import sympy as sp

EPS = sp.Symbol('eps')
ORDER = 2          # keep relative orders eps^0 .. eps^ORDER


class TS:
    __slots__ = ('v', 'c')

    def __init__(self, v, c):
        self.v = v
        self.c = list(c)[:ORDER + 1] + [sp.Integer(0)]*max(0, ORDER + 1 - len(c))

    @staticmethod
    def const(a):
        return TS(0, [sp.sympify(a)])

    @staticmethod
    def from_expr(expr, simp=sp.expand):
        """exact Laurent expansion of a rational function of EPS (other symbols as coefficients)."""
        num, den = sp.fraction(sp.together(expr))
        pn = sp.Poly(sp.expand(num), EPS); pd = sp.Poly(sp.expand(den), EPS)
        cn = pn.all_coeffs()[::-1]; cd = pd.all_coeffs()[::-1]
        vn = next(i for i, a in enumerate(cn) if a != 0); vd = next(i for i, a in enumerate(cd) if a != 0)
        N = TS(vn, cn[vn:vn + ORDER + 1]); D = TS(vd, cd[vd:vd + ORDER + 1])
        r = N*D.inv()
        return r.simplify(simp)

    def simplify(self, f=sp.expand):
        return TS(self.v, [f(a) for a in self.c])

    def __mul__(self, o):
        if not isinstance(o, TS):
            o = TS.const(o)
        c = [sp.Integer(0)]*(ORDER + 1)
        for i in range(ORDER + 1):
            for j in range(ORDER + 1 - i):
                c[i + j] += self.c[i]*o.c[j]
        return TS(self.v + o.v, c)
    __rmul__ = __mul__

    def __add__(self, o):
        if not isinstance(o, TS):
            o = TS.const(o)
        v = min(self.v, o.v)
        c = [sp.Integer(0)]*(ORDER + 1)
        for s in (self, o):
            for i, a in enumerate(s.c):
                k = i + s.v - v
                if k <= ORDER:
                    c[k] += a
        # NOTE: if the two valuations differ, the higher one loses its top orders beyond ORDER (truncation is consistent)
        return TS(v, c)
    __radd__ = __add__

    def __neg__(self):
        return TS(self.v, [-a for a in self.c])

    def __sub__(self, o):
        return self + (-o)

    def inv(self):
        a0 = self.c[0]
        assert a0 != 0
        r = [sp.Integer(1)/a0]
        for k in range(1, ORDER + 1):
            s = sum(self.c[j]*r[k - j] for j in range(1, k + 1))
            r.append(sp.expand(-s/a0))
        return TS(-self.v, r)

    def __pow__(self, n):
        if n < 0:
            return self.inv()**(-n)
        r = TS.const(1)
        for _ in range(n):
            r = r*self
        return r

    def sqrt(self, r0):
        """square root with given leading coefficient r0 (r0**2 == c0, valuation even)."""
        assert self.v % 2 == 0 and sp.simplify(r0**2 - self.c[0]) == 0
        r = [r0]
        for k in range(1, ORDER + 1):
            s = self.c[k] - sum(r[j]*r[k - j] for j in range(1, k))
            r.append(sp.expand(s/(2*r0)))
        return TS(self.v//2, r)

    def lead_normalised(self):
        """(v, [1, c1/c0, c2/c0])"""
        return self.v, [sp.cancel(a/self.c[0]) for a in self.c]
