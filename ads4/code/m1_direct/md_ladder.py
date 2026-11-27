"""md_ladder.py -- the ladder (fp-H) per-nu image polynomials, exact in sympy, for the comparisons of m1_direct.
Implements Eqs. (fp-U), (fp-zy), (fp-H) of prl/sec_fullpos.tex from the TeX (own code; nothing imported from the project).
per-nu k=2 image of state (Delta, J) on ladder b, up to the state-independent nu-factor g(nu) and the nu-independent B|a|^2:
   h2_b(nu; Delta, J) = (1-r)^3 sum_{k<=n} D_k^2 U_{b,k} [ z_{k+2}(nu) (1/(k+3) - r/(k+1)) + (2 y_{k+2}(nu) + V_{b,k} z_{k+2}(nu))/T ]
"""
import sympy as sp

s = sp.Symbol('s')   # s = nu^2 + 1/4


def fseq(M):
    f = [sp.Integer(0), sp.Integer(1)]
    for m in range(1, M):
        f.append(sp.expand((-(s - 2*m*m)*f[m] - (m - 1)**2*f[m - 1])/sp.Integer((m + 1)**2)))
    return f


def zy(M):
    f = fseq(M + 4)
    F = lambda m: f[m] if m >= 0 else sp.Integer(0)
    z = {}; y = {}
    for m in range(0, M + 1):
        z[m] = sp.expand(12*F(m) - (m - 1)*F(m - 1) + (m + 1)*F(m + 1))
        y[m] = sp.expand(sp.Rational(12 - m*m, 2)*F(m) + sp.Rational(13 - 7*m, 2)*F(m - 1) + sp.Rational(13 + 7*m, 2)*F(m + 1)
                         + sp.Rational((m - 2)**2, 4)*F(m - 2) + sp.Rational((m + 2)**2, 4)*F(m + 2))
    return z, y


def U(b, k, Delta, J):
    out = sp.Integer(1)
    T = (Delta - 1)*(Delta - 2); C = J*(J + 1)
    for j in range(k):
        u = (b + j)*(b + j + 1)
        out *= sp.Rational(1)*(C - u)/(T - u)
    return out


def h2(b, Delta, J, z, y, stripped=True):
    """per-nu k=2 ladder polynomial in s. stripped=True returns hhat = T^4/m^6 * (1-r)^3 sum ... /T^0 ... i.e. the ladder sum
    WITHOUT the (1-r)^3 factor times T (the 'stripped' normalisation of NOTE.md (3.2) uses T*sum); here we return S = sum_k [...]
    so that H = (1-r)^3 S, and the state factor B (1-r)^3 = (2 pi G L^4/9) T is nu-independent."""
    Delta = sp.nsimplify(Delta); T = (Delta - 1)*(Delta - 2); C = J*(J + 1); r = sp.Rational(1)*C/T
    n = J - b
    S = sp.Integer(0)
    for k in range(0, n + 1):
        Dk = sp.binomial(k + 3, 3)
        V = sp.Rational(4*k, k + 3) if b == 0 else 4 + sp.Rational(12, k + 1)
        S += Dk**2*U(b, k, Delta, J)*(z[k + 2]*(sp.Rational(1, k + 3) - r/(k + 1)) + (2*y[k + 2] + V*z[k + 2])/T)
    return sp.expand(S)


if __name__ == '__main__':
    z, y = zy(6)
    print('z_2 =', sp.factor(z[2]))
    print('y_2 =', sp.factor(y[2]))
    print('c(s) = 6 y2/z2 =', sp.simplify(6*y[2]/z[2]))
    nu = sp.Symbol('nu')
    for v in (1, 2, 4, 8):
        val = (6*y[2]/z[2]).subs(s, sp.Rational(v*v) + sp.Rational(1, 4))
        print(v, float(val))
