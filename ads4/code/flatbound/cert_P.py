"""Certificate P: for every spin J < JP (both channels) the exact polynomial FF_h(J,u) (degree <= J+1 in u,
rational coefficients, pessimistic Phi) is > 0 on the whole closed interval u in [0,1], i.e. for all m >= 1.

Method (exact integer arithmetic, no floating point):
  * FF(J,0) > 0 is checked exactly;
  * [0,1] is covered by dyadic intervals I=[k/2^s,(k+1)/2^s]; on each, Q(y)=2^{sd} P((k+y)/2^s) and
    R(x)=(1+x)^d Q(1/(1+x)); if R has no sign variation, Q has no root in (0,1) (Descartes' rule of signs);
    the endpoint values are checked to be nonzero exactly; otherwise I is bisected (Vincent/Collins-Akritas).
Usage: python3 cert_P.py JP [channel]
"""
import sys, time
import flint
from common import poly_u, to_fmpz_poly

def sign_var(coeffs):
    s = [c for c in coeffs if c != 0]
    return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))

def interval_poly(P, k, s):
    d = P.degree(); c = P.coeffs()
    x = flint.fmpz_poly([k, 1]); Q = flint.fmpz_poly([c[d]])
    for j in range(d - 1, -1, -1):
        Q = Q * x + flint.fmpz_poly([c[j] * (2 ** (s * (d - j)))])
    return Q

def descartes01(Q):
    d = Q.degree(); c = list(Q.coeffs()) + [0] * (d + 1 - len(Q.coeffs()))
    R = flint.fmpz_poly(list(reversed(c)))(flint.fmpz_poly([1, 1]))
    return sign_var(R.coeffs())

def positive_on_01(P, maxdepth=60):
    if P(0) <= 0:
        return False, 0
    stack = [(0, 0)]; n = 0
    while stack:
        k, s = stack.pop(); n += 1
        Q = interval_poly(P, k, s)
        if Q(0) == 0 or Q(1) == 0:
            return False, n
        if descartes01(Q) == 0:
            continue
        if s >= maxdepth:
            return False, n
        stack += [(2 * k, s + 1), (2 * k + 1, s + 1)]
    return True, n

def main():
    JP = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    chans = [sys.argv[2]] if len(sys.argv) > 2 else ['++', '+-']
    allok = True; t0 = time.time(); cnt = 0
    for h in chans:
        Js = range(0, JP, 2) if h == '++' else range(4, JP)
        for J in Js:
            P, L = to_fmpz_poly(poly_u(J, h))
            ok, n = positive_on_01(P)
            cnt += 1
            print('P %s J=%d deg=%d intervals=%d %s' % (h, J, P.degree(), n, 'PASS' if ok else 'FAIL'), flush=True)
            allok &= ok
    print('VERDICT cert_P (JP=%d, %d polynomials, %.0fs): %s' % (JP, cnt, time.time() - t0, 'PASS' if allok else 'FAIL'))

if __name__ == '__main__':
    main()
