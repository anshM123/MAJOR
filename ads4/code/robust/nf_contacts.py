"""NOFORWARD_FUNCTIONAL.md section 1: can smeared B^{(1)}_k, k >= 4 (no p = 0 terms) be combined so that ALL contact couplings cancel?

Exact computation (sympy residues for the low sides, FLINT rational linear algebra for the null space).
  * Low side of a contact monomial f = t^a (s u)^b (degree n = a + 2b; every symmetric polynomial f(s,u) = f(u,s) is a combination of
    these) in -B^{(1)}_k: (1/2) sum of residues at s = 0, p^2 of the integrand of (2.26), exactly as flatbound/cert_low.py.  It equals
    c_{k}(a,b) p^{2(n-k+4)} (homogeneity), and vanishes identically for k = 2, 3 (no poles: B_2, B_3 are contact-free).
  * Check of the closed form [D, section 1.2 of the .md]: even k: c = coefficient of z^{k/2-2} of h(s) + h(q-s), odd k: coefficient of
    z^{(k-5)/2} of (h(q-s) - h(s))/(2s-q), with h(s) = f(s,-q), q = p^2, z = s(s-q)  (up to the overall factor 1).
  * Null space: weights psi_k(p) = sum_{i<=D} c_{k,i} p^i (k = 4..K, arbitrary polynomials, no vanishing conditions imposed);
    the functional sum_k int_0^1 psi_k(-B_k|_contact) dp vanishes for the monomials of degree <= N iff the rational matrix
    R[(a,b), (k,i)] = c_k(a,b) / (i + 2(n-k+4) + 1) annihilates the coefficient vector.  rank = number of unknowns => only psi = 0.
  * Truncated EFT (only g4, g5, g6, g6' = degrees n <= 2, as in (2.7a)): null-space dimension printed for contrast.
  * Structure check of the proof's induction (top-order data enter the pair (B_{2l+4}, B_{2l+5}) with opposite-sign ratios).
Usage: python3 nf_contacts.py [K] [D] [N]
"""
import sys
import sympy as sp
import flint

s, p, q, z = sp.symbols('s p q z')

def minus_B_low_contact(k, a, b):
    T = lambda S, U: (-S - U)
    f = lambda S, U: T(S, U) ** a * (S * U) ** b
    br_even = f(s, -p**2) + f(p**2 - s, -p**2)
    br_odd = -f(s, -p**2) + f(p**2 - s, -p**2)
    if k % 2 == 0:
        e = (2*s - p**2) / (s*(s - p**2))**sp.Rational(k - 2, 2) * br_even
    else:
        e = 1 / (s*(s - p**2))**sp.Rational(k - 3, 2) * br_odd
    return sp.expand(sum(sp.residue(e, s, zz) for zz in [0, p**2]) / 2)

def closed_form(k, a, b):
    h = lambda S: (-(S) - (-q)) ** a * (S * (-q)) ** b         # f(S, -q) with t = -S - u
    if k % 2 == 0:
        F = sp.expand(h(s) + h(q - s)); ell = k // 2 - 2
    else:
        F = sp.cancel(sp.expand(h(q - s) - h(s)) / (2*s - q)); ell = (k - 5) // 2
    # F is a polynomial in s symmetric under s -> q - s, hence a polynomial in z = s(s - q)
    F = sp.Poly(F, s)
    # expand in z by repeated division
    coeffs = []
    rem = F.as_expr()
    for j in range(0, 40):
        c0 = sp.expand(rem.subs(s, 0))
        coeffs.append(c0)
        rem = sp.expand(sp.cancel((rem - c0) / (s*(s - q))))
        if rem == 0:
            break
    return sp.expand(coeffs[ell] if ell < len(coeffs) else 0)

def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    D = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    ok = True
    # (i) B2, B3 contact-free; closed form check (small degrees)
    for k in [2, 3]:
        for (a, b) in [(0, 0), (1, 0), (2, 0), (0, 1), (3, 1), (1, 2)]:
            ok &= minus_B_low_contact(k, a, b) == 0
    print('B_2, B_3 receive no contact contribution (degrees <= 5 checked):', ok)
    cf_ok = True
    for k in range(4, K + 1):
        ratios = set()
        for n in range(0, 7):
            for b in range(0, n // 2 + 1):
                a = n - 2*b
                r = minus_B_low_contact(k, a, b)
                c = closed_form(k, a, b).subs(q, p**2)
                if c == 0:
                    cf_ok &= (r == 0)
                else:
                    ratios.add(sp.simplify(r / c))
        cf_ok &= len(ratios) == 1 and list(ratios)[0] != 0
        print('  k=%d: residue / closed form = %s (same for all monomials of degree <= 6)' % (k, ratios))
    print('closed form (z-expansion coefficient) agrees with the residues up to a nonzero constant per k (n <= 6):', cf_ok)
    # (i') structure used in the proof (NOFORWARD_FUNCTIONAL.md section 1.2): the order-l data, isolated by P = w^l t^L (s = 0 family)
    #      and P = t^l w^L (s = q family), L = 8, enter only B_k with k >= 2l+4, and in the pair (B_{2l+4}, B_{2l+5}) with ratios
    #      B_{2l+5}/B_{2l+4} of opposite sign (-1/q vs +1/q up to a common factor): nonsingular 2x2 system.
    st_ok = True
    Lh = 8
    for l in range(0, 3):
        rat = []
        for (a, b) in [(Lh, l), (l, Lh)]:
            vals = {k: closed_form(k, a, b) for k in range(4, 2*l + 6)}
            st_ok &= all(vals[k] == 0 for k in range(4, 2*l + 4))
            rat.append(sp.simplify(vals[2*l + 5] / vals[2*l + 4]) if vals[2*l + 4] != 0 else None)
        st_ok &= rat[0] is not None and rat[1] is not None and sp.simplify(rat[0] + rat[1]) == 0
        print('  order l=%d: zero in B_4..B_%d; ratio B_%d/B_%d = %s (s=0 family), %s (s=q family)' % (l, 2*l + 3, 2*l + 5, 2*l + 4, rat[0], rat[1]))
    print('proof structure (top-order data, opposite-sign pair) verified for l <= 2:', st_ok)
    ok &= st_ok
    # (ii) contact data c_k(a,b): coefficient of p^{2(n-k+4)}
    rows = []; labels = []
    for n in range(0, N + 1):
        for b in range(0, n // 2 + 1):
            a = n - 2*b
            row = []
            for k in range(4, K + 1):
                e = n - k + 4
                if e < 0:
                    ck = sp.Integer(0)
                else:
                    r = closed_form(k, a, b).subs(q, p**2)
                    ck = sp.Poly(sp.expand(r), p).coeff_monomial(p**(2*e)) if r != 0 else sp.Integer(0)
                    assert sp.expand(r - ck * p**(2*e)) == 0
                for i in range(0, D + 1):
                    row.append(sp.Rational(ck) / (i + 2*e + 1) if e >= 0 else sp.Integer(0))
            rows.append(row); labels.append((n, a, b))
    nun = (K - 3) * (D + 1)
    def rank_upto(nmax):
        R = [r for r, l in zip(rows, labels) if l[0] <= nmax]
        M = flint.fmpq_mat(len(R), nun, [flint.fmpq(int(x.p), int(x.q)) for r in R for x in r])
        return len(R), M.rank()
    for nmax in sorted(set([2, 4, 8, 12, 16, 20, 24, 28, N])):
        nr, rk = rank_upto(nmax)
        print('contacts of degree n <= %2d (%3d monomials): rank %3d of %3d unknowns (psi_4..psi_%d of degree <= %d) -> null space dim %d'
              % (nmax, nr, rk, nun, K, D, nun - rk), flush=True)
    nr, rk = rank_upto(N)
    full = rk == nun
    print('VERDICT nf_contacts (K=%d, D=%d, N=%d): %s' % (K, D, N, 'PASS: only psi_k = 0 (k >= 4) cancels all contacts' if (ok and cf_ok and full) else 'FAIL/INCONCLUSIVE'))

if __name__ == '__main__':
    main()
