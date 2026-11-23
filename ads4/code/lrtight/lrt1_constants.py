"""
lrtight/lrt1_constants.py (LRINF_TIGHT.md §2.1): certified constants of Lemma LR-inf.

  * int_0^1 |psi_k| dp, k = 2, 3 [certified]: the sign changes of psi_k on (0, 1) are counted EXACTLY with a Sturm sequence over Q
    (sympy, rational arithmetic); each root is bracketed by exact rational bisection (sign of psi at rationals, exact);
    int|psi| = sum over sign-constant pieces of |Q(e_{i+1}) - Q(e_i)|, Q the exact antiderivative, evaluated in Arb at the root
    bracket balls.  Cross-check: rigorous Riemann upper/lower sums of |psi| on 4096 cells (Arb).
  * a1 = psi'(0), a2 = psi''(0), A1min = 2 a1^(2) - a1^(3), A2' = -2 a2^(2) - |a2^(3)|  (exact rationals).
  * the crude constant Phi_max = sum |c_n|/(n+1) of lsl_lrinf.py, for comparison.
Output: cache/lrt1_constants.json, cache/lrt1_constants.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import sympy as sp
from flint import arb
from lrt_common import Weight, A, CACHE, up, lo

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)


def peval(c, x):
    s = Fr(0)
    for cn in reversed(c):
        s = s*x + cn
    return s


def certified_absint(W, nbis=200):
    c = W.c
    p = sp.symbols('p')
    P = sp.Poly(sum(sp.Rational(cn.numerator, cn.denominator)*p**n for n, cn in enumerate(c)), p, domain='QQ')
    nroots = P.count_roots(sp.Rational(1, 10**12), 1 - sp.Rational(1, 10**12))    # exact (Sturm); ends: psi(0)=0, (1-p)^5
    # psi has a simple zero at 0 and a 5-fold zero at 1; check there are no other roots in (0, 1e-12] or [1-1e-12, 1):
    # psi(p)/p and psi(p)/(1-p)^5 are nonzero at the endpoints and continuous; verify exactly that no root lies in the end slivers
    Pq = sp.quo(P, sp.Poly(p*(1 - p)**5, p, domain='QQ'))
    assert sp.rem(P, sp.Poly(p*(1 - p)**5, p, domain='QQ')).is_zero
    n_all = Pq.count_roots(0, 1)                                                  # roots of X_k on [0, 1] (exact)
    # isolate with sympy's exact real-root isolation (rational intervals)
    ivs = [iv for iv in Pq.intervals() if iv[0][1] > 0 and iv[0][0] < 1]
    roots = []
    for (a, b), mult in ivs:
        a = Fr(int(a.p), int(a.q)); b = Fr(int(b.p), int(b.q))
        a = max(a, Fr(0)); b = min(b, Fr(1))
        cq = [Fr(int(v.p), int(v.q)) for v in reversed(Pq.all_coeffs())]    # X_k (psi / (p (1-p)^5)), nonzero at 0, 1
        c_save = c; c = cq
        fa, fb = peval(c, a), peval(c, b)
        if a == b:
            roots.append((a, b)); continue
        assert fa*fb < 0, 'bracket'
        for _ in range(nbis):
            mid = (a + b)/2; fm = peval(c, mid)
            if fm == 0:
                a = b = mid; break
            if fa*fm < 0:
                b, fb = mid, fm
            else:
                a, fa = mid, fm
            if b - a < Fr(1, 2**120):
                break
        roots.append((a, b))
        c = c_save
    assert len(roots) == n_all == nroots, (len(roots), n_all, nroots)
    Q = [Fr(0)] + [cn/(n + 1) for n, cn in enumerate(c)]           # exact antiderivative
    def Qb(x):
        s = arb(0)
        for qn in reversed(Q):
            s = s*x + A(qn)
        return s
    pts = [A(0)] + [A(a).union(A(b)) for a, b in sorted(roots)] + [A(1)]
    tot = arb(0)
    for i in range(len(pts) - 1):
        tot += abs(Qb(pts[i + 1]) - Qb(pts[i]))
    signed = Qb(A(1)) - Qb(A(0))
    return tot, roots, signed


def riemann_check(W, N=4096):
    """rigorous Riemann bounds: sum_cells |cell| * sup|psi(cell)| (sup from the exact Taylor expansion at the cell centre) >= int|psi| >= sum_cells |int_cell psi|."""
    c = W.c
    upper = arb(0); lower = arb(0)
    Q = [Fr(0)] + [cn/(n + 1) for n, cn in enumerate(c)]
    def ev(cc, x):
        s = arb(0)
        for qn in reversed(cc):
            s = s*x + A(qn)
        return s
    from flint import fmpq_poly, fmpq
    P = fmpq_poly([fmpq(cn.numerator, cn.denominator) for cn in c])
    h = Fr(1, 2*N)
    for i in range(N):
        m = Fr(2*i + 1, 2*N)
        T = P(fmpq_poly([fmpq(m.numerator, m.denominator), 1]))            # exact Taylor expansion about the cell centre
        sup = sum(abs(Fr(int(T[j].p), int(T[j].q)))*h**j for j in range(T.degree() + 1))
        upper += A(sup)/N
        lower += abs(ev(Q, A(Fr(i + 1, N))) - ev(Q, A(Fr(i, N))))
    return upper, lower


def main():
    pr(__doc__.split('Output')[0])
    out = {}
    W = {k: Weight(k) for k in (2, 3)}
    for k in (2, 3):
        w = W[k]
        tot, roots, signed = certified_absint(w)
        upR, loR = riemann_check(w)
        crude = sum(abs(cn)/(n + 1) for n, cn in enumerate(w.c))
        pr(f'k = {k}: interior sign changes of psi_{k} on (0,1): {len(roots)} at p in '
           + ', '.join(f'[{float(a):.15f}, {float(b):.15f}]' for a, b in roots))
        pr(f'   int_0^1 |psi_{k}| = {tot.str(20, radius=True)}   [certified: Sturm + exact bisection + Arb]')
        pr(f'   Riemann cross-check (4096 cells): {lo(loR).str(12)} <= int|psi| <= {up(upR).str(12)}')
        pr(f'   int_0^1 psi_{k} = {signed.str(15)};  crude Phi_max = sum|c_n|/(n+1) = {float(crude):.6e}  (ratio crude/int|psi| = {float(crude)/float(tot.mid()):.1f})')
        pr(f'   a1 = psi\'(0) = {float(w.a1):.10f}   a2 = psi\'\'(0) = {float(w.a2):.10e}   (exact rationals)')
        out[k] = dict(absint_mid=float(tot.mid()), absint_rad=float(tot.rad()), crude=float(crude),
                      a1=str(w.a1), a2=str(w.a2), roots=[[str(a), str(b)] for a, b in roots])
        # store a rational upper bound (safe): round the ball's upper endpoint up at 1e-15 relative
        u = float(up(tot).mid())
        ub = Fr(u) * (1 + Fr(1, 10**12))
        assert A(ub) > tot
        out[k]['absint_upper_rational'] = str(ub)
        pr(f'   rational upper bound stored: {float(ub):.15f}')
    A1min = 2*W[2].a1 - W[3].a1
    A2p = -2*W[2].a2 - abs(W[3].a2)
    pr(f'A1min = 2 a1(2) - a1(3) = {float(A1min):.10f};  A2\' = -2 a2(2) - |a2(3)| = {float(A2p):.10f}  (exact rationals)')
    out['A1min'] = str(A1min); out['A2p'] = str(A2p)
    json.dump(out, open(os.path.join(CACHE, 'lrt1_constants.json'), 'w'), indent=1)
    open(os.path.join(CACHE, 'lrt1_constants.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
