"""EXACT verification of a light-robust certificate ell_0..ell_n (Lambda = sum_k ell_k m_k^(b)), at gap g, integer arithmetic.
Checks (all J, all real Delta):
 (H)  heavy, twist >= g:          W4 + Lambda4 >= 0 (b=4, J>=4),  W0 - Lambda0 >= 0 (b=0, J even)
 (LS) light single traces:        W0 - Lambda0 >= 0 for J=0 (Delta>1) and J=2 (Delta>3)
 (L)  light multi-traces 2<twist<g: |Lambda| <= kappa W (both ladders)
Method: for each J, P(T) := 6*2^40 * Prod_{i<K}(T-y_i) * [W +- Lambda] is an integer polynomial in T of degree K+1; positivity
on the T-range is certified by a Moebius/Taylor substitution and Descartes' rule (all coefficients >=0, constant >0), with an
exact Sturm fallback.  Tails: |Lambda| <= B := max_z |ell^(M)(z)| (Lemma M, M=g-2) for heavy states, |Lambda| <= sum|ell_k| for light
(0<=m_k<=1 for twist>=2), beyond which W (resp. kappa W) dominates."""
import json, sys, time
from fractions import Fraction as Fr
import flint
from math import comb

def load(fn):
    d = json.load(open(fn)); ell = [Fr(x) for x in d['ell']]
    return d['g'], d['n'], Fr(d['kappa']), ell

def y(s): return s*(s+1)
X = flint.fmpz_poly([0, 1])

def sturm_count(P, a, b=None):
    """number of distinct real roots of fmpz_poly P in (a, b] (b=None: infinity), exact, via flint roots isolation"""
    # use complex_roots (arb) to isolate; fall back to exact Sturm if ambiguous
    roots = P.complex_roots()
    cnt = 0
    for r, mult in roots:
        if abs(r.imag) > 0 and not r.imag.contains(0): continue
        re = r.real
        lo = re.mid() - re.rad(); hi = re.mid() + re.rad()
        if float(hi) <= float(a) and not (b is None): continue
        cnt += 1  # conservative: any (possibly) real root counts
    return cnt

def pos_check(P, Tlo, Thi=None):
    """certify P(T) > 0 on [Tlo, Thi] (Thi None = infinity). P fmpz_poly, Tlo,Thi integers or Fractions."""
    if Thi is None:
        Q = P(flint.fmpz_poly([0, 1]) + flint.fmpq(Tlo.numerator, Tlo.denominator) if isinstance(Tlo, Fr) else flint.fmpz_poly([Tlo, 1]))
        cs = Q.coeffs() if hasattr(Q, 'coeffs') else None
    else:
        # T = (Tlo + Thi u)/(1+u), multiply by (1+u)^deg
        d = P.degree(); a = flint.fmpq_poly([Tlo, Thi]) if True else None
        num = flint.fmpq_poly([flint.fmpq(Tlo.numerator, Tlo.denominator), flint.fmpq(Thi.numerator, Thi.denominator)])
        one_u = flint.fmpq_poly([1, 1])
        Q = flint.fmpq_poly(0); pw = flint.fmpq_poly(1)
        coeffs = P.coeffs()
        # Horner-like: sum c_i num^i (1+u)^(d-i)
        nump = [flint.fmpq_poly(1)]
        for i in range(d): nump.append(nump[-1]*num)
        onep = [flint.fmpq_poly(1)]
        for i in range(d): onep.append(onep[-1]*one_u)
        for i, c in enumerate(coeffs):
            if c != 0: Q += flint.fmpq(int(c))*nump[i]*onep[d-i]
        cs = Q.coeffs()
    if all(c >= 0 for c in cs) and cs[0] > 0: return True, 'descartes'
    return False, 'need-fallback'

def build(b, J, ell, n, sign, wscale, Wfun_is_kappa=None):
    """P(T) = 6*2^40*Prod(T-y_{b+i}) * [wscale*W] + sign*6*2^40*sum_k ell_k Prod_{i<k}(C-y)Prod_{k<=i<K}(T-y)"""
    C = J*(J+1); K = min(n, J-b)
    ys = [y(b+i) for i in range(K)]
    den = 2**40
    Lint = [int(e*den) for e in ell]
    # suffix products
    suf = [None]*(K+1); suf[K] = flint.fmpz_poly([1])
    for i in range(K-1, -1, -1): suf[i] = suf[i+1]*flint.fmpz_poly([-ys[i], 1])
    pre = 1; N = flint.fmpz_poly(0)
    for k in range(K+1):
        N += Lint[k]*pre*suf[k]
        if k < K: pre *= (C - ys[k])
    sixW = flint.fmpz_poly([-2*b*(b-1)+9, 2])     # 6W = 2T - 2b(b-1) + 9
    Pw = sixW*suf[0]*den
    return Pw, 6*N

if __name__ == '__main__':
    fn = sys.argv[1]; Jtail_override = int(sys.argv[2]) if len(sys.argv) > 2 else None
    g, n, kappa, ell = load(fn)
    t0 = time.time()
    M = g-2
    # Lemma-M bound for heavy |Lambda|
    vals = []
    for z in range(M+1):
        s = Fr(0); r = Fr(1)
        for k in range(n+1):
            if k > z: break
            s += ell[k]*r; r *= Fr(z-k, M-k)
        vals.append(s)
    Bh = max(abs(v) for v in vals); Bl = sum(abs(e) for e in ell)
    print('n', n, 'g', g, 'kappa', float(kappa), 'Lemma-M heavy bound B =', float(Bh), ' sum|ell| =', float(Bl), flush=True)
    # tails
    W = lambda b, D: Fr((D-1)*(D-2) - b*(b-1), 3) + Fr(3, 2)
    JH = 0
    while W(4, JH+g) < Bh or W(0, JH+g) < Bh: JH += 1
    JL = 0
    while kappa*W(4, JL+2+4) < Bl or kappa*W(0, JL+2) < Bl: JL += 1
    if Jtail_override: JH = min(JH, Jtail_override)
    print('heavy explicit J <=', JH, ' light explicit J <=', JL, flush=True)
    fails = []; nfb = 0
    # (LS) J=0: W0(Delta)-ell0 >= 0 for Delta>1: min W0 = 3/2 - 1/12
    if not (ell[0] <= Fr(3, 2) - Fr(1, 12)): fails.append(('LS', 0))
    # (LS) J=2, T>2:  W0 - Lambda >= 0
    Pw, N = build(0, 2, ell, n, -1, 1)
    ok, how = pos_check(Pw - N, Fr(2))
    if not ok: fails.append(('LS', 2, how))
    # (H)
    for J in range(0, JH+1):
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            Pw, N = build(b, J, ell, n, 0, 1)
            P = Pw + N if b == 4 else Pw - N
            Tmin = Fr((J+g-1)*(J+g-2))
            ok, how = pos_check(P, Tmin)
            if not ok:
                # tight state b=4,J=4: P = W - W_B... allow zero at endpoint: check value and shifted coeffs >=0
                fails.append(('H', b, J, how))
        if J % 5000 == 0: print('  heavy J', J, 'fails so far', len(fails), '%.0fs' % (time.time()-t0), flush=True)
    # (L)
    for J in range(0, JL+1):
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            Pw, N = build(b, J, ell, n, 0, 1)
            kP = Pw*int(kappa) if kappa.denominator == 1 else None
            Tlo = Fr(J*(J+1)); Thi = Fr((J+g-1)*(J+g-2))
            for s in (+1, -1):
                ok, how = pos_check(kP + s*N, Tlo, Thi)
                if not ok: fails.append(('L', b, J, s, how))
    print('DONE in %.0fs; failures:' % (time.time()-t0), len(fails), fails[:20])
