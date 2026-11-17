"""Certificates (FLINT/Arb, exact rationals) for a no-forward functional (NOFORWARD_FUNCTIONAL.md section 4).

F = int psi2 B2 + psi3 B3 (a = 0).  Robust margin delta': FF - delta' N >= 0, N = |B2| + |B3|, via the four sign-pattern actions of
rcert.py (gens with a = 0, so no forward term).  Regions (FLAT_GRAVITON_BOUND.md section 3 lemmas):
  P   J < JP, all u in [0,1]                 exact Lemma-P polynomials, Descartes (rcert.poly_u, rcert.positive_on_01)
  M   u <= 1/4, bhat >= BA (all J)           Lemma M: bhat^4 FF >= gamma(u) bhat/(1+ds)^3 - C_R/(1-ds)^4, gamma(u) = -(W0''(0) + u W1''(0))/2
                                              (psi3(0) = 0, so the 1/bhat term is absent; no forward term), u split into sub-intervals
  Cpp ++, J >= JP, 1/4 <= u <= 1             the same Lemma M (Legendre kernel, m in [M, 2M]; bhat >= (2JP+1)/2)
  S   J >= JP, u <= 1/4, bhat <= BA          rcert.region_S (Taylor/Cauchy cells), forward part = 0
NOT covered: +-, J >= JP, 1/4 < u <= 1 (m in [M, 2M]): Lemma C introduces (1 - u p^2)^{-8}, singular at u = 1 unless the weights vanish
to order 11 at p = 1; the forward term that dominated there for F_g3/M02 is absent.  This region is checked in floating point only.
Numerics: the weights have large monomial coefficients (cancellation), so every interval evaluation of a polynomial over a p-cell is
done in CENTERED form (exact rational Taylor shift to the cell midpoint, then Horner in the small offset): rcert.ev / rcert.poly_jet are
patched accordingly.
Usage: python3 nf_cert.py NAME DELTA REGION [args]   (REGION in P, M, Cpp, S++, S+-, low)
"""
import sys, os, json, time, hashlib
from fractions import Fraction as Fr
import flint
import rcert
from rcert import A, fq, Gen, gens

HERE = os.path.dirname(os.path.abspath(__file__))
TRACKED = os.path.join(HERE, 'functionals_noforward.json')

def load(name):
    for fn in [TRACKED, os.path.join(HERE, 'cache', 'nf', 'functionals_nf.json')]:
        if os.path.exists(fn):
            D = json.load(open(fn))
            if name in D:
                d = D[name]
                psi2 = [Fr(c) for c in d['psi2_exact']]; psi3 = [Fr(c) for c in d['psi3_exact']]
                assert Fr(d['a_exact']) == 0
                X = [Fr(c) for c in d['X']]; Y = [Fr(c) for c in d['Y']]
                from nf_mkfun import expand
                assert expand(d['r'], X) == psi2 and expand(d['r'], Y) == psi3
                canon = 'psi2=' + ','.join(map(str, psi2)) + ';psi3=' + ','.join(map(str, psi3)) + ';a=0'
                print('FUNCTIONAL %s from %s ; sha256(exact coefficients) = %s' % (name, os.path.relpath(fn, HERE),
                      hashlib.sha256(canon.encode()).hexdigest()), flush=True)
                return psi2, psi3, Fr(0)
    raise KeyError(name)

# ------------------------------------------------------------------ centered interval evaluation (patched into rcert)
_shift_cache = {}
def _mid_fmpq(x):
    m, e = x.mid().man_exp()
    m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** (-e))
def _shifted(poly, c):
    key = (tuple(poly), c)
    S = _shift_cache.get(key)
    if S is None:
        P = flint.fmpq_poly([flint.fmpq(x.numerator, x.denominator) for x in poly])
        S = P(flint.fmpq_poly([flint.fmpq(c.numerator, c.denominator), 1]))
        S = [Fr(int(S.coeffs()[i].p), int(S.coeffs()[i].q)) for i in range(S.degree() + 1)] if S.degree() >= 0 else [Fr(0)]
        if len(_shift_cache) > 200000: _shift_cache.clear()
        _shift_cache[key] = S
    return S
_orig_ev = rcert.ev
def ev_centered(poly, x):
    if not isinstance(x, flint.arb) or x.rad() == 0:
        return _orig_ev(poly, x)
    c = _mid_fmpq(x); T = x - fq(c)
    r = A(0)
    for co in reversed(_shifted(poly, c)): r = r * T + fq(co)
    return r
_orig_poly_jet = rcert.poly_jet
def poly_jet_centered(coeffs, pj):
    p0 = pj.c[0]
    if not isinstance(p0, flint.arb) or p0.rad() == 0:
        return _orig_poly_jet(coeffs, pj)
    c = _mid_fmpq(p0)
    tj = rcert.Jet([p0 - fq(c)] + pj.c[1:])
    r = rcert.Jet([A(0)])
    for co in reversed(_shifted(coeffs, c)): r = r * tj + rcert.Jet([fq(co)])
    return r
rcert.ev = ev_centered
rcert.poly_jet = poly_jet_centered

# ------------------------------------------------------------------ region M (no forward term, psi3(0) = 0)
def lemmaM(g, b, ulo, uhi, NP):
    """lower bound for bhat^4 FF over u in [ulo, uhi], bhat >= b (increasing in bhat; checked at b)."""
    assert g.W0[0] == 0 and g.W1[0] == 0
    q0 = 2 * (g.W0[2] if len(g.W0) > 2 else 0); q1 = 2 * (g.W1[2] if len(g.W1) > 2 else 0)
    gam_lo = min(-q0 / 2 - q1 * ulo / 2, -q0 / 2 - q1 * uhi / 2)         # gamma(u) linear in u
    CR = rcert.CR_bound(g, ulo, uhi, NP)
    ds = 8 * fq(uhi).sqrt() / b if g.h == '+-' else A(0)
    val = fq(gam_lo) * b / (1 + ds) ** 3 - CR / (1 - ds) ** 4
    return gam_lo, CR, val

def region_M(F, delta, BA, NP, usplit):
    flint.ctx.prec = 128
    psi2, psi3, a = F
    b = fq(BA); ok_all = True
    for h in ['++', '+-']:
        for g in gens(psi2, psi3, a, delta, h, 1):
            for ulo, uhi in zip(usplit[:-1], usplit[1:]):
                gam, CR, val = lemmaM(g, b, ulo, uhi, NP)
                ok = gam > 0 and bool(val > 0)
                print('M %s u in [%s,%s]: gamma >= %.5g  C_R <= %s  bhat^4 FF >= %s at bhat=%s : %s' % (
                    g.label, ulo, uhi, float(gam), CR.str(5), val.str(5), BA, ok), flush=True)
                ok_all &= ok
    print('VERDICT nf_cert_M (delta=%s, BA=%s): %s' % (delta, BA, 'PASS' if ok_all else 'FAIL'))

def region_Cpp(F, delta, JP, NP, usplit, h='++'):
    """h, J >= JP, u in usplit range: Lemma M (++: Legendre; +-: Lemma C with (1-up^2)^{-8}), bhat = (2J+1) sqrt(u) >= (2JP+1) sqrt(ulo)."""
    flint.ctx.prec = 128
    psi2, psi3, a = F
    ok_all = True
    for g in gens(psi2, psi3, a, delta, h, 1):
        for ulo, uhi in zip(usplit[:-1], usplit[1:]):
            b = fq(Fr(2 * JP + 1)) * fq(ulo).sqrt()
            gam, CR, val = lemmaM(g, b, ulo, min(uhi, Fr(1)), NP)
            ok = gam > 0 and bool(val > 0)
            print('Cpp %s u in [%s,%s], bhat >= %s: gamma >= %.5g C_R <= %s bhat^4 FF >= %s : %s' % (
                g.label, ulo, uhi, b.str(4), float(gam), CR.str(5), val.str(5), ok), flush=True)
            ok_all &= ok
    print('VERDICT nf_cert_C%s (delta=%s, JP=%d, u in [%s,%s]): %s' % ('pp' if h == '++' else 'pm', delta, JP, usplit[0], usplit[-1], 'PASS' if ok_all else 'FAIL'))

def low(F):
    """exact: k, A, B; matter coefficients negative on 0 < m_l <= 1 (as polynomials in x = m_l^2, exact Descartes on [0,1])."""
    psi2, psi3, _ = F
    mom = rcert.mom
    k = -2 * mom(psi2, 6) + 2 * mom(psi3, 4)
    Lfin = sum(psi2[i] / (i - 1) for i in range(2, len(psi2)))
    # C0 = mu2_2 - mu3_0 - 2 x mu2_0 ;  m^4 C2 = sum_j c_j x^j  (x = m_l^2)
    C0 = [mom(psi2, 2) - mom(psi3, 0), -2 * mom(psi2, 0)]
    C2 = [6 * mom(psi2, 6) - 6 * mom(psi3, 4),
          -6 * mom(psi2, 4) - 12 * mom(psi2, 4) + 6 * mom(psi3, 2),
          mom(psi2, 2) + 12 * mom(psi2, 2) - mom(psi3, 0),
          -2 * mom(psi2, 0)]
    ok = True
    for nm, P in [('C0', C0), ('m^4 C2', C2)]:
        Q, L = rcert.fbc.to_fmpz_poly([-c for c in P])
        pos, nint = rcert.positive_on_01(Q)
        print('low: %s (poly in m_l^2) = %s ; -poly > 0 on [0,1] (exact Descartes): %s' % (nm, [str(c) for c in P], pos)); ok &= pos
    print('low: k = %s = %.10f ; A = %.6f ; B = %.6f ; B_AdS = %.6f' % (k, float(k), float(16 * psi2[1] / k), float(16 * Lfin / k),
          float(16 * Lfin / k) - float(16 * psi2[1] / k) * 0.4054651081081644))
    print('VERDICT nf_cert_low: %s' % ('PASS' if ok else 'FAIL'))

if __name__ == '__main__':
    name, delta, region = sys.argv[1], Fr(sys.argv[2]), sys.argv[3]
    F = load(name)
    if region == 'P':
        rcert.region_P(F, delta, int(sys.argv[4]) if len(sys.argv) > 4 else 400)
    elif region == 'M':
        BA = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(120)
        region_M(F, delta, BA, 2000, [Fr(0), Fr(1, 100), Fr(1, 25), Fr(1, 10), Fr(1, 4)])
    elif region == 'Cpp':
        region_Cpp(F, delta, int(sys.argv[4]) if len(sys.argv) > 4 else 400, 2000, [Fr(1, 4), Fr(3, 8), Fr(1, 2), Fr(5, 8), Fr(3, 4), Fr(7, 8), Fr(1)])
    elif region == 'Cpm':
        ups = [Fr(x) for x in sys.argv[5].split(',')]
        region_Cpp(F, delta, int(sys.argv[4]), 2000, ups, h='+-')
    elif region.startswith('S'):
        rcert.region_S(F, delta, region[1:], Fr(1), Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(120), int(sys.argv[5]) if len(sys.argv) > 5 else 400,
                       Tmax=Fr(sys.argv[6]) if len(sys.argv) > 6 else Fr(10 ** 6))
    elif region == 'low':
        low(F)
