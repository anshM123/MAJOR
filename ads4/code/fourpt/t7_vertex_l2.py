"""
<TTTT> (l = 2) at the level of the exact shadow integrand (no blocks): heavy limit of the doubly-spinning integrand
   T(z1,z2',z3,z4') = s13^3 s24^3 ((-1)^J/(4 b^(0))) L[v_a](1,2,5; Delta,J) L[v_b](3,4,5; 3-Delta,-1-J)
for Weyl-covariant conserved families v_a (shadow map = (1/b^(0)) x Weyl reflection), built numerically by the Hilbert-90 recipe
   v := c + T(Delta~) c(Delta~),   c(Delta,J) := the conserved vector with fixed values of three word coefficients (a rational family),
so that T(Delta) v(Delta) = v(Delta~) exactly.  Two different normalisation rules give two independent families v_1, v_2.
Checks at fixed (X, Y) = (m x, m y), beta fixed, m = 50 ... 5000:
 (i) helicity content: only (z1.eps_{h1})^2 (z2'.eps'_{h2})^2 with the four (h1,h2), i.e. 'same' (h2 = -h1) and 'opposite' (h2 = h1);
     longitudinal / mixed components -> 0 like 1/m;
 (ii) factorisation T^{ab}_h = V_a(h1,h2) V_b(h3,h4) and parity V(h1,h2) = V(-h1,-h2);
 (iii) Weyl antisymmetry of the O(1/m) term: the (a,b,h) <-> (b,a,h') symmetrised ratio approaches its limit like 1/m^2.
"""
import os, sys, functools, random
import mpmath as mp
import t6_structures_l2 as T6
from t6_structures_l2 import WS, FQ, LM, BR, cons_space
mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from formulas import log_b                     # noqa: E402
from vertex2 import frame, eps, inv, mdot, msub
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

RULES = [('AACC', 'HH', 'BBEE'), ('BBEE', 'HH', 'AACC')]      # (word set to 1, two words set to 0)

import numpy as np
_ns = {}
GAPS = []
def family(D, Jv, rule):
    key = (mp.nstr(D, 20), int(Jv))
    if key not in _ns:
        M = T6.divergence_rows_n(float(D), int(Jv), [float(v) for v in T6.X1], [float(v) for v in T6.X2], [float(v) for v in T6.X5],
                                 [float(v) for v in T6.N5], [[float(v) for v in z] for z in T6.OTHERS])
        sc = np.abs(M).max(axis=0)
        U, sv, Vt = np.linalg.svd(M/sc)
        GAPS.append(float(sv[-3]/sv[-4]))                      # the three smallest: 2 conserved + the 3d relation
        _ns[key] = [Vt[-k]/sc for k in (1, 2, 3)]
    ns = _ns[key]
    i1, i0a, i0b = (WS.index(w) for w in rule)
    A = np.array([[v[i1] for v in ns], [v[i0a] for v in ns], [v[i0b] for v in ns]])
    coef = np.linalg.solve(A, np.array([1.0, 0, 0]))
    return [mp.mpf(float(sum(coef[k]*ns[k][i] for k in range(3)))) for i in range(len(WS))]

def weyl_family(D, Jv, rule):
    c = family(D, Jv, rule); ct = family(3 - D, -1 - Jv, rule)
    # T(Delta~) maps coefficients at Delta~ to Delta:  c_w(Delta~) b^(0)/b_w(Delta~)
    v = [c[i] + ct[i]/BR[w](3 - D, -1 - Jv) for i, w in enumerate(WS)]
    vt = [ct[i] + c[i]/BR[w](D, Jv) for i, w in enumerate(WS)]             # = v at the reflected labels
    return v, vt

def struct(coeffs, D, Jv, xa, xb, x5, za, zb, n):
    a = T6.inv_args(xa, xb, x5, n, za, zb, D, Jv)
    return LM(*a), sum(c*FQ[w](*a) for c, w in zip(coeffs, WS))

def integrand(cL, cRt, D, Jv, x, y, z1, z2p, z3, z4p):
    e = [mp.mpf(1), 0, 0]; n = [mp.mpf(1), mp.mpf(1), 0]
    x2n = mdot(x, x)
    x1 = [-v for v in y]; x2 = [v/x2n for v in x]; x3 = list(y); x4 = [-v/x2n for v in x]
    z2 = inv(x, z2p); z4 = inv(x, z4p)
    lL, pL = struct(cL, D, Jv, x1, x2, e, z1, z2, n)
    lR, pR = struct(cRt, 3 - D, -1 - Jv, x3, x4, e, z3, z4, n)
    s13 = -mdot(msub(x1, x3), msub(x1, x3)); s24 = -mdot(msub(x2, x4), msub(x2, x4))
    return lL + lR + 3*(mp.log(s13) + mp.log(s24)) - mp.log(4) - log_b(D, Jv, 0), pL*pR      # sign (-1)^J dropped (J non-integer allowed)

def helicity(cL, cRt, D, Jv, X, Y):
    Ap, Am = D - mp.mpf(3)/2, Jv + mp.mpf(1)/2
    m = mp.sqrt(Ap*Ap - Am*Am)
    x = [v/m for v in X]; y = [v/m for v in Y]
    P = [Ap, Am, 0]; Pb = [Ap, -Am, 0]
    ph = [v/m for v in P]; pbh = [v/m for v in Pb]
    fL = frame(pbh, ph); fR = frame(ph, pbh)
    ax = mp.sqrt(-mdot(x, x)); ay = mp.sqrt(-mdot(y, y))
    lpw = mp.log(2) + 6*mp.log(4) + 6*(mp.log(ax) + mp.log(ay)) + 2*mdot(P, x) + 2*mdot(Pb, y) - mp.log(4) - log_b(D, Jv, 0)
    res = {}
    for h in [(a, b, c, d) for a in (1, -1) for b in (1, -1) for c in (1, -1) for d in (1, -1)]:
        # degree 2 in each polarisation: pick the coefficient of (z.eps_h)^2 by z = conj(eps_h) = eps_{-h}  (eps_h.eps_{-h} = 1, eps_h^2 = 0)
        z1 = eps(*fL, -h[0]); z3 = eps(*fL, -h[2]); z2p = eps(*fR, -h[1]); z4p = eps(*fR, -h[3])
        ls, v = integrand(cL, cRt, D, Jv, x, y, z1, z2p, z3, z4p)
        res[h] = v*mp.exp(ls - lpw)
    # non-helicity content: Fourier modes of z1 on the null circle pbarhat + cos th e_par + sin th e_perp (other polarisations fixed):
    # the leading tensor has only modes q1 = +-2; leakage = max(|q1 = 0, +-1|)/max(|q1 = +-2|)
    th = [2*mp.pi*k/5 for k in range(5)]
    vals = []
    for t in th:
        z1 = [pbh[i] + mp.cos(t)*fL[0][i] + mp.sin(t)*fL[1][i] for i in range(3)]
        ls, v = integrand(cL, cRt, D, Jv, x, y, z1, eps(*fR, 1), eps(*fL, 1), eps(*fR, -1))
        vals.append(v*mp.exp(ls - lpw))
    modes = {q: sum(vals[k]*mp.exp(-1j*q*th[k]) for k in range(5))/5 for q in (-2, -1, 0, 1, 2)}
    leak = max(abs(modes[q]) for q in (-1, 0, 1))/max(abs(modes[2]), abs(modes[-2]))
    return res, leak, m

if __name__ == "__main__":
    pr("# <TTTT> heavy shadow integrand with Hilbert-90 Weyl-covariant families (fourpt/t7_vertex_l2.py)")
    pr("# fixed beta exactly: A_+ = s A_+0, A_- = s A_-0 (J = A_- - 1/2 non-integer is allowed for the 3-pt functions), s = 1, 2, 4, 8")
    X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
    Ap0, Am0 = mp.mpf('48.75'), mp.mpf('20.5')
    HS = [(a, b, c, d) for a in (1, -1) for b in (1, -1) for c in (1, -1) for d in (1, -1)]
    data = []
    import pickle
    for sc_ in (1, 2, 4, 8):
        D = mp.mpf(3)/2 + sc_*Ap0; Jv = sc_*Am0 - mp.mpf(1)/2
        fn = os.path.join(HERE, 'cache', 't7_scale%d.pkl' % sc_)
        cached_ = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else None
        fams = None if cached_ else [weyl_family(D, Jv, r) for r in RULES]
        R = {}
        for a in range(2):
            for b in range(2):
                if cached_:
                    res, leak, m = cached_[(a, b)]
                else:
                    res, leak, m = helicity(fams[a][0], fams[b][1], D, Jv, X, Y)
                    pickle.dump({**(pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}), (a, b): (res, leak, m)}, open(fn, 'wb'))
                R[(a, b)] = res
                big = max(abs(v) for v in res.values())
                pr(f"  m={mp.nstr(m, 6)} (a,b)=({a+1},{b+1}): non-helicity leakage {mp.nstr(leak, 3)};  helicity components/max: " +
                   " ".join(f"{''.join('+' if s_ > 0 else '-' for s_ in h)}:{mp.nstr((v/big).real, 5)}" for h, v in res.items()))
        dev = 0
        for h in HS:
            lhs = R[(0, 1)][h]*R[(1, 0)][(h[2], h[3], h[0], h[1])]
            rhs = R[(0, 0)][(h[0], h[1], h[0], h[1])]*R[(1, 1)][(h[2], h[3], h[2], h[3])]
            dev = max(dev, abs(lhs/rhs - 1))
        par = max(abs(R[(a, b)][h]/mp.conj(R[(a, b)][tuple(-x for x in h)]) - 1) for a in range(2) for b in range(2) for h in HS)
        pr(f"  m={mp.nstr(m, 6)}: factorisation T^12_(h12,h34) T^21_(h34,h12)/(T^11_(h12,h12) T^22_(h34,h34)) - 1, max_h: {mp.nstr(dev, 3)};  parity |T_h/conj(T_-h) - 1| max: {mp.nstr(par, 3)};  null-space gap sv_3/sv_4 (max) = {max(GAPS):.1e}")
        data.append((m, R))
    # Richardson test of the 1/m structure, no closed form needed: ratios to a reference component
    pr("# Richardson ratios (g(2s)-g(s))/(g(4s)-g(2s)) for s = 1, 2:  -> 2 for an O(1/m) correction, -> 4 for O(1/m^2)")
    href = (1, -1, -1, 1)
    for (a, b) in ((0, 1), (0, 0)):
        for h in [(1, 1, -1, -1), (1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, 1)]:
            ht = (h[2], h[3], h[0], h[1]); hrt = (href[2], href[3], href[0], href[1])
            g = [R[(a, b)][h]/R[(a, b)][href] for (m, R) in data]
            f = [(R[(a, b)][h] + R[(b, a)][ht])/(R[(a, b)][href] + R[(b, a)][hrt]) for (m, R) in data]
            rg = [(g[i+1] - g[i])/(g[i+2] - g[i+1]) for i in range(2)]
            rf = [(f[i+1] - f[i])/(f[i+2] - f[i+1]) for i in range(2)]
            pr(f"  (a,b)=({a+1},{b+1}) h={h}: unsymmetrised g = {[mp.nstr(x.real, 8) for x in g]}, ratios {[mp.nstr(x.real, 4) for x in rg]};"
               f"  Weyl-symmetrised f = {[mp.nstr(x.real, 8) for x in f]}, ratios {[mp.nstr(x.real, 4) for x in rf]}")
    open(os.path.join(HERE, 'out_t7_vertex_l2.txt'), 'w').write("\n".join(out) + "\n")
