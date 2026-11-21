"""Result Q, Stage 5 (RESULT_Q_STAGE5.md §1): EXACT leading-Regge heavy image of the raw SPINNING k = 5 functional Psi_{5,nu}[o_H q]
(j' = 1, f = f-_{5;2}, the C^_3 (1,3)/(2,4) structures), by the stage4_hj / v4_hclosed method extended to RANK-2 equivariant tensors with
the parity-odd (eps-carrying) j' = 1 sides of q1_c3_regge_kernel.

Kernel.  Lemma B (Stage 3): f-_{5;2} = f-_{3;2} R, R -> 16 eta^2 r^2 at leading Regge.  The k = 3 kernel (q1) is (|x||y|)^{a-2} K0, per side
|x|^{b-3} with b = 8 and no gradient; the k = 5 kernel is 16 r^2 eta^2 x (k = 3 kernel) = 16 (x.y)^2 x (k = 3 kernel), i.e. per side
|x|^{b-3} x_mu x_nu with the SAME b = 8.  Hence
    L_5(p, pb) = +16 < d^2 Ahat(p, z0), d^2 Bhat(pb, z0) >     (full contraction with eta^{-1}; z0-integral over the celestial circle)
relative to L_3 = Ahat Bhat (q1's normalisation; the common sign of f- is dropped in both, so the relative k = 3 / k = 5 sign is kept).
Side monomials: q1's A3_monomials('C3') (o_H q, with the z0-dependent w-vectors w_e = eta.cross(z0, e)) and B_monomials_j(1) (j' = 1
celestial side), with sympy coefficients (lambdify intercepted, as v4_hclosed._monos_sym).  Each monomial coef |x|^c (-z0.x)^-lam
(e.x)^i (f.x)^j (w.x)^k transforms to K(c,lam) (e.d)^i (f.d)^j (w.d)^k F_{c,lam}; the extra x_mu x_nu is d_mu d_nu (lower-index tensor).
Derivative patterns: truncated Taylor expansion of S^-lam R^kap in 3 directional slots (e, f, w) + 2 coordinate slots.
Rank-2 decomposition at the reference point (p = (1,0,0), z0 = (1,0,1), zh = z0/(-z0.p)): basis {pp, pz, zz, eta} (n-even) and
{pn, zn} (n-odd, n = s eps(p, zh, .)).  For the parity-odd o_H q / j' = 1 sides the n-odd parts need NOT vanish (unlike stage 4); they
are kept.  Contraction in the Gram basis {p, pb, zh_A}: n_A = G^-1 (0, sV, 0), n_B = G^-1 (-s r V, 0, 0), V = eps(p, zh_A, pb),
V^2 = -det G; terms odd in V integrate to zero over the z0 circle (reflection fixing p, pb); V^2 terms are polynomial.
Then sum C^k r^s -> sum_j R_j P_{l+j} (stage4_hj.to_harmonics).
Outputs R_j(nu), j in {0, +-2} (and any +-1, asserted absent), per helicity class: (h1, h2) = (+,-) [++] and (+,+) [+-], and the
dual-line / light-ray-line FLIP configurations.  Checks: rank 0 vs q1's numerical c3 (exact ratio = known nu-only factor); rank-2 scalar
(j' = 0, no spin) vs the stage-4 proxy a^(2)_j(l) E_10(nu - ij); n-odd parts.
Output: cache/stage5_hj.txt, cache/stage5_hj.pkl.
"""
import os, sys, io, time, pickle, itertools, functools
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
H7 = os.path.join(HERE, '..', 'h7i14'); SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, HERE); sys.path.insert(0, H7); sys.path.insert(0, SF)
import numpy as np
import sympy as sp
import v3_regge_kernel as V
import v4_hclosed as H4
import q1_c3_regge_kernel as Q1
import stage4_hj as S4

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
OUTP = os.path.join(CACHE, 'stage5_hj.pkl')
d = V._dl
nu_s = H4.nu_s
I = sp.I
lam_s, kap_s = sp.symbols('Lam Kap')
ETA = sp.diag(-1, 1, 1)
md = S4.md
BCODE = 8                                   # radial weight of the k = 3 (and k = 5) per-side transform (q1)


# ------------------------------------------------------------------ symbolic monomials (lambdify intercepted)
def _sym(fun, *args):
    orig = Q1.sp.lambdify
    Q1.sp.lambdify = lambda a, expr, *x, **k: expr
    try:
        return fun.__wrapped__(*args)
    finally:
        Q1.sp.lambdify = orig


@functools.lru_cache(maxsize=None)
def monos(kind):
    """list of (coef, k, m, i, j, wkey).  kind: 'A3' (o_H q), 'B1' (j' = 1 light-ray side), 'Asc', 'Bsc' (scalar, j' = 0)."""
    if kind == 'A3':
        return _sym(Q1.A3_monomials, 'C3')
    if kind == 'B1':
        return _sym(Q1.B_monomials_j, 1, True)
    if kind == 'Asc':
        return [(c, k, m, i, j, None) for (c, k, m, i, j) in H4._monos_sym(False, 'A')]
    if kind == 'Bsc':
        return [(c, k, m, i, j, None) for (c, k, m, i, j) in H4._monos_sym(False, 'B')]
    raise KeyError(kind)


# ------------------------------------------------------------------ Taylor patterns with 3 directional slots
def _pmul(A, B, caps, rank):
    out = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            if any(k[i] > caps[i] for i in range(3)) or sum(k[3:]) > rank:
                continue
            out[k] = out.get(k, 0) + va*vb
    return out


def taylor_patterns3(e, f, w, z, rank):
    """dict (i, j, k, s0, s1, s2) -> (e.d)^i (f.d)^j (w.d)^k d0^s0 d1^s1 d2^s2 [S^-lam R^kap]/[S^-lam R^kap] at p0 = (1,0,0),
    i, j <= 2, k <= 1, s0 + s1 + s2 = rank (lower-index derivatives w.r.t. p^mu)."""
    caps = (2, 2, 1)
    vecs = [e, f, w, (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    def lin(v):
        P = {}
        for a in range(6):
            key = [0]*6; key[a] = 1
            c = sp.nsimplify(md(v, vecs[a])) if False else md(v, vecs[a])
            if c != 0:
                P[tuple(key)] = c
        return P
    p0 = (1, 0, 0)
    sig = {k: -v for k, v in lin(z).items()}
    hp = lin(p0)
    hh = {}
    for a in range(6):
        for b in range(6):
            c = md(vecs[a], vecs[b])
            if c != 0:
                key = [0]*6; key[a] += 1; key[b] += 1
                hh[tuple(key)] = hh.get(tuple(key), 0) + c
    rho = {k: -2*v for k, v in hp.items()}
    for k, v in hh.items():
        rho[k] = rho.get(k, 0) - v
    Nmax = 5 + rank
    def binser(X, alpha):
        out = {(0,)*6: sp.Integer(1)}; term = {(0,)*6: sp.Integer(1)}; coef = sp.Integer(1)
        for n in range(1, Nmax + 1):
            coef = coef*(alpha - n + 1)/n
            term = _pmul(term, X, caps, rank)
            if not term:
                break
            for k, v in term.items():
                out[k] = out.get(k, 0) + coef*v
        return out
    F = _pmul(binser(sig, -lam_s), binser(rho, kap_s), caps, rank)
    pats = {}
    for k, v in F.items():
        if sum(k[3:]) != rank:
            continue
        fac = 1
        for x in k:
            fac *= sp.factorial(x)
        pats[k] = sp.expand(v*fac)
    return pats


# reference (X-exponent k, s-exponent m) per monomial family: all K(c, lam) of a family differ from K(c_ref, lam_ref) by integer shifts
# of both Gamma arguments (asserted).  N(nu) := 2 pi K(c_A^ref, lam_A^ref) K(c_B^ref, lam_B^ref), c = (b-3) + k, lam = -m.
REF = {'A3': (4 - d, d - 5), 'B1': (d + 1, -d - 3), 'Asc': (1 - d, -(1 - d)), 'Bsc': (d, -d)}


def Nref(nu, kindA='A3', kindB='B1', b=BCODE):
    """numerical common factor 2 pi K_A^ref K_B^ref (v3 logKc), complex."""
    dl = 0.5 + 1j*np.asarray(nu, complex)
    out = 2*np.pi
    for kind in (kindA, kindB):
        k, m = REF[kind]
        c = (b - 3) + sp.lambdify(d, k, 'numpy')(dl); lam = -sp.lambdify(d, m, 'numpy')(dl)
        out = out*np.exp(V.logKc(c, lam))
    return out


_s2 = 1/sp.sqrt(2)
EPS = {+1: (0, _s2, I*_s2), -1: (0, _s2, -I*_s2)}     # eps_h at p = (1,0,0) (frames(0), parallel frame, positive orientation)


def wvec_sym(z, e):
    c = (z[1]*e[2] - z[2]*e[1], z[2]*e[0] - z[0]*e[2], z[0]*e[1] - z[1]*e[0])
    return (-c[0], c[1], c[2])


def side_tensor(kind, side, rank, ha, hb, b=BCODE):
    """symmetric rank-`rank` tensor d^rank (side transform)/K^ref at the reference celestial point, dict (s0,s1,s2) -> expr(delta).
    The pair's polarisations: e = eps_{ha}, f = eps_{hb} (preserving: hb = -ha; flip: hb = ha)."""
    z = (1, 0, 1)
    e, f = EPS[ha], EPS[hb]
    kref, mref = REF[kind]
    if side == 'A':
        wv = {'WE': wvec_sym(z, e), 'WF': wvec_sym(z, f), None: (0, 0, 0)}
    else:
        wv = {'W1': wvec_sym(z, e), 'W3': wvec_sym(z, f), None: (0, 0, 0)}
    cargs = (md(z, e), md(z, f), md(e, f))
    cref = (b - 3) + kref; lref = -mref
    pats = {wk: taylor_patterns3(e, f, w, z, rank) for wk, w in wv.items()}
    T = {}
    for (coef, kf, mf, i, j, wk) in monos(kind):
        k = sp.sympify(kf); m = sp.sympify(mf)
        c = (b - 3) + k; lam = -m; kap = (2*lam - c - 3)/2
        if side == 'A':
            cv = sp.sympify(coef).subs({V._a2: cargs[0], V._a4: cargs[1], V._e24: cargs[2]})
        else:
            cv = sp.sympify(coef).subs({V._A1: cargs[0], V._A3: cargs[1], V._E13: cargs[2]})
        n1 = sp.simplify(c/2 - cref/2); n2 = sp.simplify((c + 3)/2 - lam - (cref + 3)/2 + lref)
        dmu = sp.simplify((c + 3 - lam) - (cref + 3 - lref))
        assert n1.is_integer and n2.is_integer and dmu.is_integer, (n1, n2, dmu)
        Kr = sp.Integer(2)**dmu*H4._gratio(1 + cref/2, n1)*H4._gratio((cref + 3)/2 - lref, n2)
        kw = 0 if wk is None else 1
        for key, pv in pats[wk].items():
            if key[0] != i or key[1] != j or key[2] != kw:
                continue
            ms = key[3:]
            T[ms] = T.get(ms, 0) + cv*Kr*pv.subs({lam_s: lam, kap_s: kap})
    return {k: sp.factor(sp.simplify(v)) for k, v in T.items()}


def basis2():
    P = (-1, 0, 0); Z = (-1, 0, 1); N = (0, 1, 0)
    return [('pp', S4._outer([P, P])), ('pz', S4._outer([P, Z])), ('zz', S4._outer([Z, Z])), ('ee', S4._sym_eta([], ETA)),
            ('pn', S4._outer([P, N])), ('zn', S4._outer([Z, N]))]


def decompose(T, rank):
    if rank == 0:
        return {'1': sp.factor(sp.simplify(T.get((0, 0, 0), 0)))}
    Bs = basis2()
    comps = sorted(set(tuple(sorted(idx)) for idx in itertools.product(range(3), repeat=2)))
    def cnt(idx):
        c = [0, 0, 0]
        for a in idx:
            c[a] += 1
        return tuple(c)
    M = sp.Matrix([[B.get(idx, 0) for _, B in Bs] for idx in comps])
    bb = sp.Matrix([T.get(cnt(idx), 0) for idx in comps])
    sol = M.solve(bb)
    return {nm: sp.factor(sp.simplify(sol[k])) for k, (nm, _) in enumerate(Bs)}


C_s, r_s, V_s = sp.symbols('C r Vol')


def contract2(cA, cB, rank):
    """< T_A(p, zA), T_B(pb, zB) > in the Gram basis {p, pb, zA}; zB = r zA; n-odd vectors via V = eps(p, zA, pb). Returns (even part
    as a polynomial in C, r, 1/r after V^2 -> -det G, coefficient of V^1 (must integrate to 0))."""
    G = sp.Matrix([[-1, -C_s, -1], [-C_s, -1, -1/r_s], [-1, -1/r_s, 0]])
    if rank == 0:
        return sp.expand(cA['1']*cB['1']), 0
    Gi = sp.simplify(G.inv())
    p = (1, 0, 0); pb = (0, 1, 0); zA = (0, 0, 1); zB = (0, 0, r_s)
    nA = tuple(sp.simplify(x) for x in Gi*sp.Matrix([0, V_s, 0]))
    nB = tuple(sp.simplify(x) for x in Gi*sp.Matrix([-r_s*V_s, 0, 0]))
    def build(coef, P, Z, N):
        Bs = {'pp': S4._outer([P, P]), 'pz': S4._outer([P, Z]), 'zz': S4._outer([Z, Z]), 'ee': S4._sym_eta([], Gi),
              'pn': S4._outer([P, N]), 'zn': S4._outer([Z, N])}
        T = {}
        for nm, B in Bs.items():
            cf = coef.get(nm, 0)
            if cf == 0:
                continue
            for k, v in B.items():
                T[k] = T.get(k, 0) + cf*v
        return T
    TA = build(cA, p, zA, nA); TB = build(cB, pb, zB, nB)
    tot = 0
    for ia, va in TA.items():
        for ib, vb in TB.items():
            g = G[ia[0], ib[0]]*G[ia[1], ib[1]]
            if g != 0:
                tot += va*vb*g
    tot = sp.expand(sp.cancel(sp.together(tot)))
    detG = sp.simplify(G.det())
    P = sp.Poly(tot, V_s)
    even = 0; odd = 0
    for (k,), c in P.terms():
        if k % 2 == 0:
            even += c*(-detG)**(k//2)
        else:
            odd += c*V_s**k
    return sp.expand(sp.cancel(sp.together(even))), sp.simplify(odd)


def Rj(kindA, kindB, rank, h1, h2, flipA=False, flipB=False, b=BCODE, log=None):
    """R_j(nu) of the heavy image sum_j N R_j P_{l+j}; (h1, h2) = helicities on the (1,3) and (2,4) lines."""
    TA = side_tensor(kindA, 'A', rank, h2, h2 if flipA else -h2, b)
    TB = side_tensor(kindB, 'B', rank, h1, h1 if flipB else -h1, b)
    cA = decompose(TA, rank); cB = decompose(TB, rank)
    if log:
        log(f"      side A coefs: { {k: str(v) for k, v in cA.items() if v != 0} }")
        log(f"      side B coefs: { {k: str(v) for k, v in cB.items() if v != 0} }")
    pref = 16 if rank == 2 else 1
    even, odd = contract2(cA, cB, rank)
    H = S4.to_harmonics(pref*even) if even != 0 else {}
    sub = {d: sp.Rational(1, 2) + I*nu_s}
    return {j: sp.factor(sp.simplify(v.subs(sub))) for j, v in H.items()}, cA, cB, odd


def build():
    out = {}
    t0 = time.time()
    for rank in (0, 2):
        for cls, (h1, h2) in (('pp', (+1, -1)), ('pm', (+1, +1)), ('mm', (-1, +1)), ('mp', (-1, -1))):
            R, cA, cB, odd = Rj('A3', 'B1', rank, h1, h2)
            out[f'C{3 + rank}_{cls}'] = R
            out[f'C{3 + rank}_{cls}_sides'] = ({k: v for k, v in cA.items()}, {k: v for k, v in cB.items()})
            out[f'C{3 + rank}_{cls}_odd'] = odd
            print(f"  rank {rank} {cls} done [{time.time() - t0:.0f}s]", flush=True)
    # flips (dual line: z2' = z4'; light-ray line: z1 = z3), rank 0 and 2
    for rank in (0, 2):
        for nm, fa, fb in (('flipA', True, False), ('flipB', False, True), ('flipAB', True, True)):
            R, cA, cB, odd = Rj('A3', 'B1', rank, +1, -1, flipA=fa, flipB=fb)
            out[f'C{3 + rank}_{nm}'] = R
            out[f'C{3 + rank}_{nm}_sides'] = (cA, cB)
            out[f'C{3 + rank}_{nm}_odd'] = odd
        print(f"  flips rank {rank} done [{time.time() - t0:.0f}s]", flush=True)
    # scalar control (j' = 0, no spin), rank 0 and 2
    for rank in (0, 2):
        R, cA, cB, odd = Rj('Asc', 'Bsc', rank, +1, -1)
        out[f'sc{rank}'] = R
        out[f'sc{rank}_sides'] = (cA, cB)
    print(f"  scalar done [{time.time() - t0:.0f}s]", flush=True)
    ser = {}
    for k, D in out.items():
        if isinstance(D, dict):
            ser[k] = ('d', {j: sp.srepr(v) for j, v in D.items()})
        elif isinstance(D, tuple):
            ser[k] = ('t', tuple({j: sp.srepr(v) for j, v in x.items()} for x in D))
        else:
            ser[k] = ('e', sp.srepr(D))
    pickle.dump(ser, open(OUTP, 'wb'))
    return out


def load():
    if not os.path.exists(OUTP):
        return build()
    ser = pickle.load(open(OUTP, 'rb'))
    ev = lambda s: eval(s, {**sp.__dict__})
    out = {}
    for k, (t, D) in ser.items():
        if t == 'd':
            out[k] = {j: ev(v) for j, v in D.items()}
        elif t == 't':
            out[k] = tuple({j: ev(v) for j, v in x.items()} for x in D)
        else:
            out[k] = ev(D)
    return out


def numfuncs(key='C5_pp'):
    return {j: sp.lambdify(nu_s, v, 'numpy') for j, v in load()[key].items()}


def a2_coef(j, l):
    c = {0: 1}
    for _ in range(2):
        n = {}
        for jj, v in c.items():
            m = l + jj
            n[jj + 1] = n.get(jj + 1, 0) + v*(m + 1)/(2*m + 1)
            n[jj - 1] = n.get(jj - 1, 0) + v*m/(2*m + 1)
        c = n
    return c.get(j, 0)


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t00 = time.time()
    D = build()
    pr(f"# stage5_hj: exact R_j of the leading-Regge heavy image of Psi_{{5,nu}}[o_H q] (j' = 1): image = N sum_j R_j P_{{l+j}}(cosh beta),"
       f" N = 2 pi K_A^ref K_B^ref x (dropped j'=1 celestial moment <|s|^(1-delta)>), b = {BCODE}  [{time.time() - t00:.0f}s]")
    for k in sorted(D):
        if k.endswith('_odd'):
            pr(f"  {k}: V-odd part of the contraction (integrates to 0): {D[k]}")
    # sides n-odd content
    for k in ('C5_pp_sides', 'C5_pm_sides'):
        cA, cB = D[k]
        pr(f"  {k}: side A coefs { {n: str(v) for n, v in cA.items() if v != 0} }")
        pr(f"  {k}: side B coefs { {n: str(v) for n, v in cB.items() if v != 0} }")
    pr("(K3) rank 0 (k = 3): shifts present and R_0 per class:")
    for cls in ('pp', 'pm', 'mm', 'mp'):
        pr(f"   C3 {cls}: {D['C3_' + cls]}")
    pr("(K5) rank 2 (k = 5): shifts present and R_j per class:")
    for cls in ('pp', 'pm', 'mm', 'mp'):
        pr(f"   C5 {cls}: shifts {sorted(D['C5_' + cls])}")
        for j in sorted(D['C5_' + cls]):
            pr(f"      R_{j:+d} = {D['C5_' + cls][j]}")
    pr("(H) helicity: R_j(+-) + R_j(++) (sigma_h <=> 0), R_j(+-) - R_j(++) (blind <=> 0), total-flip (mm vs pp, mp vs pm):")
    for rk in (3, 5):
        pp, pm, mm, mp = (D[f'C{rk}_{c}'] for c in ('pp', 'pm', 'mm', 'mp'))
        js = sorted(set(pp) | set(pm))
        pr(f"   k = {rk}: sum   " + ", ".join(f"j={j}: {sp.simplify(pm.get(j, 0) + pp.get(j, 0))}" for j in js))
        pr(f"   k = {rk}: diff  " + ", ".join(f"j={j}: {sp.simplify(pm.get(j, 0) - pp.get(j, 0))}" for j in js))
        pr(f"   k = {rk}: mm-pp " + ", ".join(f"j={j}: {sp.simplify(mm.get(j, 0) - pp.get(j, 0))}" for j in js) +
           ";  mp-pm " + ", ".join(f"j={j}: {sp.simplify(mp.get(j, 0) - pm.get(j, 0))}" for j in js))
    pr("(F) flips (z2' = z4': flipA; z1 = z3: flipB): R_j")
    for rk in (3, 5):
        for nm in ('flipA', 'flipB', 'flipAB'):
            pr(f"   k = {rk} {nm}: {D[f'C{rk}_{nm}']};  sides A nonzero: {[n for n, v in D[f'C{rk}_{nm}_sides'][0].items() if v != 0]}, "
               f"B nonzero: {[n for n, v in D[f'C{rk}_{nm}_sides'][1].items() if v != 0]}")
    pr("(R) ratio R_j^{k=5}(++)/R_0^{k=3}(++):")
    for j in sorted(D['C5_pp']):
        pr(f"   j = {j:+d}: {sp.factor(sp.simplify(D['C5_pp'][j]/D['C3_pp'][0]))}")
    pr("    symmetry R_{-j}(nu) - R_j(-nu): " + ", ".join(f"j={j}: {sp.simplify(D['C5_pp'].get(-j, 0) - D['C5_pp'].get(j, 0).subs(nu_s, -nu_s))}"
                                                    for j in (0, 2)))
    # scalar control vs proxy
    import mpmath as mp
    E10 = lambda v: mp.gamma((10 - mp.mpf(1)/2 - 1j*v)/2)**2*mp.gamma((10 - mp.mpf(1)/2 + 1j*v)/2)**2
    pr(f"(C2) rank-2 scalar (j' = 0, no spin, b = {BCODE}) R_j vs the stage-4 proxy a^(2)_j(l) E_10(nu - ij) (ratio to j = 0, minus 1); "
       f"shifts {sorted(D['sc2'])}; rank-0 scalar shifts {sorted(D['sc0'])}")
    sc = {j: sp.lambdify(nu_s, v, 'mpmath') for j, v in D['sc2'].items()}
    for nv in (0.3, 1.7, 5.2):
        l = -mp.mpf(1)/2 + 1j*nv
        rr = []
        for j in sorted(sc):
            ex = sc[j](mp.mpf(nv))/sc[0](mp.mpf(nv))
            px = a2_coef(j, l)*E10(nv - 1j*j)/(a2_coef(0, l)*E10(nv))
            rr.append(f"j={j}: {complex(ex/px - 1):.1e}")
        pr(f"   nu = {nv}: " + ", ".join(rr))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'stage5_hj.txt'), 'w').write(buf.getvalue())
