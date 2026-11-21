"""Result Q, Stage 4 (RESULT_Q_STAGE4.md §1): EXACT leading-Regge heavy-image coefficients h_j(nu), j = +-1, +-3, of the raw subtracted
k = 4 CLS-D functional Psi_{4,nu}[EE] (and of the scalar skeleton), by the v4_hclosed / u2 method extended to RANK-3 equivariant tensors.

Kernel.  Stage 12 §1.2 / Stage 3 Lemma B: f+_{4;2} = f+_{2;2} R, R -> 16 eta^2 r^2 at leading Regge.  The k = 2 eta-kernel is
eta K0 = -(x.y)/(|x||y|) K0 (per side |x|^{b-3} x^mu, b = a - 1 = 6, a = 7); the k = 4 kernel is
16 r^2 eta^3 K0 = -16 (x.y)^3/(|x||y|) K0, i.e. per side |x|^{b-3} x^mu x^nu x^rho with the SAME b = 6.  Hence
    L_4(p, pb) = -16 < d^3 Ahat(p, z0), d^3 Bhat(pb, z0) >      (full contraction with eta^{-1}; z0-integral over the celestial circle)
with Ahat, Bhat the plane-wave transforms of v3_regge_kernel (same monomials, same K(c, lam) constants as v4_hclosed).
Rank-3 equivariance: at one celestial point (p = (1,0,0), z0 = (1,0,1)) d^3 Ahat is a symmetric covariant 3-tensor; it is decomposed on
{p^3, p^2 zh, p zh^2, zh^3, eta p, eta zh} (+ the 4 n-odd tensors, checked to vanish), zh = z0/(-z0.p); covariance + homogeneity in z0 fix
it everywhere.  Contraction with the other side (Gram algebra in the basis p, pb, zh_A; zh_B = r zh_A, r = (-z0.p)/(-z0.pb)) gives
sum C^k r^s alpha^-lam gamma^-lam', whose z0-integral is 2 pi C^k P_{l+s}(C) (I(a) of v4_hclosed), C = cosh beta; C^k P by the Legendre
recurrence.  Result: L_4 = 2 pi K_A^ref K_B^ref sum_j R_j(nu) P_{l+j}(cosh beta), R_j rational in nu (exact).
Rank 1 with the same code reproduces v4_hclosed (q^2) and u2 (EE) exactly (check printed).  Scalar: the ratios R_j/R_j' reproduce
stage12_rawkernel's a^(3)_j(l) E_9(nu - i j) (check printed).
Output: cache/stage4_hj.txt, cache/stage4_hj.pkl (srepr of R_j for 'sc', 'EE', 'q2' (rank 3) and rank-1 checks).  ~ minutes.
"""
import os, sys, io, time, pickle, itertools
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
H7 = os.path.join(HERE, '..', 'h7i14'); SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, H7); sys.path.insert(0, SF)
import numpy as np
import sympy as sp
import v3_regge_kernel as V
import v4_hclosed as H4

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
OUTP = os.path.join(CACHE, 'stage4_hj.pkl')
d = V._dl
nu_s = H4.nu_s
I = sp.I
lam_s, kap_s = sp.symbols('Lam Kap')
ETA = sp.diag(-1, 1, 1)


def md(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


# ------------------------------------------------------------------ derivative patterns by truncated Taylor expansion
def _pmul(A, B, caps):
    out = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            if any(k[i] > caps[i] for i in range(2)) or sum(k[2:]) > caps[2]:
                continue
            out[k] = out.get(k, 0) + va*vb
    return out


def taylor_patterns(e, f, z, rank):
    """dict (i, j, s0, s1, s2) -> (e.d)^i (f.d)^j d0^s0 d1^s1 d2^s2 [S^-lam R^kap] / [S^-lam R^kap] at p0 = (1,0,0)
    (derivatives w.r.t. the contravariant components p^mu, i.e. lower-index tensors), for i, j <= 2, s0+s1+s2 = rank."""
    caps = (2, 2, rank)
    # h = a e + b f + sum s_mu E_mu ; monomial exponent key (ia, ib, s0, s1, s2)
    def lin(vec):          # vec . h as a polynomial (dict)
        P = {}
        P[(1, 0, 0, 0, 0)] = md(vec, e)
        P[(0, 1, 0, 0, 0)] = md(vec, f)
        for mu in range(3):
            key = [0, 0, 0, 0, 0]; key[2 + mu] = 1
            E = [0, 0, 0]; E[mu] = 1
            P[tuple(key)] = md(vec, E)
        return {k: v for k, v in P.items() if v != 0}
    p0 = (1, 0, 0)
    sig = {k: -v for k, v in lin(z).items()}                          # S = 1 + sig
    hp = lin(p0)
    # h.h
    hh = {}
    basis = {(1, 0, 0, 0, 0): e, (0, 1, 0, 0, 0): f, (0, 0, 1, 0, 0): (1, 0, 0), (0, 0, 0, 1, 0): (0, 1, 0), (0, 0, 0, 0, 1): (0, 0, 1)}
    for (k1, v1), (k2, v2) in itertools.product(basis.items(), basis.items()):
        k = tuple(x + y for x, y in zip(k1, k2)); c = md(v1, v2)
        if c != 0:
            hh[k] = hh.get(k, 0) + c
    rho = {k: -2*v for k, v in hp.items()}
    for k, v in hh.items():
        rho[k] = rho.get(k, 0) - v
    Nmax = 4 + rank
    def binser(X, alpha):
        out = {(0, 0, 0, 0, 0): sp.Integer(1)}; term = {(0, 0, 0, 0, 0): sp.Integer(1)}; coef = sp.Integer(1)
        for n in range(1, Nmax + 1):
            coef = coef*(alpha - n + 1)/n
            term = _pmul(term, X, caps)
            for k, v in term.items():
                out[k] = out.get(k, 0) + coef*v
        return out
    F = _pmul(binser(sig, -lam_s), binser(rho, kap_s), caps)
    pats = {}
    for k, v in F.items():
        if sum(k[2:]) != rank:
            continue
        fac = sp.factorial(k[0])*sp.factorial(k[1])*sp.factorial(k[2])*sp.factorial(k[3])*sp.factorial(k[4])
        pats[k] = sp.expand(v*fac)
    return pats


# ------------------------------------------------------------------ side tensors
def monos_scalar(side):
    return H4._monos_sym(False, side)


def monos_spinB():
    return H4._monos_sym(True, 'B')


def monos_struct(name):
    import u2_hclosed_struct as U2
    return U2.monos_struct_sym(name)


def side_tensor(monos, side, rank, hel=+1):
    """symmetric rank-`rank` tensor d^rank (side transform)/K^ref at the reference celestial point, as dict multi-index -> expr(delta).
    side 'A' (pair 2,4: e = eps_{h}, f = eps_{-h}, h = -hel) or 'B' (pair 1,3: e = eps_{hel}, f = eps_{-hel})."""
    z = (1, 0, 1)
    s2 = 1/sp.sqrt(2)
    em = (0, s2, -I*s2); ep = (0, s2, I*s2)
    if side == 'A':
        e, f = (em, ep) if hel == +1 else (ep, em)
        kref, mref = 1 - d, -(1 - d)
    else:
        e, f = (ep, em) if hel == +1 else (em, ep)
        kref, mref = d, -d
    cargs = (md(z, e), md(z, f), md(e, f))
    pats = taylor_patterns(e, f, z, rank)
    cref = 3 + kref; lref = -mref
    T = {}
    for (coef, kf, mf, i, j) in monos:
        k = sp.sympify(kf); m = sp.sympify(mf)
        c = 3 + k; lam = -m; kap = (2*lam - c - 3)/2
        if side == 'A':
            cv = sp.sympify(coef).subs({V._a2: cargs[0], V._a4: cargs[1], V._e24: cargs[2]})
        else:
            cv = sp.sympify(coef).subs({V._A1: cargs[0], V._A3: cargs[1], V._E13: cargs[2]})
        n1 = sp.simplify(c/2 - cref/2); n2 = sp.simplify((c + 3)/2 - lam - (cref + 3)/2 + lref)
        dmu = sp.simplify((c + 3 - lam) - (cref + 3 - lref))
        assert n1.is_integer and n2.is_integer and dmu.is_integer, (n1, n2, dmu)
        Kr = sp.Integer(2)**dmu*H4._gratio(1 + cref/2, n1)*H4._gratio((cref + 3)/2 - lref, n2)
        for key, pv in pats.items():
            if key[0] != i or key[1] != j:
                continue
            ms = key[2:]
            T[ms] = T.get(ms, 0) + cv*Kr*pv.subs({lam_s: lam, kap_s: kap})
    return {k: sp.factor(sp.simplify(v)) for k, v in T.items()}


# ------------------------------------------------------------------ decomposition on the covariant basis
def _outer(vs):
    """symmetrised outer product (sum over distinct orderings) of a list of 3-vectors -> rank-len(vs) sympy Array (as dict)."""
    r = len(vs); out = {}
    for perm in set(itertools.permutations(range(r))):
        for idx in itertools.product(range(3), repeat=r):
            val = sp.Integer(1)
            for a, b in zip(idx, perm):
                val *= vs[b][a]
            if val != 0:
                out[idx] = out.get(idx, 0) + val
    return out


def _sym_eta(vecs, G):
    """symmetrised eta (x) vecs (rank 2 + len(vecs)), with eta given as a symmetric matrix G (indices in the chosen basis)."""
    r = 2 + len(vecs); out = {}
    slots = list(itertools.combinations(range(r), 2))
    for (a, b) in slots:
        rest = [k for k in range(r) if k not in (a, b)]
        for perm in set(itertools.permutations(range(len(vecs)))):
            for idx in itertools.product(range(3), repeat=r):
                val = G[idx[a], idx[b]]
                for k, pk in zip(rest, perm):
                    val = val*vecs[pk][idx[k]]
                if val != 0:
                    out[idx] = out.get(idx, 0) + val
    return out


def basis_tensors(rank, P, Z, N, G):
    """the covariant basis at rank 1 or 3. Returns list of (name, dict)."""
    if rank == 1:
        return [('p', _outer([P])), ('z', _outer([Z])), ('n', _outer([N]))]
    B = [('ppp', _outer([P, P, P])), ('ppz', _outer([P, P, Z])), ('pzz', _outer([P, Z, Z])), ('zzz', _outer([Z, Z, Z])),
         ('ep', _sym_eta([P], G)), ('ez', _sym_eta([Z], G)),
         ('ppn', _outer([P, P, N])), ('pzn', _outer([P, Z, N])), ('zzn', _outer([Z, Z, N])), ('nnn', _outer([N, N, N]))]
    return B


def decompose(T, rank):
    """coefficients of T (dict multi-index (s0,s1,s2) counts -> value of the derivative) on the covariant basis at the reference point
    (lower-index components: p = (-1,0,0), z = (-1,0,1), n = (0,1,0), eta = diag(-1,1,1))."""
    P = (-1, 0, 0); Z = (-1, 0, 1); Nn = (0, 1, 0)
    Bs = basis_tensors(rank, P, Z, Nn, ETA)
    # T given by derivative counts: component T_{mu nu rho} with multiset of indices
    comps = sorted(set(tuple(sorted(idx)) for idx in itertools.product(range(3), repeat=rank)))
    def cnt(idx):
        c = [0, 0, 0]
        for a in idx:
            c[a] += 1
        return tuple(c)
    rows = []; rhs = []
    for idx in comps:
        rows.append([B.get(idx, 0) for _, B in Bs])
        rhs.append(T.get(cnt(idx), 0))
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    sol = M.solve(b) if M.shape[0] == M.shape[1] else M.solve_least_squares(b)
    return {nm: sp.factor(sp.simplify(sol[k])) for k, (nm, _) in enumerate(Bs)}


# ------------------------------------------------------------------ contraction and z0-integral
C_s, r_s = sp.symbols('C r')


def contract(coefA, coefB, rank):
    """sum_{alpha,beta} a_alpha b_beta <B_alpha(p, zA), B_beta(pb, zB)> in the Gram basis {p, pb, zA}; zB = r zA. Returns poly in C, r, 1/r.
    (n-odd parts must vanish; asserted by the caller.)"""
    G = sp.Matrix([[-1, -C_s, -1], [-C_s, -1, -1/r_s], [-1, -1/r_s, 0]])
    Gi = sp.simplify(G.inv())
    p = (1, 0, 0); pb = (0, 1, 0); zA = (0, 0, 1); zB = (0, 0, r_s)
    def build(coef, P, Z):
        if rank == 1:
            Bs = {'p': _outer([P]), 'z': _outer([Z])}
        else:
            Bs = {'ppp': _outer([P, P, P]), 'ppz': _outer([P, P, Z]), 'pzz': _outer([P, Z, Z]), 'zzz': _outer([Z, Z, Z]),
                  'ep': _sym_eta([P], Gi), 'ez': _sym_eta([Z], Gi)}
        T = {}
        for nm, B in Bs.items():
            cf = coef.get(nm, 0)
            if cf == 0:
                continue
            for k, v in B.items():
                T[k] = T.get(k, 0) + cf*v
        return T
    TA = build(coefA, p, zA); TB = build(coefB, pb, zB)
    tot = 0
    for ia, va in TA.items():
        for ib, vb in TB.items():
            g = sp.Integer(1)
            for a, b in zip(ia, ib):
                g = g*G[a, b]
                if g == 0:
                    break
            if g != 0:
                tot += va*vb*g
    return sp.expand(sp.simplify(tot))


def to_harmonics(expr, maxshift=4):
    """sum c_{k,s} C^k r^s  ->  sum_j R_j P_{l+j}  (z0-integral: alpha^{-lam+s} gamma^{-lam'-s} -> P_{l+s}; C P_mu recurrence).
    Returns dict j -> rational(delta)."""
    l = d - 1
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    dpow = sp.Poly(den, r_s).degree() if den.has(r_s) else 0
    dd = sp.simplify(den/r_s**dpow)
    assert not dd.has(r_s) and not dd.has(C_s), dd
    poly = sp.Poly(sp.expand(num), C_s, r_s)
    out = {}
    for (k, sr), c in poly.terms():
        s = sr - dpow
        c = c/dd
        cur = {s: sp.Integer(1)}
        for _ in range(k):
            new = {}
            for j, v in cur.items():
                mu = l + j
                new[j + 1] = new.get(j + 1, 0) + v*(mu + 1)/(2*mu + 1)
                new[j - 1] = new.get(j - 1, 0) + v*mu/(2*mu + 1)
            cur = new
        for j, v in cur.items():
            out[j] = out.get(j, 0) + c*v
    return {j: sp.factor(sp.simplify(v)) for j, v in out.items() if sp.simplify(v) != 0}


def Rj(monA, monB, rank, hel=+1, verbose=None):
    TA = side_tensor(monA, 'A', rank, hel); TB = side_tensor(monB, 'B', rank, hel)
    cA = decompose(TA, rank); cB = decompose(TB, rank)
    odd = [nm for nm in cA if 'n' in nm and sp.simplify(cA[nm]) != 0] + [nm for nm in cB if 'n' in nm and sp.simplify(cB[nm]) != 0]
    if verbose:
        verbose(f"      side A coefficients: { {k: str(v) for k, v in cA.items() if v != 0} }")
        verbose(f"      side B coefficients: { {k: str(v) for k, v in cB.items() if v != 0} }")
    assert not odd, ('n-odd components', odd)
    pref = -16 if rank == 3 else -1
    L = contract(cA, cB, rank)
    H = to_harmonics(pref*L)
    sub = {d: sp.Rational(1, 2) + I*nu_s}
    return {j: sp.factor(sp.simplify(v.subs(sub))) for j, v in H.items()}, cA, cB


def build():
    out = {}
    t0 = time.time()
    monAsc, monBsc = monos_scalar('A'), monos_scalar('B')
    monAEE = monos_struct('EE'); monAq2 = monos_struct('q2'); monB = monos_spinB()
    out['sc1'] = Rj(monAsc, monBsc, 1)[0]
    out['q21'] = Rj(monAq2, monB, 1)[0]
    out['EE1'] = Rj(monAEE, monB, 1)[0]
    print(f"  rank 1 done [{time.time() - t0:.0f}s]", flush=True)
    out['sc3'] = Rj(monAsc, monBsc, 3)[0]
    print(f"  rank 3 scalar done [{time.time() - t0:.0f}s]", flush=True)
    out['EE3'] = Rj(monAEE, monB, 3)[0]
    print(f"  rank 3 EE done [{time.time() - t0:.0f}s]", flush=True)
    out['EE3m'] = Rj(monAEE, monB, 3, hel=-1)[0]
    print(f"  rank 3 EE (other helicity) done [{time.time() - t0:.0f}s]", flush=True)
    out['q23'] = Rj(monAq2, monB, 3)[0]
    print(f"  rank 3 q2 done [{time.time() - t0:.0f}s]", flush=True)
    pickle.dump({k: {j: sp.srepr(v) for j, v in D.items()} for k, D in out.items()}, open(OUTP, 'wb'))
    return out


def load():
    if not os.path.exists(OUTP):
        return build()
    D = pickle.load(open(OUTP, 'rb'))
    return {k: {j: eval(v, {**sp.__dict__}) for j, v in E.items()} for k, E in D.items()}


def numfuncs(key='EE3'):
    """{j: callable(nu)} for the rational R_j of `key` (h_j = N_com R_j, N_com even, positive, as v4_hclosed)."""
    return {j: sp.lambdify(nu_s, v, 'numpy') for j, v in load()[key].items()}


def cfun(v):
    return 32*(4*v**2 + 1)*(4*v**2 + 9)/((4*v**2 + 25)*(4*v**2 + 81))


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t00 = time.time()
    D = build()
    pr(f"# stage4_hj: exact R_j (h_j = 2 pi K_A^ref K_B^ref R_j = N_com R_j) of the leading-Regge heavy image sum_j h_j P_{{l+j}}(cosh beta)"
       f"  [{time.time() - t00:.0f}s]")
    # ---- rank-1 checks
    ref = H4.load()
    import u2_hclosed_struct as U2
    refEE = U2.load()['EE']
    pr("(C1) rank 1 vs v4_hclosed / u2 (h_-  <-> j = -1, h_+ <-> j = +1):")
    pr(f"   scalar: j=-1 diff {sp.simplify(D['sc1'].get(-1, 0) - ref['RmS'])}, j=+1 diff {sp.simplify(D['sc1'].get(1, 0) - ref['RpS'])}, other j: {[j for j in D['sc1'] if abs(j) != 1]}")
    pr(f"   q2    : j=-1 diff {sp.simplify(D['q21'].get(-1, 0) - ref['RmT'])}, j=+1 diff {sp.simplify(D['q21'].get(1, 0) - ref['RpT'])}")
    pr(f"   EE    : j=-1 diff {sp.simplify(D['EE1'].get(-1, 0) - refEE[0])}, j=+1 diff {sp.simplify(D['EE1'].get(1, 0) - refEE[1])}")
    # ---- rank-3 scalar vs stage12_rawkernel proxy: R_j / R_1 = a_j E_9(nu - i j) / (a_1 E_9(nu - i))
    import mpmath as mp
    E9 = lambda v: mp.gamma((9 - mp.mpf(1)/2 - 1j*v)/2)**2*mp.gamma((9 - mp.mpf(1)/2 + 1j*v)/2)**2
    mu = sp.Symbol('mu')
    def a3(j, l):
        coeffs = {0: 1}
        for _ in range(3):
            new = {}
            for jj, c in coeffs.items():
                m = l + jj
                new[jj + 1] = new.get(jj + 1, 0) + c*(m + 1)/(2*m + 1)
                new[jj - 1] = new.get(jj - 1, 0) + c*m/(2*m + 1)
            coeffs = new
        return coeffs[j]
    pr("(C2) rank-3 scalar R_j vs the stage-12 proxy a^(3)_j(l) E_9(nu - i j) (ratio to j = -1, minus 1):")
    sc = {j: sp.lambdify(nu_s, v, 'mpmath') for j, v in D['sc3'].items()}
    pr(f"   scalar rank-3 shifts present: {sorted(D['sc3'])}")
    for nv in (0.3, 1.7, 5.2):
        l = -mp.mpf(1)/2 + 1j*nv
        rr = []
        for j in sorted(sc):
            ex = sc[j](mp.mpf(nv))/sc[-1](mp.mpf(nv))
            px = a3(j, l)*E9(nv - 1j*j)/(a3(-1, l)*E9(nv + 1j))
            rr.append(f"j={j}: {complex(ex/px - 1):.1e}")
        pr(f"   nu = {nv}: " + ", ".join(rr))
    # ---- EE vs scalar per harmonic: is R_j^EE / R_j^sc = c4(nu - i j) for one function c4?
    pr("(E) EE rank 3: R_j^EE / R_j^sc (exact):")
    rat = {}
    for j in sorted(D['EE3']):
        rat[j] = sp.factor(sp.simplify(D['EE3'][j]/D['sc3'][j]))
        pr(f"   j = {j:+d}: {rat[j]}")
    pr("   test c4(nu') := ratio at the shifted harmonic: ratio_j(nu) == F(nu - i j) for one F?  F := ratio_{-1}(nu + ... ) shifted:")
    F = rat[-1].subs(nu_s, nu_s - I)                     # ratio_{-1}(nu) = F(nu + i)
    for j in sorted(rat):
        pr(f"     j = {j:+d}: ratio_j(nu) - F(nu - i j) = {sp.simplify(rat[j] - F.subs(nu_s, nu_s - I*j))}")
    pr(f"   F(nu) = {sp.factor(sp.simplify(F))};   k = 2 EE factor c(nu) = 32(4nu^2+1)(4nu^2+9)/((4nu^2+25)(4nu^2+81));  F/c = {sp.factor(sp.simplify(F/cfun(nu_s)))}")
    pr("(H) helicity blindness: R_j^EE(hel = -1) - R_j^EE(hel = +1): " +
       ", ".join(f"j={j}: {sp.simplify(D['EE3m'].get(j, 0) - D['EE3'][j])}" for j in sorted(D['EE3'])))
    pr("(S) EE rank 3 exact R_j:")
    for j in sorted(D['EE3']):
        pr(f"   R_{j:+d} = {D['EE3'][j]}")
    pr("    symmetry R_{-j}(nu) = R_j(-nu): " + ", ".join(f"j={j}: {sp.simplify(D['EE3'][-j] - D['EE3'][j].subs(nu_s, -nu_s))}" for j in (1, 3)))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'stage4_hj.txt'), 'w').write(buf.getvalue())
