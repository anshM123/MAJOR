"""
Parity-odd conserved <T T O_{Delta,J}> structures (d = 3) for the Taylor-mode <TTTT> waves (PARITY_ODD_HEAVY.md sect. 4).

Odd operators (embedding space; 'O', 'Q', 'R', 'S' applied FIRST to the scalar seed, i.e. rightmost in a word; 'U', 'V' last; see odd_words):
   'O' = eps(Z_i, Z_j, P_i, P_j, d/dP_i)   [seed shift (0, 1): P_j-degree +1]
   'Q' = eps(Z_i, Z_j, P_i, P_j, d/dP_j)   [seed shift (1, 0)]
(eps contracted with the inverse metric: eps(..., d/dP) sigma(P.X) = eps(..., X)).  Acting on the Z-independent seed the result is
gauge invariant (eps(Z_i + a P_i, ...) = eps(Z_i, ...)), tangent to the cone, and has spin (1, 1); the parity-even CPPR operators applied
afterwards keep it a valid structure.  Odd words for l = 2:  w + 'O' and w + 'Q', w in the l = 1 even words {AC, AE, BC, BE, H}.
Even operators act on the pseudo-scalars eps(Z1, Z2, P1, P2, X) only through Euler terms (every other replacement repeats a slot).

This module: (i) the odd words' three-point functions in the Euclidean frame at infinity (O_1(x), O_2(0), O_5(inf, zeta)), symbolic in
x, a = z1, b = z2, zeta;  (ii) CHL's parity-odd helicity structures in position space (odd part of arXiv:2102.08160 eq. Tpm,
X^O_opp = i[T^{(-1,-1)} - T^{(1,1)}], X^O_same = i[T^{(-1,1)} - T^{(1,-1)}] in the (s1, s2) labels of chl_xspace), with the SAME dropped
constant as chl_xspace;  (iii) numerical decomposition of X^O into the odd words (and of X^E into the even words, as a check against the
cached production coefficients cache/xwords_*.npz).
"""
import os, sys, functools, pickle, itertools
import sympy as sp, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from embed2 import _mk_engine, op, PAIRS3, Delta, J, D1, D2
from structures2 import words as even_words
import chl_xspace as CX
from ft_helicity import AX, BX, WX, XX, AB, AW, BW
CACHE = os.path.join(HERE, 'cache')
L = 2
EX1, EX2, EP = sp.symbols('EX1 EX2 EP')           # eps(e1, x, e3), eps(e2, x, e3), eps(e1, e2, e3)   (3d, Euclidean)

def odd_words(l=2):
    """'O', 'Q' applied first (rightmost);  'U' = eps(Z1, P1, P2, Z2, d/dZ2) [spin 1 +1], 'V' = eps(Z2, P2, P1, Z1, d/dZ1) [spin 2 +1]
    applied last (leftmost) to an even structure of spins (l-1, l) resp. (l, l-1)."""
    if l != 2: raise NotImplementedError
    out = [w + o for w in ['AC', 'AE', 'BC', 'BE', 'H'] for o in ('O', 'Q')]
    out += ['U' + a + c for a in 'AB' for c in ('CC', 'CE', 'EE')] + ['U' + 'H' + c for c in 'CE']
    out += ['V' + a + c for a in ('AA', 'AB', 'BB') for c in 'CE'] + ['V' + 'H' + a for a in 'AB']
    out += ['H' + c + 'R' for c in 'CE'] + ['H' + a + 'S' for a in 'AB']
    out += [a + c + 'R' for a in 'AB' for c in ('CC', 'CE', 'EE')] + [a + c + 'S' for a in ('AA', 'AB', 'BB') for c in 'CE']
    return out

def odd_words_min(l=2):
    """O/Q words (rank 8) + the four cheapest second-order words"""
    return [w + o for w in ['AC', 'AE', 'BC', 'BE', 'H'] for o in ('O', 'Q')] + ['HCR', 'HER', 'HAS', 'HBS']

def shifts_odd(word):
    s1 = word.count('A') + word.count('E') + word.count('H') + word.count('Q') + word.count('U') + word.count('V')
    s2 = word.count('B') + word.count('C') + word.count('H') + word.count('O') + word.count('U') + word.count('V')
    return s1, s2

# ------------------------------------------------------------------ (i) odd words, three points
ORDER3 = ['Z1', 'Z2', 'P1', 'P2', 'P5', 'Z5']
EPS3 = {}                                           # missing vector -> symbol eps(ORDER3 without it, in ORDER3 order)
for _m in ORDER3:
    EPS3[_m] = sp.Symbol('Eno' + _m)
E5, EZ = EPS3['Z5'], EPS3['P5']                     # eps(Z1, Z2, P1, P2, P5), eps(Z1, Z2, P1, P2, Z5)
EPS_SLOTS = {q: tuple(v for v in ORDER3 if v != m) for m, q in EPS3.items()}
SYM_OF = {v: q for q, v in EPS_SLOTS.items()}

def eps_canon(slots):
    """(sign, symbol) of eps(slots) (5 distinct names from ORDER3), 0 if repeated"""
    if len(set(slots)) < 5: return 0, None
    idx = [ORDER3.index(v) for v in slots]
    sign = 1
    for i in range(5):
        for j in range(i + 1, 5):
            if idx[i] > idx[j]: sign = -sign
    return sign, SYM_OF[tuple(sorted(slots, key=ORDER3.index))]

def _engine3():
    e, S = _mk_engine(PAIRS3)
    old = e.dsym
    for q in EPS_SLOTS: e.syms[q] = None
    def dsym(q, V, W):
        if q in EPS_SLOTS:
            sl = list(EPS_SLOTS[q])
            if W not in sl: return 0
            sl[sl.index(W)] = V
            sg, qq = eps_canon(sl)
            return sg*qq if qq is not None else 0
        return old(q, V, W)
    e.dsym = dsym
    return e, S

def odd_op3(e, S, Phi, which):
    """(1/K) eps(Z1, Z2, P1, P2, d/dP_d) (K Phi), d = 1 ('O') or 2 ('Q'); Phi may not contain Z1, Z2 dependence (applied first)"""
    Pd = 'P1' if which == 'O' else 'P2'
    out = 0
    tot = sp.together(sp.log(sp.Symbol('dummy')))   # placeholder (unused)
    for q, pr in e.syms.items():
        if q in EPS_SLOTS or pr is None: continue
        A, B = pr
        if Pd not in (A, B): continue
        other = B if A == Pd else A
        fac = 2 if str(q).startswith('s') else 1
        sg, qq = eps_canon(['Z1', 'Z2', 'P1', 'P2', other])
        if qq is None: continue
        eps = sg*qq
        der = sp.diff(Phi, q) + sp.diff(e.logK, q)*Phi
        out += fac*eps*der
    return sp.expand(out)

def rot_op3(e, S, Phi, which):
    """'U': eps(Z1, P1, P2, Z2, d/dZ2) Phi = sum_q dPhi/dq eps(Z1,P1,P2,Z2, partner of Z2 in q);  eps(Z1,P1,P2,Z2,X) = eps(Z1,Z2,P1,P2,X).
       'V': eps(Z2, P2, P1, Z1, d/dZ1) Phi;  eps(Z2,P2,P1,Z1,X) = eps(Z1,Z2,P1,P2,X) (even permutation).  K has no Z1, Z2."""
    Zd = 'Z2' if which == 'U' else 'Z1'
    out = 0
    for q, pr in e.syms.items():
        if q in EPS_SLOTS or pr is None: continue
        A, B = pr
        if Zd not in (A, B): continue
        other = B if A == Zd else A
        sl = ['Z1', 'P1', 'P2', 'Z2', other] if which == 'U' else ['Z2', 'P2', 'P1', 'Z1', other]
        sg, qq = eps_canon(sl)
        if qq is None: continue
        out += sg*qq*sp.diff(Phi, q)
    return sp.expand(out)

def two_op3(e, S, which):
    """(1/K) eps(Z_s, P1, P2, d/dP1, d/dP2) K on the seed (Phi = 1), s = 1 ('R') or 2 ('S')"""
    Zs = 'Z1' if which == 'R' else 'Z2'
    out = 0
    items = [(q, pr) for q, pr in e.syms.items() if q not in EPS_SLOTS and pr is not None]
    K = sp.exp(e.logK)
    for q1, (A1, B1) in items:
        if 'P1' not in (A1, B1): continue
        X = B1 if A1 == 'P1' else A1; f1 = 2 if str(q1).startswith('s') else 1
        for q2, (A2, B2) in items:
            if 'P2' not in (A2, B2): continue
            Y = B2 if A2 == 'P2' else A2; f2 = 2 if str(q2).startswith('s') else 1
            sg, qq = eps_canon([Zs, 'P1', 'P2', X, Y])
            if qq is None: continue
            d2 = sp.diff(e.logK, q1, q2) + sp.diff(e.logK, q1)*sp.diff(e.logK, q2)     # (1/K) d^2 K
            out += f1*f2*sg*qq*d2
    return sp.expand(out)

@functools.lru_cache(None)
def word3_odd(word):
    """Q_w = (w applied to S3(D1+sh1, D2+sh2)) / S3(D1+sh1, D2+sh2), symbolic in the 3-pt symbols, E5, EZ, Delta, J, D1, D2"""
    e, S = _engine3()
    v5 = (S['z51']*S['s25'] - S['z52']*S['s15'])/2
    sh1, sh2 = shifts_odd(word)
    d1, d2 = D1 + sh1, D2 + sh2
    al = (d1 + d2 - Delta + J)/2; be = (Delta + J + d1 - d2)/2; ga = (Delta + J - d1 + d2)/2
    e.logK = J*sp.log(v5) - al*sp.log(S['s12']) - be*sp.log(S['s15']) - ga*sp.log(S['s25'])
    if word[-1] in 'OQ':
        Phi = odd_op3(e, S, sp.Integer(1), word[-1])
        rest = word[:-1]
    elif word[-1] in 'RS':
        Phi = two_op3(e, S, word[-1])
        rest = word[:-1]
    else:
        Phi = sp.Integer(1); rest = word
    for kind in reversed(rest):
        if kind in 'UV':
            Phi = rot_op3(e, S, Phi, kind)
        else:
            Phi = op(e, Phi, kind, 1, 2)
    return sp.together(Phi), S, al

def odd_words_inf(l=2):
    """F_w / [(-zeta.x/2)^J |x|^{-2 al_odd}] in the Euclidean frame at infinity, as expressions in AX, BX, WX, XX, AB, AW, BW, E5, EZ
    (E5 = eps(Z1,Z2,P1,P2,P5), EZ = eps(Z1,Z2,P1,P2,Z5) evaluated numerically later).  Delta1 = Delta2 = l + 1."""
    fn = os.path.join(CACHE, 'odd3_inf_l%d_v4.pkl' % l)
    out = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    for w in odd_words(l):
        if w in out: continue
        Q, S, al = word3_odd(w)
        isub = {S['s12']: -XX, S['s15']: -1, S['s25']: -1, S['z12']: -AX, S['z21']: BX, S['z15']: 0, S['z25']: 0,
                S['z51']: WX, S['z52']: 0, S['z1z2']: AB, S['z1z5']: AW, S['z2z5']: BW}
        out[w] = sp.together(Q.subs({D1: l + 1, D2: l + 1}).xreplace(isub))
        pickle.dump(out, open(fn, 'wb'))
        print('   odd word', w, 'built', flush=True)
    return out

def eps_all_inf(x, a, b, zeta):
    dd = lambda u, v: sum(u[i]*v[i] for i in range(3))
    V = {'P1': [1, dd(x, x)] + list(x), 'P2': [1, 0, 0, 0, 0], 'P5': [0, 1, 0, 0, 0], 'Z1': [0, 2*dd(x, a)] + list(a),
         'Z2': [0, 0] + list(b), 'Z5': [0, 0] + list(zeta)}
    return {q: mp.det(mp.matrix([V[v] for v in EPS_SLOTS[q]])) for q in EPS_SLOTS}

def eps5_inf(x, a, b, zeta):
    """E5, EZ in the Euclidean frame at infinity: P1 = (1, x^2, x), P2 = (1, 0, 0), P5 = (0, 1, 0), Z1 = (0, 2x.a, a), Z2 = (0, 0, b),
    Z5 = (0, 0, zeta); eps = det of the components in the basis (+, -, mu)"""
    dd = lambda u, v: sum(u[i]*v[i] for i in range(3))
    P1 = [1, dd(x, x)] + list(x); P2 = [1, 0, 0, 0, 0]; P5 = [0, 1, 0, 0, 0]
    Z1 = [0, 2*dd(x, a)] + list(a); Z2 = [0, 0] + list(b); Z5 = [0, 0] + list(zeta)
    return mp.det(mp.matrix([Z1, Z2, P1, P2, P5])), mp.det(mp.matrix([Z1, Z2, P1, P2, Z5]))

# ------------------------------------------------------------------ (ii) CHL odd helicity structures, position space
def _deriv_x(expr, which):
    if which == 'a':
        rules = {AX: 0, BX: AB, WX: AW, XX: 2*AX, EX1: 0, EX2: -EP}
    else:
        rules = {AX: AB, BX: 0, WX: BW, XX: 2*BX, EX1: EP, EX2: 0}
    return sum(sp.diff(expr, v)*r for v, r in rules.items() if r != 0)

def _deriv_eps(expr, k):
    """eps(e_k, d_x, e3) on a function of the invariants (no EX yet)"""
    rules = {AX: (0 if k == 1 else -EP), BX: (EP if k == 1 else 0), XX: 2*(EX1 if k == 1 else EX2)}
    assert not expr.has(EX1) and not expr.has(EX2)
    return sum(sp.diff(expr, v)*r for v, r in rules.items() if r != 0)

def momentum_structure_odd(l, kind):
    """odd part:  i [T(s) - T(s')]  with (s, s') = ((-1,-1), (1,1)) for opp, ((-1,1), (1,-1)) for same.
    Returns list of (k, a, b, kk, s2) : coefficient for  B_k P1^a P2^b P3^(J+kk) PP^(s2/2)  (s2 = twice the PP power, odd)."""
    P1, P2, P3, PP = CX.P1, CX.P2, CX.P3, CX.PP
    A1 = AW - P1*P3/PP; A2 = BW - P2*P3/PP
    G11 = sp.Matrix([[0, P1, AW], [P1, PP, P3], [AW, P3, 0]]).det()
    G22 = sp.Matrix([[0, P2, BW], [P2, PP, P3], [BW, P3, 0]]).det()
    b1, b2 = sp.symbols('b1 b2')
    (sa, sb) = ((-1, -1), (1, 1)) if kind == 'opp' else ((-1, 1), (1, -1))
    term = lambda s: sp.expand((A1 + s[0]*sp.I*b1)**l*(A2 + s[1]*sp.I*b2)**l)
    tot = sp.expand(sp.I*(term(sa) - term(sb)))
    poly = sp.Poly(tot, b1, b2)
    out = {}
    for (i, j), c in poly.terms():
        assert (i + j) % 2 == 1
        if i % 2 == 1:
            k, rest = 1, G11**((i - 1)//2)*G22**(j//2)
        else:
            k, rest = 2, G11**(i//2)*G22**((j - 1)//2)
        # b1^i b2^j = B_k * rest * PP^{-(i+j)/2}
        red = sp.expand(c*rest*P3**(-2*l))
        for t in sp.Add.make_args(red):
            pw = t.as_powers_dict()
            a_, b_, kk, s = int(pw.get(P1, 0)), int(pw.get(P2, 0)), int(pw.get(P3, 0)), int(pw.get(PP, 0))
            cc = t/(P1**a_*P2**b_*P3**kk*PP**s)
            key = (k, a_, b_, kk, 2*s - (i + j))
            out[key] = out.get(key, 0) + cc
    return out

def xspace_odd(l, Dv, Jv, kind, dps=40):
    """x-space X^O_kind (CHL frame), up to the SAME common factor as chl_xspace.xspace_structure (reference s0, n0 = J):
    expression in AX, BX, WX, XX, AB, AW, BW, EX1, EX2, EP with mpmath-evaluated constants (Dv may be irrational)."""
    mp.mp.dps = dps
    Dv = mp.mpf(Dv)
    beta = 2*(2*l + 1) - (Dv + Jv)
    s0 = (beta - 3)/2
    terms = momentum_structure_odd(l, kind)
    acc = 0
    for (k, a, b, kk, s2), c in terms.items():
        n = Jv + kk; ds = mp.mpf(s2)/2; ss = s0 + ds
        Kr = [1, 1j, -1, -1j][kk % 4]*mp.mpf(2)**(2*ds + kk)*mp.gamma(ss + n + mp.mpf(3)/2)/mp.gamma(s0 + Jv + mp.mpf(3)/2)*mp.gamma(-s0)/mp.gamma(-ss)
        Kr = sp.Float(mp.re(Kr), dps) + sp.I*sp.Float(mp.im(Kr), dps)
        ex_xx = -ss - mp.mpf(3)/2 - n
        # f = Kr WX^n XX^{ex_xx}; B_k -> -i eps(e_k, d, e3); p.e1 -> -i e1.d, p.e2 -> -i e2.d.  Track XX power symbolically: XX^(ex_xx) = XX^t
        t = sp.Symbol('t')
        f = WX**n*XX**t
        f = -sp.I*_deriv_eps(f, k)
        for _ in range(a): f = -sp.I*_deriv_x(f, 'a')
        for _ in range(b): f = -sp.I*_deriv_x(f, 'b')
        f = sp.expand(f)
        acc += sp.expand(c*Kr*f.subs(t, sp.Float(ex_xx, dps)))
    return acc

# ------------------------------------------------------------------ (iii) numerical decomposition
def _rand_null(rng):
    """random complex null 3-vector (Euclidean)"""
    u = [mp.mpf(rng.uniform(-1, 1)) for _ in range(3)]; v = [mp.mpf(rng.uniform(-1, 1)) for _ in range(3)]
    # make v orthogonal to u and of the same norm: n = u + i v with u.v = 0, |u| = |v|
    uu = sum(x*x for x in u); uv = sum(x*y for x, y in zip(u, v))
    v = [y - uv/uu*x for x, y in zip(u, v)]
    vv = sum(y*y for y in v); sc = mp.sqrt(uu/vv)
    return [x + 1j*y*sc for x, y in zip(u, v)]

def _eps3(u, v, w):
    return mp.det(mp.matrix([u, v, w]))

def _inv(x, a, b, z):
    d = lambda u, v: sum(p*q for p, q in zip(u, v))
    return {AX: d(a, x), BX: d(b, x), WX: d(z, x), XX: d(x, x), AB: d(a, b), AW: d(a, z), BW: d(b, z)}

def configs(n, seed=11):
    import random
    rng = random.Random(seed); out = []
    for _ in range(n):
        x = [mp.mpf(rng.uniform(-1, 1)) for _ in range(3)]
        out.append((x, _rand_null(rng), _rand_null(rng), _rand_null(rng)))
    return out

def _num(expr, sub):
    f = sp.lambdify(list(sub.keys()), expr, 'mpmath')
    return f(*sub.values())

def _rows(Dv, Jv, Tf, ws, fW, al, cfgs):
    syms = [AX, BX, WX, XX, AB, AW, BW]
    rows, rhs = [], []
    for (x, a, b, z) in cfgs:
        I = _inv(x, a, b, z); y = [-c for c in x]; Iy = _inv(y, a, b, z)
        ea = eps_all_inf(y, a, b, z)
        rows.append([(-Iy[WX]/2)**Jv*Iy[XX]**(-(2*(L + 1) + sum(shifts_odd(w)) - mp.mpf(Dv) + Jv)/2)*fW[w](*[Iy[s] for s in syms], *[ea[q] for q in EPS3.values()], mp.mpf(Dv), Jv) for w in ws])
        rhs.append(Tf(*[I[s] for s in syms], _eps3(a, x, z), _eps3(b, x, z), _eps3(a, b, z)))
    return mp.matrix(rows), mp.matrix(rhs)

def decompose_odd(Dv, Jv, kind, npts=60, dps=60, words=None):
    """real coefficients c_r with  X^O_kind(x) = i sum_w c_r[w] Fn_w(-x)  (CHL frame x, our frame y = -x), as in t8_kappa_exact.decompose.
    Least squares on npts random configurations (real unknowns, [Re; Im] stacking), SVD with the exact null space dropped; validated on
    fresh configurations.  Returns (words, c (complex, = i c_r), fit residual, fresh residual, rank)."""
    mp.mp.dps = dps
    T = xspace_odd(L, Dv, Jv, kind, dps)
    ws = words or odd_words(L); ex = odd_words_inf(L)
    al = (2*(L + 1) + 3 - mp.mpf(Dv) + Jv)/2
    syms = [AX, BX, WX, XX, AB, AW, BW]
    Tf = sp.lambdify(syms + [EX1, EX2, EP], T, 'mpmath')
    fW = {w: sp.lambdify(syms + list(EPS3.values()) + [Delta, J], ex[w], 'mpmath') for w in ws}
    A, bvec = _rows(Dv, Jv, Tf, ws, fW, al, configs(npts, seed=11))
    nr = A.rows; Ar = mp.matrix(2*nr, A.cols); br = mp.matrix(2*nr, 1)
    for i in range(nr):
        t = -1j*bvec[i]
        br[i] = mp.re(t); br[nr + i] = mp.im(t)
        for j in range(A.cols):
            Ar[i, j] = mp.re(A[i, j]); Ar[nr + i, j] = mp.im(A[i, j])
    sc = [mp.norm(Ar.column(j)) for j in range(Ar.cols)]
    for j in range(Ar.cols):
        for i in range(Ar.rows): Ar[i, j] = Ar[i, j]/sc[j]
    U, S, V = mp.svd_r(Ar)
    smax = max(S[i] for i in range(len(S)))
    keep = [i for i in range(len(S)) if S[i] > smax*mp.mpf(10)**(-(dps*2)//3)]
    cr = mp.matrix(Ar.cols, 1)
    for i in keep:
        coef = (U.column(i).T*br)[0]/S[i]
        for j in range(Ar.cols): cr[j] += coef*V[i, j]
    for j in range(Ar.cols): cr[j] = cr[j]/sc[j]
    c = [1j*cr[j] for j in range(Ar.cols)]
    res = mp.norm(A*mp.matrix(c) - bvec)/mp.norm(bvec)
    A2, b2 = _rows(Dv, Jv, Tf, ws, fW, al, configs(8, seed=1234))
    res2 = mp.norm(A2*mp.matrix(c) - b2)/mp.norm(b2)
    return ws, c, res, res2, len(keep)

def check_even(Dv, Jv, kind, npts=20, dps=40):
    """residual of the cached production coefficients (xwords) against CHL's even structure, with the same frame conventions"""
    mp.mp.dps = dps
    from tmoments_l2 import helicity_words
    from ft_helicity import R_inv
    hw = helicity_words(Dv, Jv)[kind]
    T = CX.xspace_structure(L, sp.Rational(Dv), sp.Integer(Jv), kind)
    ws = even_words(L)
    al = (4*L + 2 - mp.mpf(Dv) + Jv)/2
    syms = [AX, BX, WX, XX, AB, AW, BW]
    fT = sp.lambdify(syms, T, 'mpmath')
    Rw = {w: sp.lambdify(syms, R_inv(L, {w: 1}, sp.Rational(Dv), sp.Integer(Jv)), 'mpmath') for w in ws}
    num = den = 0
    for (x, a, b, z) in configs(npts):
        I = _inv(x, a, b, z); y = [-c for c in x]; Iy = _inv(y, a, b, z)
        pre = (-Iy[WX]/2)**Jv*Iy[XX]**(-al)
        lhs = sum(complex(hw[w])*pre*Rw[w](*[Iy[s] for s in syms]) for w in ws)
        tv = fT(*[I[s] for s in syms])
        num += abs(lhs - tv)**2; den += abs(tv)**2
    return mp.sqrt(num/den)

# ------------------------------------------------------------------ conservation of the odd words (physical space, Lorentzian as structures2)
def conservation_odd(Dval, Jval, npt=2):
    """rows: divergence equations at x1 and x2 at rational points; columns: odd words.  Returns (matrix, words)."""
    import random
    from structures2 import _phys, _div_spin2, eta
    ws = odd_words(L)
    _, S, _ = word3_odd(ws[0])
    X1, X2, Z1, Z2, sub = _phys(S)
    x5 = sp.Matrix([sp.Rational(9, 4), sp.Rational(-2, 3), sp.Rational(3, 7)]); n = sp.Matrix([5, 3, 4])
    def emb_P(x): return [1, -x[0]**2 + x[1]**2 + x[2]**2] + list(x)
    def emb_Z(x, z): return [0, 2*(-x[0]*z[0] + x[1]*z[1] + x[2]*z[2])] + list(z)
    P1, P2, P5 = emb_P(X1), emb_P(X2), emb_P(x5); Zv1, Zv2, Zv5 = emb_Z(X1, Z1), emb_Z(X2, Z2), emb_Z(x5, n)
    e5 = sp.Matrix([Zv1, Zv2, P1, P2, P5]).det(method='berkowitz'); ez = sp.Matrix([Zv1, Zv2, P1, P2, Zv5]).det(method='berkowitz')
    cols = []
    rng = random.Random(3)
    pts = [({X1[i]: sp.Rational(rng.randint(-9, 9), 7) for i in range(3)} | {X2[i]: sp.Rational(rng.randint(-9, 9), 5) for i in range(3)},
            {Z1[i]: sp.Rational(rng.randint(-9, 9), 4) for i in range(3)} | {Z2[i]: sp.Rational(rng.randint(-9, 9), 3) for i in range(3)})
           for _ in range(npt)]
    allrows = None
    for w in ws:
        Q, S_, al = word3_odd(w)
        sh1, sh2 = shifts_odd(w)
        d1, d2 = 3 + sh1, 3 + sh2
        alw = (d1 + d2 - Delta + J)/2; be = (Delta + J + d1 - d2)/2; ga = (Delta + J - d1 + d2)/2
        v5 = (S['z51']*S['s25'] - S['z52']*S['s15'])/2
        logM = (J*sp.log(v5) - alw*sp.log(S['s12']) - be*sp.log(S['s15']) - ga*sp.log(S['s25']))
        num = {Delta: Dval, J: Jval, D1: 3, D2: 3}
        F = Q.subs(num).subs(sub).subs({E5: e5, EZ: ez})
        lm = logM.subs(num).subs(sub)
        eqs = []
        for (xp, zp) in pts:
            for (X, Z, Zo) in ((X1, Z1, Z2), (X2, Z2, Z1)):
                dv = _div_spin2(F, lm, X, Z)
                for d in dv:
                    d = d.subs(xp)
                    for zz in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]):
                        eqs.append(sp.nsimplify(d.subs({Zo[i]: zz[i] for i in range(3)}).subs(zp)))
        cols.append(eqs)
    return sp.Matrix([[cols[k][i] for k in range(len(ws))] for i in range(len(cols[0]))]), ws
