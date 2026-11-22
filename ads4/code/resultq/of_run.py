"""
of_run (RESULT_Q_FERMION_HEFT.md, task 1): production run of the O(1/C_T) spinning check of pointwise Lemma A (swapped form, D_LR domain,
CLS-D T-shift) on the exact connected <TTTT> of the free massless Dirac fermion in d = 3 (engine of_fermion.py).
Output: cache/of_check.txt (or cache/of_check_<sections>.txt).   OMP_NUM_THREADS=1 python3 of_run.py [section ...]
sections: val kin ee tens corner pair ctrl gff mix   (default: all)
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import og_free as F
import of_fermion as Fm

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []
T0 = time.time()
CONN = ('1234', '1243', '1324')
SW = {'1234': '1234', '1243': '1324', '1324': '1243'}      # image of each cycle class under 1 <-> 3
VT = Fm.VERT_T
NS, GM, GI = F.null_basis()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def rnull(rng):
    v = rng.normal(size=2); v /= np.linalg.norm(v); s = rng.uniform(0.5, 2.0)
    return np.array([s, s*v[0], s*v[1]])


def norm_ratio(p=1.5, vert=VT, seed=0, n=6):
    """<O O>/flat structure at n random spacelike configurations: (mean, max relative spread)."""
    rng = np.random.default_rng(seed)
    rs = []
    for _ in range(n):
        th = rng.uniform(0, 2*np.pi); x = rng.uniform(0.7, 1.5)*np.array([rng.uniform(-0.5, 0.5), np.cos(th), np.sin(th)])
        z = F.zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2); w = F.zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2)
        rs.append(Fm.twopt(x, z, w, p=p, vert=vert)/Fm.flat_TT(x, z, w, Delta=2*p))
    rs = np.array(rs)
    return complex(np.mean(rs)), float(np.max(np.abs(rs - rs[0]))/abs(rs[0]))


CT_D, _ = norm_ratio()
C2 = CT_D.real**2          # Dirac <TT> = CT_D x (v4c C_T = 1 structure); results /C2 are in v4c C_T = 1 MFT units


def combine(fun, structs=('EE',)):
    """fun(pol2, pol4) -> {class: value}; returns {struct: {class: contraction}} for the rotation-invariant (2,4) structures
    EE = T4^{mn} T2_{mn}, '0000' = T4_{00} T2_{00}, '0m0m' = T4_{0m} T2_0^m (null-basis decomposition, Gram inverse)."""
    K = {'EE': GM, '0000': np.ones((5, 5)), '0m0m': np.array([[F.dot(a, b) for b in NS] for a in NS])}
    M = {}
    for a in range(5):
        for b in range(5):
            for c, v in fun(NS[b], NS[a]).items():
                M.setdefault(c, np.zeros((5, 5), complex))[a, b] = v
    return {s: {c: complex(np.sum((GI @ Mc @ GI)*K[s])) for c, Mc in M.items()} for s in structs}


def links_pol(t1, t3, z1, z3, P, y13, p=1.5, mixed_phase=1.0):
    """gauge links with ARBITRARY polarisations P[i] (controls); positions x1 = -t1 z1, x3 = t3 z3; exact (13) base.
    mixed_phase: constant factor on every mixed link (control kappa = i)."""
    t1 = np.asarray(t1, complex); t3 = np.asarray(t3, complex)
    sh = np.broadcast(t1, t3).shape
    t1 = np.broadcast_to(t1, sh); t3 = np.broadcast_to(t3, sh)
    X = {1: -t1[..., None]*z1, 2: np.broadcast_to(-F.E0, sh + (3,)).astype(complex),
         3: t3[..., None]*z3, 4: np.broadcast_to(F.E0, sh + (3,)).astype(complex)}
    sigma = {pr_: (-1 if pr_ in Fm.MIXED else 1) for pr_ in Fm.PAIRS}
    wover = {}
    for (i, j), tt in (((1, 2), t1), ((3, 4), t3)):
        x = X[i] - X[j]; w = 2*tt - 1 + 0j
        wover[(i, j)] = (w, w**-p, F.dot(P[i], x), F.dot(P[j], x), float(F.dot(P[i], P[j])), Fm.slash(x))
    x13 = X[1] - X[3]
    w = 2*y13*t1*t3; r = (2*y13)**-p*t1**-p*t3**-p
    wover[(1, 3)] = (w, r, F.dot(P[1], x13), F.dot(P[3], x13), float(F.dot(P[1], P[3])), Fm.slash(x13))
    Lk = Fm.links_generic(X, P, Fm.PAIRS, sigma, p, wover)
    if mixed_phase != 1.0:
        for pr_ in Fm.MIXED:
            Lk[pr_] = {mn: mixed_phase*v for mn, v in Lk[pr_].items()}
    return Lk


def E_ctrl(k, phi1, phi3, p2, p4, n=80, q=0, vert=VT, refl=False, classes=CONN):
    """E+ with controls: weight t^(3+q); vertex `vert`; refl: polarisations of the ^+ operators 1, 2 reflected (y -> -y)."""
    z1 = F.zvec(phi1); z3 = F.zvec(phi3); y13 = 2*np.sin((phi1 - phi3)/2)**2
    R = np.diag([1.0, 1.0, -1.0])
    P = {1: R @ z1 if refl else z1, 2: R @ p2 if refl else p2, 3: z3, 4: p4}
    t, w = F.contour(n); T1 = t[:, None]; T3 = t[None, :]
    Lk = links_pol(T1, T3, z1, z3, P, y13)
    G = Fm.correlator(Lk, P, classes, vert)
    T1b, T3b = np.broadcast_arrays(T1, T3)
    h0 = T1b**(3 + q)*T3b**(3 + q)*F.fk(k, T1b, T3b, z1, z3, y13=y13)
    W = w[:, None]*w[None, :]
    return {c: complex(np.sum(W*h0*G[c])) for c in classes}


def resid(A, B, sg):
    sc = max(abs(v) for v in list(A.values()) + list(B.values()))
    return abs(sum(A.values()) + sg*sum(B.values()))/sc


# ----------------------------------------------------------------------------------------------------------------------- (1)
def sec_val():
    import v4d_core as D, v4c_core as V
    pr('(1) VALIDATION of the spinor T-shift rule and of the engine')
    m, sp = norm_ratio()
    pr(f'  (1a) Dirac <T(x,z)T(0,w)> / [((I z).w)^2/(x^2)^3], 6 random spacelike configs: {m.real:.12f} (spread {sp:.0e})  -> T primary, '
       f'C2 = {C2:.6g}')
    for Dp in (0.8, 1.2, 1.5, 2.0, 2.5):
        s1 = norm_ratio(p=Dp + 0.5, vert=[(1.0, 1, 0), (-1.0, 0, 1)], seed=2)[1]
        sp_other = [norm_ratio(p=Dp + 0.5, vert=[(1.0, 1, 0), (a, 0, 1)], seed=2)[1] for a in (-0.9, -0.5, 0.0)]
        pr(f'      GFF fermion Dpsi = {Dp}: O = psibar zslash z.d psi + a (z.d psibar) zslash psi; 2-pt conformality spread at a = -1: {s1:.0e}; '
           f'at a = -0.9, -0.5, 0: {", ".join(f"{x:.2f}" for x in sp_other)}  -> primary iff a = -1')
    rng = np.random.default_rng(5)
    pr('  (1b) 2-pt <T(x3) T(x1^+)>, x1 - x3 past-timelike, mixed-link rule (base -x^2 > 0, coefficient +1) vs the CLS-D embedding')
    pr('       continuation H13^2/|X13|^{Delta+2} (= -1 x flat formula; v4c "-1 per mixed T link"):')
    for Dp in (1.0, 0.8, 1.7):
        p = Dp + 0.5; ratios = []
        cn, _ = norm_ratio(p=p, vert=[(1.0, 1, 0), (-1.0, 0, 1)], seed=3)
        for _ in range(5):
            v = rng.normal(size=2); v *= rng.uniform(0.05, 0.9)/np.linalg.norm(v)
            x = np.array([-1.0, v[0], v[1]])*rng.uniform(0.5, 2)
            z, w = rnull(rng), rnull(rng)
            a = Fm.twopt(x, z, w, p=p, vert=[(1.0, 1, 0), (-1.0, 0, 1)], sg=-1)
            X2 = F.dot(x, x); Iz = z - 2*x*F.dot(x, z)/X2
            emb = (F.dot(Iz, w)*X2)**2/abs(X2)**(2*p + 2)
            ratios.append(a/(cn*emb))
        pr(f'      Dpsi = {Dp}: rule / embedding over 5 random configs: {np.round(np.real(ratios), 13)}  (max |Im| {np.max(np.abs(np.imag(ratios))):.0e})')
    pr('  (1c) 3-pt functions with ^+ points in the gauge configuration (real t1, t3 > 1/2): Wick sum with the rule vs the CLS-D')
    pr('       continuation N_flat / [prod_{same} (x^2)^{p+1} prod_{mixed} |x^2|^{p+1}] (N_flat: the rational numerator, any consistent branch);')
    pr('       generic null polarisations at all points; control kappa = i (phase i per mixed link, the "no T-shift" branch) in brackets:')
    for Dp in (1.0, 0.8, 1.6):
        p = Dp + 0.5; e = p + 1; vert = [(1.0, 1, 0), (-1.0, 0, 1)]
        out = []
        for trial in range(3):
            t1, t3 = rng.uniform(0.6, 2, 2); f1, f3 = rng.uniform(-np.pi, np.pi, 2)
            z1 = F.zvec(f1); z3 = F.zvec(f3)
            X = {1: (-t1*z1).astype(complex), 2: (-F.E0).astype(complex), 3: (t3*z3).astype(complex), 4: F.E0.astype(complex)}
            # generic (not light-ray) polarisations at 1 and 3: with light-ray ones the triangles containing (13) vanish identically
            # (kinematic identity of section 2) and the comparison would be 0/0
            P = {1: rnull(rng), 2: rnull(rng), 3: rnull(rng), 4: rnull(rng)}
            for idx, mixed in (((1, 2, 3), {(1, 3), (2, 3)}), ((1, 3, 4), {(1, 3), (1, 4)})):
                pairs = [(idx[0], idx[1]), (idx[0], idx[2]), (idx[1], idx[2])]

                def three(sig, ph=1.0):
                    Lk = Fm.links_generic(X, P, pairs, {q_: sig.get(q_, 1) for q_ in pairs}, p)
                    for q_ in pairs:
                        if q_ in mixed and ph != 1.0:
                            Lk[q_] = {mn: ph*v for mn, v in Lk[q_].items()}
                    a_, b_, c_ = idx
                    return Fm.loop(Lk, (a_, b_, c_), P, vert) + Fm.loop(Lk, (a_, c_, b_), P, vert)
                s2 = {q_: F.dot(X[q_[0]] - X[q_[1]], X[q_[0]] - X[q_[1]]).real for q_ in pairs}
                N = three({})*np.prod([(s2[q_] + 0j)**e for q_ in pairs])
                cont = N/np.prod([abs(s2[q_])**e if q_ in mixed else (s2[q_] + 0j)**e for q_ in pairs])
                rule = three({q_: -1 for q_ in mixed}); ctl = three({q_: -1 for q_ in mixed}, ph=1j)
                out.append(f'{idx}: {(rule/cont).real:.13f} [{(ctl/cont).real:+.3f}]')
        pr(f'      Dpsi = {Dp}: ' + '; '.join(out))
    pr('  (1d) contour swapped form on the DISCONNECTED (MFT) Dirac terms vs the exact v4d_core closed forms (EE, k = 4, /C2), n = 100:')
    for d in (0.3, 1.3, 2.2, 3.0):
        r = F.EE_combine(lambda p2, p4: Fm.Eplus(4, d/2, -d/2, p2, p4, n=100, classes=('s', 't', 'u')))
        c = np.cos(d)
        pr(f'      dphi = {d}: s = {r["s"].real/C2:+.12e} (exact {-D.RtEE(c):+.12e}, rel {abs(r["s"]/C2 + D.RtEE(c))/abs(D.RtEE(c)):.1e});  '
           f't rel {abs(r["t"]/C2 - D.RtEE(c))/abs(D.RtEE(c)):.1e};  |u|/|t| = {abs(r["u"])/abs(r["t"]):.1e} (exact 0)')
    pr('  (1e) random null (2,4) polarisations, MFT s, t terms vs v4c_core.edge_residue_full (independent Cauchy-FFT code):')
    lo = lambda v: v*np.array([-1, 1, 1])
    for tr in range(3):
        p2, p4 = rnull(rng), rnull(rng); a, b = rng.uniform(-np.pi, np.pi, 2)
        r = Fm.Eplus(4, a, b, p2, p4, n=80, classes=('s', 't', 'u'))
        Ts = V.edge_residue_full(np.array([a]), np.array([b]), 4, 's')[0]
        Tt = V.edge_residue_full(np.array([a]), np.array([b]), 4, 't')[0]
        cs = np.einsum('abcd,a,b,c,d->', Ts, lo(p4), lo(p4), lo(p2), lo(p2))
        ct = np.einsum('abcd,a,b,c,d->', Tt, lo(p4), lo(p4), lo(p2), lo(p2))
        pr(f'      trial {tr}: |s/(C2 R_s) - 1| = {abs(r["s"]/C2/cs - 1):.1e},  |t/(C2 R_t) - 1| = {abs(r["t"]/C2/ct - 1):.1e},  '
           f'|u|/|t| = {abs(r["u"]/r["t"]):.1e}')
    pr('  (1f) connected class 1234, generic component, (phi1, phi3) = (0.7, -1.1): contour radius 0.4 / 0.5 / 0.75 and n = 60 / 100:')
    rng2 = np.random.default_rng(9); p2, p4 = rnull(rng2), rnull(rng2)
    ref = Fm.Eplus(4, 0.7, -1.1, p2, p4, n=100, classes=('1234',))['1234']
    for rad in (0.4, 0.75):
        v = Fm.Eplus(4, 0.7, -1.1, p2, p4, n=100, rad=rad, classes=('1234',))['1234']
        pr(f'      rad = {rad}: 1234 = {v.real/C2:+.13e}  rel. to rad 0.5: {abs(v/ref - 1):.1e}')
    v = Fm.Eplus(4, 0.7, -1.1, p2, p4, n=60, classes=('1234',))['1234']
    pr(f'      rad = 0.5: 1234 = {ref.real/C2:+.13e};  n = 60 vs 100: {abs(v/ref - 1):.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (2)
def sec_kin():
    pr('(2) KINEMATIC IDENTITY in the light-ray gauge: classes 1243 and 1324 contain the (13) link sandwiched between the light-ray')
    pr('    vertices, zslash3 [a zslash1 + b zslash3] zslash1 = 0 (x13 = -t1 z1 - t3 z3 in span(z1, z3), zslash^2 = z^2 = 0), so they vanish')
    pr('    IDENTICALLY for the fermion (any t1, t3, any pol2, pol4).  Raw (un-integrated) correlator at random complex (t1, t3):')
    rng = np.random.default_rng(12)
    mx = {c: 0.0 for c in CONN}
    for _ in range(40):
        t1, t3 = rng.normal(size=2) + 1j*rng.normal(size=2); a, b = rng.uniform(-np.pi, np.pi, 2)
        Lk, P = Fm.links_gauge(np.array([t1]), np.array([t3]), F.zvec(a), F.zvec(b), rnull(rng), rnull(rng), 2*np.sin((a - b)/2)**2)
        G = Fm.correlator(Lk, P)
        for c in CONN:
            mx[c] = max(mx[c], abs(G[c][0])/abs(G['1234'][0]))
    pr(f'    max over 40 points of |class|/|1234|: ' + ', '.join(f'{c}: {v:.1e}' for c, v in mx.items()) + '  (roundoff)')
    pr('    => for the fermion the whole Lemma-A content sits in class 1234 (links 12 [^+^+], 23, 41 [mixed], 34): 1234(a,b) +- 1234(b,a) = 0.')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (3)
def sec_ee():
    pr('(3) ROTATION-INVARIANT (2,4) STRUCTURES (EE = T4.T2, 0000, 0m0m), k = 4, connected, /C2.  For these E(a,b) = E(b,a) by rotation +')
    pr('    reflection, so Lemma A predicts every class to vanish (1234 by Lemma A; 1243, 1324 by the identity of (2)).')
    pr('    reference scale: max over the 25 basis components of |1234| at the same dphi (in brackets).')
    for d in (1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 3.0):
        r = combine(lambda p2, p4: Fm.Eplus(4, d/2, -d/2, p2, p4, n=100), structs=('EE', '0000', '0m0m'))
        scale = max(abs(Fm.Eplus(4, d/2, -d/2, NS[ib], NS[ia], n=100, classes=('1234',))['1234'])
                    for ia in range(5) for ib in range(5))
        row = '  '.join(f'{s}: ' + ' '.join(f'{abs(v[c])/C2:.0e}' for c in CONN) for s, v in r.items())
        pr(f'    dphi = {d:7.0e}: |1234| |1243| |1324| per structure  {row}   [{scale/C2:.2e}]   max ratio {max(abs(v[c]) for v in r.values() for c in ("1234",))/scale:.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (4)
def sec_tens():
    pr('(4) ALL (2,4) tensor components (25 null-basis pairs), connected part, k = 2, 3, 4, 5.  Lemma A: E(a,b) + (-1)^k E(b,a) = 0.')
    for (a, b) in ((0.7, -1.1), (2.6, 0.4), (-2.2, 1.9)):
        for k in (2, 3, 4, 5):
            sg = 1 if k % 2 == 0 else -1
            res, wrong, scl, pcl = [], [], [], {c: [] for c in CONN}
            for ia in range(5):
                for ib in range(5):
                    A = Fm.Eplus(k, a, b, NS[ib], NS[ia], n=80); B = Fm.Eplus(k, b, a, NS[ib], NS[ia], n=80)
                    sc = max(abs(v) for v in list(A.values()) + list(B.values()))
                    tA, tB = sum(A.values()), sum(B.values())
                    res.append(abs(tA + sg*tB)); wrong.append(abs(tA - sg*tB)); scl.append(sc)
                    for c in CONN:
                        pcl[c].append(abs(A[c] + sg*B[SW[c]]))
            S = max(scl)
            nz = [w/s for w, s in zip(wrong, scl) if s > 1e-12*S]
            pr(f'    ({a}, {b}) k = {k}: max |residual|/S = {max(res)/S:.1e}  (S = {S/C2:.3e});  per class ' +
               ' '.join(f'{c}:{max(v)/S:.0e}' for c, v in pcl.items()) +
               f';  wrong symmetry (median over {len(nz)} non-vanishing): {np.median(nz):.2f};  {25 - len(nz)} components vanish identically')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (5)
def sec_corner():
    pr('(5) COLLINEAR CORNER, class 1234 (the only non-zero class), two generic basis components, Phi = (phi1 + phi3)/2 = 0.3:')
    for (ia, ib) in ((0, 2), (3, 1)):
        for d in (1e-4, 1e-3, 1e-2, 0.1, 1.0):
            a, b = 0.3 + d/2, 0.3 - d/2
            A = Fm.Eplus(4, a, b, NS[ib], NS[ia], n=100, classes=('1234',))['1234']
            B = Fm.Eplus(4, b, a, NS[ib], NS[ia], n=100, classes=('1234',))['1234']
            pr(f'    comp ({ia},{ib}) dphi = {d:.0e}: 1234(a,b) = {A.real/C2:+.10e}, dphi*1234 = {A.real/C2*d:+.8e}, '
               f'|1234(a,b) + 1234(b,a)|/|1234| = {abs(A + B)/abs(A):.1e}')
    pr('    => each family ~ dphi^{-1} (the f_4 factor; no (13) link in 1234), i.e. dphi^{+3/2} with the j\'=0 kernel: no corner singularity')
    pr('       at all at O(1/C_T) for the fermion (scalar: dphi^{-3} per family from the twist-1 (13) channel, RESULT_Q_OG_CHECK §3.6).')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (6)
def sec_pair():
    from v4c_wu import Kfun
    pr('(6) PER-nu LIGHT SIDE, swapped form, nu outside, j\'=0 kernel 2 A13^{(2+delta)/2} K_nu(dphi) (as RESULT_Q_OG_CHECK §3.3).')
    pr('    EE: identically 0 class by class (section 3), so I_nu[EE] = 0 for every nu (trivial).  Non-trivial version: a generic constant')
    pr('    tensor component (n_a at 4, n_b at 2), Phi = 0.3 fixed, I_nu[X] = int_{1e-4}^{pi} d(dphi) ker(dphi) X(Phi + dphi/2, Phi - dphi/2),')
    pr('    X = 1234(a,b), 1234(b,a) separately and their sum (Lemma A: sum = 0); 64 Gauss nodes in log dphi; class 1243 = 1324 = 0 by (2).')
    x, wx = np.polynomial.legendre.leggauss(64)
    lo, hi = np.log(1e-4), np.log(np.pi)
    ld = lo + (hi - lo)*(x + 1)/2; wd = wx*(hi - lo)/2
    dph = np.exp(ld); wd = wd*dph
    for (ia, ib) in ((0, 2), (3, 1)):
        Av, Bv = [], []
        for d in dph:
            a, b = 0.3 + d/2, 0.3 - d/2
            Av.append(Fm.Eplus(4, a, b, NS[ib], NS[ia], n=80, classes=('1234',))['1234']/C2)
            Bv.append(Fm.Eplus(4, b, a, NS[ib], NS[ia], n=80, classes=('1234',))['1234']/C2)
        Av, Bv = np.array(Av), np.array(Bv)
        for nu in (1.0, 3.0):
            dl = 0.5 + 1j*nu
            K = np.array([Kfun(d, nu) for d in dph])
            ker = 2*(2 - 2*np.cos(dph))**((2 + dl)/2)*K*wd
            IA, IB, It = np.sum(ker*Av), np.sum(ker*Bv), np.sum(ker*(Av + Bv))
            Iabs = np.sum(np.abs(ker)*np.abs(Av))
            pr(f'    comp ({ia},{ib}) nu = {nu}: I[1234(a,b)] = {IA.real:+.6e}{IA.imag:+.6e}i  I[1234(b,a)] = {IB.real:+.6e}{IB.imag:+.6e}i  '
               f'I[sum] = {It.real:+.1e}{It.imag:+.1e}i  |I[sum]|/|I[1234(a,b)]| = {abs(It)/abs(IA):.1e}  |I[sum]|/int|ker||X| = {abs(It)/Iabs:.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (7)
def sec_ctrl():
    pr('(7) CONTROLS (connected part, random null (2,4) polarisations, two celestial pairs; residual = |E(a,b) +- E(b,a)|/max|class|)')
    rng = np.random.default_rng(7)
    sets = [(rnull(rng), rnull(rng), *rng.uniform(-np.pi, np.pi, 2)) for _ in range(2)]
    for (p2, p4, a, b) in sets:
        pr(f'    pair (phi1, phi3) = ({a:.4f}, {b:.4f})')
        rows = []
        for lab, kw in (('none (a = -1, weight t^3)', {}), ('(a) weight t^4', dict(q=1)), ('(a) weight t^5', dict(q=2)),
                        ('(b) non-primary vertex a = 0', dict(vert=[(1.0, 1, 0), (0.0, 0, 1)])),
                        ('(b) non-primary vertex a = -1/2', dict(vert=[(1.0, 1, 0), (-0.5, 0, 1)])),
                        ('(f) reflected frame at 1^+, 2^+ (T(x^+, P z))', dict(refl=True))):
            out = []
            for k in (4, 5):
                sg = 1 if k == 4 else -1
                A = E_ctrl(k, a, b, p2, p4, **kw); B = E_ctrl(k, b, a, p2, p4, **kw)
                out.append(f'k={k}: {resid(A, B, sg):.1e}')
            rows.append(f'      {lab:48s} ' + '   '.join(out))
        for r_ in rows:
            pr(r_)
        for k in (4, 5):
            sg = 1 if k == 4 else -1
            A = Fm.Eplus(k, a, b, p2, p4); B = Fm.Eplus(k, b, a, p2, p4)
            wgt = {'1234': -1, '1243': -1, '1324': 1}
            sc = max(abs(v) for v in list(A.values()) + list(B.values()))
            e = abs(sum(wgt[c]*A[c] for c in CONN) + sg*sum(wgt[c]*B[c] for c in CONN))/sc
            w_ = abs(sum(A.values()) - sg*sum(B.values()))/sc
            pr(f'      (e) no T-shift (phase i per mixed link -> class weights (-1,-1,+1)), k = {k}: {e:.1e}  [DEGENERATE: only 1234 is non-zero];'
               f'  wrong symmetry: {w_:.2f}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (8)
def sec_gff():
    pr('(8) GENERIC EXPONENTS: O = [psibar psi]_{0,2} of a generalised free Dirac fermion, Delta_O = 2 Dpsi + 1, J = 2, vertex a = -1,')
    pr('    links slash(x) (sigma x^2)^{-Dpsi - 1/2}, light-transform weight t^{2 Dpsi + 1}; connected part; two random tensor/angle sets.')
    rng = np.random.default_rng(11)
    sets = [(rnull(rng), rnull(rng), *rng.uniform(-np.pi, np.pi, 2)) for _ in range(2)]
    vert = [(1.0, 1, 0), (-1.0, 0, 1)]
    for Dp in (0.8, 1.2, 1.5, 2.0, 2.5):
        p = Dp + 0.5; row = []
        for (p2, p4, a, b) in sets:
            for k in (4, 5):
                sg = 1 if k == 4 else -1
                A = Fm.Eplus(k, a, b, p2, p4, n=100, p=p, vert=vert); B = Fm.Eplus(k, b, a, p2, p4, n=100, p=p, vert=vert)
                row.append(f'k={k}: |E|={abs(sum(A.values())):.2e}, res {resid(A, B, sg):.0e}')
        pr(f'    Dpsi = {Dp}: ' + ';  '.join(row))
    pr('    d/dDpsi at Dpsi = 1.5 (central difference, h = 1e-3; log(x_ij^2) in every link):')
    p2, p4, a, b = sets[0]
    for k in (4, 5):
        sg = 1 if k == 4 else -1
        out = []
        for h in (+1e-3, -1e-3):
            p = 2.0 + h
            A = sum(Fm.Eplus(k, a, b, p2, p4, n=100, p=p, vert=vert).values()); B = sum(Fm.Eplus(k, b, a, p2, p4, n=100, p=p, vert=vert).values())
            out.append((A, B))
        dA = (out[0][0] - out[1][0])/2e-3; dB = (out[0][1] - out[1][1])/2e-3
        pr(f'      k = {k}: dE(a,b)/dDpsi = {dA.real:+.6e}{dA.imag:+.1e}i, dE(b,a)/dDpsi = {dB.real:+.6e}, '
           f'|d(E(a,b) {"+" if sg > 0 else "-"} E(b,a))/dDpsi| / |dE(a,b)/dDpsi| = {abs(dA + sg*dB)/abs(dA):.0e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (9)
def sec_mix():
    pr('(9) t4 = 0 MIXTURE (equal C_T per real scalar and per Majorana fermion: n_s = n_M gives the Einstein-like <TTT>, t4 = 0).')
    pr('    E_mix = E_scalar/C2_s + E_Majorana/C2_M per field, with Majorana = Dirac/2 and C2_M = (CT_D/2)^2 (v4c C_T = 1 units each):')
    rng = np.random.default_rng(21)
    C2s = 2.25; C2M = (CT_D.real/2)**2
    for _ in range(3):
        p2, p4 = rnull(rng), rnull(rng); a, b = rng.uniform(-np.pi, np.pi, 2)
        for k in (4, 5):
            sg = 1 if k == 4 else -1
            Sa = sum(F.Eplus(k, a, b, p2, p4).values())/C2s; Sb = sum(F.Eplus(k, b, a, p2, p4).values())/C2s
            Fa = sum(Fm.Eplus(k, a, b, p2, p4).values())/2/C2M; Fb = sum(Fm.Eplus(k, b, a, p2, p4).values())/2/C2M
            mixa, mixb = Sa + Fa, Sb + Fb
            sc = max(abs(Sa), abs(Fa), abs(Sb), abs(Fb))
            pr(f'    ({a:+.3f}, {b:+.3f}) k = {k}: scalar {Sa.real:+.4e}, Majorana {Fa.real:+.4e} (per field, C_T units); '
               f'|mix(a,b) {"+" if sg > 0 else "-"} mix(b,a)|/scale = {abs(mixa + sg*mixb)/sc:.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


if __name__ == '__main__':
    secs = sys.argv[1:] or ['val', 'kin', 'ee', 'tens', 'corner', 'pair', 'ctrl', 'gff', 'mix']
    pr('# of_run: O(1/C_T) spinning check of pointwise Lemma A (swapped form, D_LR, CLS-D T-shift) on the exact free Dirac fermion '
       'connected <TTTT> (of_fermion.py)')
    for s in secs:
        {'val': sec_val, 'kin': sec_kin, 'ee': sec_ee, 'tens': sec_tens, 'corner': sec_corner, 'pair': sec_pair, 'ctrl': sec_ctrl,
         'gff': sec_gff, 'mix': sec_mix}[s]()
    pr(f'# total {time.time() - T0:.0f} s')
    open(os.path.join(CACHE, 'of_check' + ('' if not sys.argv[1:] else '_' + '_'.join(secs)) + '.txt'), 'w').write('\n'.join(L) + '\n')
