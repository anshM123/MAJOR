"""
og_run (RESULT_Q_OG_CHECK.md): production run of the O(1/C_T) spinning check of the pointwise Lemma A (swapped form, D_LR) on the exact
free-scalar connected <TTTT> (engine: og_free.py).  Output: cache/og_check.txt.  Runtime ~ 10 min, 1 core, < 1 GB.
    OMP_NUM_THREADS=1 python3 og_run.py [section ...]      (sections: val ee tens ctrl pair gff; default all)
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import og_free as F

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []
C2 = 2.25            # (C_T^free / C_T^{v4c})^2: <TT>_free = 3/2 x the C_T = 1 structure of v4c_core
SW = {'1234': '1234', '1243': '1324', '1324': '1243'}      # image of each cycle class under 1 <-> 3
CONN = ('1234', '1243', '1324')
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def rnull(rng):
    v = rng.normal(size=2); v /= np.linalg.norm(v); s = rng.uniform(0.5, 2.0)
    return np.array([s, s*v[0], s*v[1]])


def E_custom(k, phi1, phi3, p2, p4, n=80, q=0, fname=None, classes=CONN):
    """E+ with optional light-transform weight t^(3+q) and generic subtraction factor names f0, f1 (controls)."""
    z1 = F.zvec(phi1); z3 = F.zvec(phi3); y13 = 2*np.sin((phi1 - phi3)/2)**2
    t, w = F.contour(n); T1 = t[:, None]; T3 = t[None, :]
    Lk = F.links(T1, T3, z1, z3, p2, p4, y13=y13); G = F.correlator(Lk, classes)
    T1b, T3b = np.broadcast_arrays(T1, T3)
    if fname is None:
        f = F.fk(k, T1b, T3b, z1, z3, y13=y13)
    else:
        f4 = F.fk(2, T1b, T3b, z1, z3, y13=y13)          # f2 = v' - u'
        f = {'f1': f4**2, 'f0': f4**3}[fname]
    h0 = T1b**(3 + q)*T3b**(3 + q)*f
    W = w[:, None]*w[None, :]
    return {c: complex(np.sum(W*h0*G[c])) for c in classes}


# ----------------------------------------------------------------------------------------------------------------------- (1)
def sec_val():
    import v4d_core as D, v4c_core as V
    pr('(1) VALIDATION')
    rng = np.random.default_rng(1)
    rat = []
    for _ in range(5):
        X = rng.normal(size=3); X[0] *= 0.3
        z = F.zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2); w = F.zvec(rng.uniform(0, 6))*rng.uniform(0.5, 2)
        X2 = F.dot(X, X); r = X2**-0.5
        Dm = lambda m, n: F.DFUN[(m, n)](r, X2, F.dot(z, X), F.dot(w, X), F.dot(z, w), 1.0)
        tot = sum(c1*c2*(Dm(a, c)*Dm(b, d) + Dm(a, d)*Dm(b, c)) for c1, (a, b) in F.VERTEX for c2, (c, d) in F.VERTEX)
        Iz = z - 2*X*F.dot(X, z)/X2
        rat.append(tot/(F.dot(Iz, w)**2/X2**3))
    pr(f'  (1a) <T(x,z)T(y,w)>_free / [((I z).w)^2/(x^2)^3] at 5 random spacelike configurations: {np.round(rat, 12)}  '
       f'(conformal primary, C = 3/2 => C2 = 2.25)')
    pr('  (1b) contour swapped form on the DISCONNECTED (MFT) free-field terms vs the exact v4d_core/v4c_mft closed forms (EE, k = 4, /C2):')
    for d in (0.3, 1.3, 2.2, 3.0):
        r = F.EE_combine(lambda p2, p4: F.Eplus(4, d/2, -d/2, p2, p4, n=80, classes=('s', 't', 'u')))
        c = np.cos(d)
        pr(f'      dphi = {d}: s = {r["s"].real/C2:+.12e} (exact {-D.RtEE(c):+.12e}, rel {abs(r["s"]/C2 + D.RtEE(c))/abs(D.RtEE(c)):.1e});  '
           f't = {r["t"].real/C2:+.12e} (rel {abs(r["t"]/C2 - D.RtEE(c))/abs(D.RtEE(c)):.1e});  |u|/|t| = {abs(r["u"])/abs(r["t"]):.1e} (exact 0)')
    pr('  (1c) random null (2,4) polarisations, MFT s, t terms vs v4c_core.edge_residue_full (Cauchy-FFT, independent code):')
    lo = lambda v: v*np.array([-1, 1, 1])
    for tr in range(3):
        p2, p4 = rnull(rng), rnull(rng); a, b = rng.uniform(-np.pi, np.pi, 2)
        r = F.Eplus(4, a, b, p2, p4, n=80, classes=('s', 't', 'u'))
        Ts = V.edge_residue_full(np.array([a]), np.array([b]), 4, 's')[0]
        Tt = V.edge_residue_full(np.array([a]), np.array([b]), 4, 't')[0]
        cs = np.einsum('abcd,a,b,c,d->', Ts, lo(p4), lo(p4), lo(p2), lo(p2))
        ct = np.einsum('abcd,a,b,c,d->', Tt, lo(p4), lo(p4), lo(p2), lo(p2))
        pr(f'      trial {tr}: |s/(C2 R_s) - 1| = {abs(r["s"]/C2/cs - 1):.1e},  |t/(C2 R_t) - 1| = {abs(r["t"]/C2/ct - 1):.1e},  '
           f'|u|/|t| = {abs(r["u"]/r["t"]):.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (2)
def sec_ee():
    pr('(2) CONNECTED (O(1/C_T)) free-scalar <TTTT>, EE (2,4) structure, k = 4: E+(dphi) per cycle class (units: v4c C_T = 1 MFT, i.e. /C2).')
    pr('    Lemma A for EE (E+(a,b) = E+(b,a) by rotation + reflection) predicts E+ == 0 identically: class 1234 == 0, 1243 == -1324.')
    pr('    dphi      | 1234              | 1243                 | 1324                 | SUM                 | |SUM|/|1243| | n=60 vs 100 drift of 1243')
    for d in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        r = F.EE_combine(lambda p2, p4: F.Eplus(4, d/2, -d/2, p2, p4, n=100))
        r60 = F.EE_combine(lambda p2, p4: F.Eplus(4, d/2, -d/2, p2, p4, n=60, classes=('1243',)))
        s = sum(r.values())
        pr(f'    {d:8.1e}  | {r["1234"].real/C2:+.3e}        | {r["1243"].real/C2:+.12e} | {r["1324"].real/C2:+.12e} | '
           f'{s.real/C2:+.2e}{s.imag/C2:+.2e}i | {abs(s/r["1243"]):.1e}      | {abs(r60["1243"]/r["1243"] - 1):.1e}')
    pr('    small-dphi law of one class: 1243 * dphi^3 /C2 -> ' +
       ', '.join(f'{(F.EE_combine(lambda p2, p4: F.Eplus(4, d/2, -d/2, p2, p4, n=80, classes=("1243",)))["1243"].real/C2*d**3):.6e}'
                 for d in (1e-4, 1e-3, 1e-2)) + '  (dphi = 1e-4, 1e-3, 1e-2): class ~ dphi^-3; with the celestial kernel dphi^{5/2}: dphi^{-1/2},')
    pr('    = V4DPRIME homogeneity lemma for a twist-1 (13) exchange (corner integrand dphi^{tau - 3/2}, tau = 1): integrable.')
    pr('    contour-deformation (Cauchy) check at dphi = 1: radius 0.4 / 0.5 / 0.75 of the circle through t = 0:')
    for rad in (0.4, 0.5, 0.75):
        r = F.EE_combine(lambda p2, p4: F.Eplus(4, 0.5, -0.5, p2, p4, n=100, rad=rad))
        pr(f'      rad = {rad}: 1243 = {r["1243"].real/C2:+.12e}, 1324 = {r["1324"].real/C2:+.12e}, 1234 = {abs(r["1234"])/C2:.1e}, '
           f'|SUM|/|1243| = {abs(sum(r.values())/r["1243"]):.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (3)
def sec_tens():
    pr('(3) ALL (2,4) tensor components (25 null-basis pairs n_a at 4, n_b at 2; span of all constant (2,4) tensors), connected part.')
    pr('    Lemma A: k even (2, 4): E(a,b) + E(b,a) = 0;  k odd (3, 5): E(a,b) - E(b,a) = 0.  Per class: X(a,b) +- SW[X](b,a) = 0.')
    ns, Gm, Gi = F.null_basis()
    for (a, b) in ((0.7, -1.1), (2.6, 0.4)):
        for k in (2, 3, 4, 5):
            sg = 1 if k % 2 == 0 else -1
            res, wrong, scl, pcl = [], [], [], {c: [] for c in CONN}
            for ia in range(5):
                for ib in range(5):
                    A = F.Eplus(k, a, b, ns[ib], ns[ia], n=80); B = F.Eplus(k, b, a, ns[ib], ns[ia], n=80)
                    sc = max(abs(v) for v in list(A.values()) + list(B.values()))
                    tA, tB = sum(A.values()), sum(B.values())
                    res.append(abs(tA + sg*tB)); wrong.append(abs(tA - sg*tB)); scl.append(sc)
                    for c in CONN:
                        pcl[c].append(abs(A[c] + sg*B[SW[c]]))
            S = max(scl)          # global scale: residuals are normalised by the largest class value over all 25 components
            nz = [w/s for w, s in zip(wrong, scl) if s > 1e-12*S]
            pr(f'    (phi1, phi3) = ({a}, {b}), k = {k}: max over 25 components of |Lemma-A residual|/S = {max(res)/S:.1e}  '
               f'(S = max|class| over components = {S/C2:.2e});  per-class max/S: ' + ' '.join(f'{c}:{max(v)/S:.1e}' for c, v in pcl.items()) +
               f';  wrong-symmetry |E(a,b) -+ E(b,a)|/max|class|, median over the {len(nz)} non-vanishing components: {np.median(nz):.2f};  '
               f'{25 - len(nz)} components vanish identically (max|class| < 1e-12 S)')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (4)
def sec_ctrl():
    pr('(4) CONTROLS (random null (2,4) polarisations, one celestial pair; residual = |E(a,b) +- E(b,a)| / max|class|, connected part)')
    rng = np.random.default_rng(7)
    p2, p4 = rnull(rng), rnull(rng); a, b = rng.uniform(-np.pi, np.pi, 2)
    pr(f'    pair (phi1, phi3) = ({a:.4f}, {b:.4f})')

    def resid(A, B, sg):
        sc = max(abs(v) for v in list(A.values()) + list(B.values()))
        return abs(sum(A.values()) + sg*sum(B.values()))/sc
    pr('  (4a) light-transform weight t^(3+q) (conformal weight Delta + J - 2 = 3 is q = 0):')
    for q in (0, 1, 2):
        for k in (4, 5):
            sg = 1 if k == 4 else -1
            pr(f'      q = {q}, k = {k}: residual = {resid(E_custom(k, a, b, p2, p4, q=q), E_custom(k, b, a, p2, p4, q=q), sg):.1e}')
    pr('  (4b) non-primary "T" (improvement xi != 1/8: microcausal and Wightman-analytic, but not a conformal primary):')
    keep = list(F.VERTEX)
    for xi in (0.0, 0.125, 0.25):
        F.VERTEX[:] = [(1 - 2*xi, (1, 1)), (-2*xi, (0, 2))]
        for k in (4, 5):
            sg = 1 if k == 4 else -1
            pr(f'      xi = {xi}, k = {k}: residual = {resid(E_custom(k, a, b, p2, p4), E_custom(k, b, a, p2, p4), sg):.1e}')
    F.VERTEX[:] = keep
    pr('  (4c) Regge probe: subtractions with fewer powers (f1 = (v\'-u\')^2: "k = 1"; f0 = (v\'-u\')^3: "k = 0"); CLS-D guarantee needs k > J0 (= 1 here):')
    for fname, sg in (('f0', 1), ('f1', -1)):
        for n in (40, 80, 160, 240):
            A = E_custom(None, a, b, p2, p4, n=n, fname=fname); B = E_custom(None, b, a, p2, p4, n=n, fname=fname)
            tA, tB = sum(A.values()), sum(B.values())
            pr(f'      {fname}, n = {n:3d}: E(a,b) = {tA.real:+.9e}{tA.imag:+.2e}i   E(b,a) = {tB.real:+.9e}{tB.imag:+.2e}i   '
               f'residual = {resid(A, B, sg):.1e}')
    for n in (40, 240):
        A = F.Eplus(4, a, b, p2, p4, n=n); B = F.Eplus(4, b, a, p2, p4, n=n)
        tA, tB = sum(A.values()), sum(B.values())
        pr(f'      f4 (k = 4), n = {n:3d}: E(a,b) = {tA.real:+.12e}   E(b,a) = {tB.real:+.12e}   residual = {resid(A, B, 1):.1e}')
    pr('  (4d) endpoint (Regge) exponents of |t1^3 t3^3 f G_conn| at the light-ray endpoint t -> 0 (x1^+ -> x0), along t = -i tau:')
    z1 = F.zvec(a); z3 = F.zvec(b); y13 = 2*np.sin((a - b)/2)**2
    taus = np.array([1e-3, 1e-4, 1e-5, 1e-6])
    for lab, path in (('single (t3 = 0.3 - 0.4i fixed)', lambda s: (-1j*s, 0.3 - 0.4j + 0*s)), ('joint t1 = t3', lambda s: (-1j*s, -1j*s))):
        row = []
        for nm in ('f0', 'f1', 2, 4):
            v = []
            for s in taus:
                t1, t3 = path(np.array([s])); t1 = np.atleast_1d(t1); t3 = np.atleast_1d(t3)
                G = F.correlator(F.links(t1, t3, z1, z3, p2, p4, y13=y13))
                f2 = F.fk(2, t1, t3, z1, z3, y13=y13)
                f = {'f0': f2**3, 'f1': f2**2}[nm] if isinstance(nm, str) else F.fk(nm, t1, t3, z1, z3, y13=y13)
                v.append(abs((t1**3*t3**3*f*sum(G.values()))[0]))
            row.append(f'{nm if isinstance(nm, str) else "f" + str(nm)}: tau^{np.polyfit(np.log(taus), np.log(v), 1)[0]:+.2f}')
        pr(f'      {lab}: ' + '   '.join(row) + '   (single endpoint integrable iff exponent > -1)')
    pr('  (4e) no T-shift (operator 1 at the patch point x1 instead of x1^+, so 1 and 3 are timelike; every mixed link gets the constant')
    pr('       phase (-1)^{-1/2}: class weights (1234, 1243, 1324) -> (-1, -1, +1)):')
    for k in (4, 5):
        sg = 1 if k == 4 else -1
        A = F.Eplus(k, a, b, p2, p4, n=80); B = F.Eplus(k, b, a, p2, p4, n=80)
        wgt = {'1234': -1, '1243': -1, '1324': 1}
        sc = max(abs(v) for v in list(A.values()) + list(B.values()))
        pr(f'      k = {k}: residual = {abs(sum(wgt[c]*A[c] for c in CONN) + sg*sum(wgt[c]*B[c] for c in CONN))/sc:.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (5)
def sec_pair():
    from v4c_wu import Kfun
    pr('(5) per-nu light side in the swapped form (nu outside), EE, j\'=0 celestial kernel: I_nu[X] = int d(dphi) A13^{(2+delta)/2} K_nu(dphi) X(dphi)')
    pr('    (v4c_mft part 3 normalisation; both signs of dphi), for X = each cycle class and the total; dphi in [1e-4, pi] (Gauss-Legendre in')
    pr('    log dphi, 64 nodes); the tail below 1e-4 of a single class is ~ 2 x 1e-2 x |coefficient| (dphi^{-1/2} law) and is NOT included.')
    x, wx = np.polynomial.legendre.leggauss(64)
    lo, hi = np.log(1e-4), np.log(np.pi)
    ld = lo + (hi - lo)*(x + 1)/2; wd = wx*(hi - lo)/2
    dph = np.exp(ld); wd = wd*dph
    vals = {c: [] for c in CONN}
    for d in dph:
        r = F.EE_combine(lambda p2, p4: F.Eplus(4, d/2, -d/2, p2, p4, n=80))
        for c in CONN:
            vals[c].append(r[c]/C2)
    vals = {c: np.array(v) for c, v in vals.items()}
    tot = sum(vals.values())
    for nu in (1.0, 3.0):
        dl = 0.5 + 1j*nu
        K = np.array([Kfun(d, nu) for d in dph])
        ker = 2*(2 - 2*np.cos(dph))**((2 + dl)/2)*K*wd
        I = {c: np.sum(ker*v) for c, v in vals.items()}
        It = np.sum(ker*tot)
        Iabs = np.sum(np.abs(ker)*np.abs(vals['1243']))
        pr(f'    nu = {nu}: I[1234] = {abs(I["1234"]):.2e}   I[1243] = {I["1243"].real:+.6e}{I["1243"].imag:+.6e}i   '
           f'I[1324] = {I["1324"].real:+.6e}{I["1324"].imag:+.6e}i   I[total] = {It.real:+.2e}{It.imag:+.2e}i   '
           f'|I[total]|/|I[1243]| = {abs(It)/abs(I["1243"]):.1e}   |I[total]|/int|ker||1243| = {abs(It)/Iabs:.1e}')
    pr(f'    max over the 64 nodes of |total|/|1243| = {np.max(np.abs(tot)/np.abs(vals["1243"])):.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (6)
def sec_gff():
    import og_gff as G
    pr('(6) GENERIC EXPONENTS: O = [phi phi]_{0,2} of a generalised free scalar of dimension Dp (Delta_O = 2Dp + 2, J = 2), connected part,')
    pr('    light-transform weight t^(2Dp+2); a0/a1 = -Dp/(Dp+1) (conformal 2-pt function; Dp = 1/2 is the free T). All links (x^2)^(-Dp - n).')
    rng = np.random.default_rng(11)
    sets = [(rnull(rng), rnull(rng), *rng.uniform(-np.pi, np.pi, 2)) for _ in range(2)]
    for Dp in (0.6, 0.9, 1.3, 2.0):
        a0, spread = G.primary_coef(Dp)
        row = []
        for (p2, p4, a, b) in sets:
            for k in (4, 5):
                sg = 1 if k == 4 else -1
                A = G.Eplus(k, a, b, p2, p4, Dp, a0, n=100); B = G.Eplus(k, b, a, p2, p4, Dp, a0, n=100)
                sc = max(abs(v) for v in list(A.values()) + list(B.values()))
                row.append(f'k={k}: |E|={abs(sum(A.values())):.2e}, res {abs(sum(A.values()) + sg*sum(B.values()))/sc:.0e}')
        pr(f'    Dp = {Dp}: a0/a1 = {a0:+.6f} (2-pt conformality spread {spread:.0e});  ' + ';  '.join(row))
    pr('    d/dDp at Dp = 0.9 (central difference, h = 1e-3): every link acquires log(x_ij^2) terms (the anomalous-dimension structure):')
    p2, p4, a, b = sets[0]
    for k in (4, 5):
        sg = 1 if k == 4 else -1
        out = []
        for h in (+1e-3, -1e-3):
            Dp = 0.9 + h; a0, _ = G.primary_coef(Dp)
            A = sum(G.Eplus(k, a, b, p2, p4, Dp, a0, n=100).values()); B = sum(G.Eplus(k, b, a, p2, p4, Dp, a0, n=100).values())
            out.append((A, B))
        dA = (out[0][0] - out[1][0])/2e-3; dB = (out[0][1] - out[1][1])/2e-3
        pr(f'      k = {k}: dE(a,b)/dDp = {dA.real:+.6e}, dE(b,a)/dDp = {dB.real:+.6e}, |d(E(a,b) {"+" if sg > 0 else "-"} E(b,a))/dDp| / |dE(a,b)/dDp| = '
           f'{abs(dA + sg*dB)/abs(dA):.0e}')
    pr(f'  [{time.time() - T0:.0f}s]')


if __name__ == '__main__':
    secs = sys.argv[1:] or ['val', 'ee', 'tens', 'ctrl', 'pair', 'gff']
    pr('# og_run: O(1/C_T) spinning check of pointwise Lemma A (swapped form, D_LR) on the exact free-scalar connected <TTTT> (og_free.py)')
    for s in secs:
        {'val': sec_val, 'ee': sec_ee, 'tens': sec_tens, 'ctrl': sec_ctrl, 'pair': sec_pair, 'gff': sec_gff}[s]()
    pr(f'# total {time.time() - T0:.0f} s')
    open(os.path.join(CACHE, 'og_check' + ('' if not sys.argv[1:] else '_' + '_'.join(secs)) + '.txt'), 'w').write('\n'.join(L) + '\n')
