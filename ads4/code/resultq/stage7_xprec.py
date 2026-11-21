"""Result Q, Stage 7 (RESULT_Q_STAGE7.md §4): extended-precision realisation of the stored de-shift weights.

Stage 6 §2.4: the stored double-precision weights satisfy W(t) = Phi(t) only for t <~ 10 (C4 value), 14 (C5), because a row residual eps
of the cosh(jt)-normalised Abel system becomes |W - Phi| ~ eps cosh(j t).  Here the SAME discrete systems (same nu grid, same rows, same
gauge = the min-norm solution on the kept singular subspace of the double-precision SVD, tol 1e-11) are solved by ITERATIVE REFINEMENT:
residuals in exact fixed-point integer arithmetic (P bits; the trigonometric matrices are DCT/DST-type, entries cos(pi q/(2 Nt)) from a
table of 4 Nt values), corrections by the double-precision pseudo-inverse.  The fixed point is the min-norm solution of U_k^T (b - M Y) = 0
to P bits; the residual left in the dropped left-singular directions is reported (it is what limits large t).
Then W(t) is evaluated OFF the collocation grid in mpmath (t <= 40) and the realised images are bounded: for Phi of fixed sign the
relative image error at beta is <= sup_{t <= beta} |W/Phi - 1| (positive Abel kernel).
Pieces: C4 value (Phi = 1) and derivative (Phi = -t^2/2) of functionals_q4_deshift.json (EE3, dnu = 0.02, N = 2048);
        C5 (Phi = 1, general complex gauge, 2N unknowns, 2N rows) of functionals_q5_deshift.json (dnu = 0.02, N = 2048).
FINDING (RESULT_Q_STAGE7.md §4): on the stored grid the exact solution reproduces the stored weights' OFF-grid failure (it is truncation
at nu_max = 41, not round-off).  The kept weights therefore use extended grids (env XP_DNU, XP_N; XP_TOL = correction tolerance):
   C4: XP_DNU=0.04 XP_N=3072 (nu_max 123): W = Phi off-grid to <= 4e-15 at t = 30 (min-norm gauge on that grid);
   C5: XP_DNU=0.04 XP_N=2048 XP_TOL=1e-15: stalls (no spectral gap), |W - 1| = 5e-7 / 8e-4 at t = 15.3 / 20.2.
Usage: [XP_DNU=.. XP_N=.. XP_TOL=..] python3 stage7_xprec.py [c4v c4d c5]  -> cache/stage7_xprec_<piece><gridtag>.pkl,
       cache/stage7_xprec.txt (appended);  python3 stage7_xprec.py save  -> functionals_q7_xprec.json (the kept weights above)
"""
import os, sys, io, time, pickle, json
os.environ.setdefault('OMP_NUM_THREADS', '2'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import sympy as sp
import stage4_hj as HJ4
import stage5_hj as HJ5
import stage4_deshift as DS4
import stage5_deshift as DS3

CACHE = os.path.join(HERE, 'cache')
GRID = (os.environ.get('XP_DNU', '0.02'), int(os.environ.get('XP_N', '2048')))
GTAG = ('' if GRID == ('0.02', 2048) else f"_d{GRID[0].replace('.', '')}_N{GRID[1]}") + (f"_tol{os.environ['XP_TOL']}" if 'XP_TOL' in os.environ else '')
P = int(os.environ.get('XP_BITS', '384'))
mp.mp.prec = P + 64
TXT = os.path.join(CACHE, 'stage7_xprec.txt')


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True)
    with open(TXT, 'a') as f:
        f.write(s + "\n")


def to_int(x):
    return int(mp.nint(mp.mpf(x)*mp.mpf(2)**P))


def from_int(i, scale=1):
    return mp.ldexp(mp.mpf(i), -P*scale)


class Trig:
    """C_k(v) = sum_m cos(pi k (2m+1)/(2 Nt)) v_m and S_k(v) (sin), k < Nt, m < N, in exact integers (table scaled 2^P, v scaled 2^P)."""
    def __init__(self, Nt, N):
        self.Nt, self.N, L = Nt, N, 4*Nt
        self.L = L
        q = [mp.pi*mp.mpf(i)/(2*Nt) for i in range(L)]
        self.tc = np.array([to_int(mp.cos(x)) for x in q], dtype=object)
        self.ts = np.array([to_int(mp.sin(x)) for x in q], dtype=object)
        self.odd = 2*np.arange(N, dtype=np.int64) + 1

    def apply(self, vecs, block=256):
        """vecs: list of (kind 'c'|'s', object int vector). Returns list of object int arrays (scale 2^{2P}), length Nt."""
        out = [np.empty(self.Nt, dtype=object) for _ in vecs]
        for k0 in range(0, self.Nt, block):
            ks = np.arange(k0, min(k0 + block, self.Nt), dtype=np.int64)
            idx = np.outer(ks, self.odd) % self.L
            Tc = Ts = None
            for i, (kind, v) in enumerate(vecs):
                if kind == 'c':
                    if Tc is None:
                        Tc = self.tc[idx]
                    out[i][ks] = Tc.dot(v)
                else:
                    if Ts is None:
                        Ts = self.ts[idx]
                    out[i][ks] = Ts.dot(v)
        return out


def mpvec_int(xs):
    return np.array([to_int(x) for x in xs], dtype=object)


def mul_int(a, b):
    """fixed-point product of two int vectors (scale 2^P each) -> scale 2^P (floor)."""
    return np.array([(x*y) >> P for x, y in zip(a, b)], dtype=object)


# ------------------------------------------------------------------------------------------------ the two systems in extended precision
class XC4:
    def __init__(self, dnu_str='0.02', N=2048):
        self.N = N; self.dnu = mp.mpf(dnu_str); self.Nt = N
        nu_s = HJ4.nu_s
        F = {j: sp.lambdify(nu_s, v, 'mpmath') for j, v in HJ4.load()['EE3'].items()}
        nu = [(m + mp.mpf(1)/2)*self.dnu for m in range(N)]
        self.nu = nu
        R = {j: [F[j](mp.mpc(x)) for x in nu] for j in F}
        A1 = [(R[1][m] + R[-1][m]).real for m in range(N)]; A3 = [(R[3][m] + R[-3][m]).real for m in range(N)]
        B1 = [(-1j*(R[1][m] - R[-1][m])).real for m in range(N)]; B3 = [(-1j*(R[3][m] - R[-3][m])).real for m in range(N)]
        n = [sum(abs(R[j][m]) for j in R) for m in range(N)]
        self.R, self.A1, self.A3, self.B1, self.B3, self.n = R, A1, A3, B1, B3, n
        self.ia1, self.ia3 = mpvec_int([A1[m]/n[m] for m in range(N)]), mpvec_int([A3[m]/n[m] for m in range(N)])
        self.ib1, self.ib3 = mpvec_int([B1[m]/n[m] for m in range(N)]), mpvec_int([B3[m]/n[m] for m in range(N)])
        self.t = [mp.pi*k/(N*self.dnu) for k in range(N)]
        self.ch3 = [mp.cosh(3*t) for t in self.t]
        self.r1 = [mp.cosh(t)/c for t, c in zip(self.t, self.ch3)]
        self.s1 = [mp.sinh(t)/c for t, c in zip(self.t, self.ch3)]
        self.th = [mp.tanh(3*t) for t in self.t]
        self.trig = Trig(N, N)
        self.D = DS4.Deshift4('EE3', float(self.dnu), N)      # double-precision system (same grid to 1e-16) for the corrections

    def matvec(self, Yint):
        v = self.trig.apply([('c', mul_int(self.ia3, Yint)), ('c', mul_int(self.ia1, Yint)),
                             ('s', mul_int(self.ib3, Yint)), ('s', mul_int(self.ib1, Yint))])
        out = []
        for k in range(self.Nt):
            c3, c1, s3, s1 = (from_int(x[k], 2) for x in v)
            out.append(self.dnu*(c3 + self.r1[k]*c1 - self.th[k]*s3 - self.s1[k]*s1))
        return out

    def rhs(self, Phi):
        return [Phi(t)/c for t, c in zip(self.t, self.ch3)]

    def W(self, chi, t):
        t = mp.mpf(t); ch, sh, ch3, sh3 = mp.cosh(t), mp.sinh(t), mp.cosh(3*t), mp.sinh(3*t)
        tot = mp.mpf(0)
        for m in range(self.N):
            c, s = mp.cos(self.nu[m]*t), mp.sin(self.nu[m]*t)
            tot += chi[m]*(ch*self.A1[m]*c - sh*self.B1[m]*s + ch3*self.A3[m]*c - sh3*self.B3[m]*s)
        return tot*self.dnu


class XC5:
    def __init__(self, dnu_str='0.02', N=2048):
        self.N = N; self.dnu = mp.mpf(dnu_str); self.Nt = 2*N
        nu_s = HJ5.nu_s
        F = {j: sp.lambdify(nu_s, v, 'mpmath') for j, v in HJ5.load()['C5_pp'].items()}
        nu = [(m + mp.mpf(1)/2)*self.dnu for m in range(N)]
        self.nu = nu
        R = {j: [F[j](mp.mpc(x)) for x in nu] for j in F}
        Ns = [R[-2][m] + R[0][m] + R[2][m] for m in range(N)]
        hh = {j: [R[j][m]/Ns[m] for m in range(N)] for j in R}
        self.hh = hh
        a0 = [hh[0][m].real for m in range(N)]
        a2 = [(hh[2][m] + mp.conj(hh[-2][m]))/2 for m in range(N)]
        n = [abs(a0[m]) + 2*abs(a2[m]) for m in range(N)]
        self.n = n
        Sp = [hh[2][m] + hh[-2][m] for m in range(N)]; Dm = [hh[2][m] - hh[-2][m] for m in range(N)]
        self.c0R = mpvec_int([hh[0][m].real/n[m] for m in range(N)]); self.c0I = mpvec_int([-hh[0][m].imag/n[m] for m in range(N)])
        self.cAR = mpvec_int([Sp[m].real/n[m] for m in range(N)]); self.cAI = mpvec_int([-Sp[m].imag/n[m] for m in range(N)])
        self.cBR = mpvec_int([Dm[m].imag/n[m] for m in range(N)]); self.cBI = mpvec_int([Dm[m].real/n[m] for m in range(N)])
        Nt = self.Nt
        self.t = [mp.pi*k/(Nt*self.dnu) for k in range(Nt)]
        self.ch2 = [mp.cosh(2*t) for t in self.t]
        self.r0 = [1/c for c in self.ch2]
        self.th = [mp.tanh(2*t) for t in self.t]
        self.trig = Trig(Nt, N)
        self.D = DS3.Deshift3S('C5_pp', float(self.dnu), N, general=True)

    def matvec(self, Yint):
        YR, YI = Yint[:self.N], Yint[self.N:]
        add = lambda x, y: np.array([a + b for a, b in zip(x, y)], dtype=object)
        v0 = add(mul_int(self.c0R, YR), mul_int(self.c0I, YI))
        vA = add(mul_int(self.cAR, YR), mul_int(self.cAI, YI))
        vB = add(mul_int(self.cBR, YR), mul_int(self.cBI, YI))
        c0, cA, sB = self.trig.apply([('c', v0), ('c', vA), ('s', vB)])
        return [self.dnu*(self.r0[k]*from_int(c0[k], 2) + from_int(cA[k], 2) - self.th[k]*from_int(sB[k], 2)) for k in range(self.Nt)]

    def rhs(self, Phi):
        return [Phi(t)/c for t, c in zip(self.t, self.ch2)]

    def W(self, c, t):
        t = mp.mpf(t); chs = {j: mp.cosh(j*t) for j in (0, 2, -2)}; shs = {j: mp.sinh(j*t) for j in (0, 2, -2)}
        tot = mp.mpf(0)
        for m in range(self.N):
            cc, ss = mp.cos(self.nu[m]*t), mp.sin(self.nu[m]*t)
            acc = 0
            for j in (0, 2, -2):
                acc += self.hh[j][m]*(cc*chs[j] + 1j*ss*shs[j])
            tot += (c[m]*acc).real
        return tot*self.dnu


# ------------------------------------------------------------------------------------------------ iterative refinement
def refine(X, Y0, Phi, tol=float(os.environ.get('XP_TOL', '1e-11')), maxit=40, tag=''):
    U, s, Vh = X.D.svd()
    keep = s > tol*s[0]
    Uk, sk, Vk = U[:, keep], s[keep], Vh[keep]
    Udrop = U[:, ~keep]
    b = X.rhs(Phi)
    tt = np.array([float(t) for t in X.t])
    big = [float(1/c) if c != 0 else 0 for c in (X.ch3 if hasattr(X, 'ch3') else X.ch2)]
    Yint = np.array([int(round(y*2.0**60)) << (P - 60) for y in Y0], dtype=object)
    hist = []
    for it in range(maxit):
        t0 = time.time()
        r = [bb - mm for bb, mm in zip(b, X.matvec(Yint))]
        rf = np.array([float(x) for x in r])
        # |W - Phi| on the collocation grid = |r_k| cosh(j t_k)
        amp = np.array([float(abs(x)) for x in r])/np.array(big)
        rel = []
        for T in (10, 15, 20, 25, 30, 40):
            sel = tt <= T
            rel.append(float(np.max(amp[sel])))
        rdrop = float(np.linalg.norm(Udrop.T @ rf)) if Udrop.shape[1] else 0.0
        hist.append((it, float(np.max(np.abs(rf))), rel, rdrop))
        pr(f"  [{tag}] it {it:2d}: max|r| = {np.max(np.abs(rf)):.2e}; max|W - Phi| on grid t <= 10/15/20/25/30/40: "
           + " ".join(f"{x:.1e}" for x in rel) + f"; |U_drop^T r| = {rdrop:.1e}  [{time.time() - t0:.0f}s]")
        if it > 0 and (np.max(np.abs(rf)) < 2.0**(-P + 40) or (len(hist) > 3 and hist[-1][1] > 0.5*hist[-3][1])):
            break
        dY = Vk.T @ ((Uk.T @ rf)/sk)
        # rescale to keep relative precision of the tiny corrections
        sc = max(np.max(np.abs(dY)), 1e-300)
        e = int(np.floor(np.log2(sc))) - 52
        dYi = np.array([int(round(x*2.0**(-e))) for x in dY], dtype=object)
        Yint = np.array([y + ((d << (P + e)) if P + e >= 0 else (d >> -(P + e))) for y, d in zip(Yint, dYi)], dtype=object)
    return Yint, hist


def offgrid_W(X, chi, Phi, ts):
    out = []
    for t in ts:
        W = X.W(chi, t); ph = Phi(mp.mpf(t))
        out.append((t, float(W), float(abs(W/ph - 1)) if ph != 0 else float(abs(W))))
    return out


def run(piece):
    t00 = time.time()
    dnu_s, N = GRID
    stored = (dnu_s == '0.02' and N == 2048)
    if piece in ('c4v', 'c4d'):
        X = XC4(dnu_s, N)
        Phi = (lambda t: mp.mpf(1) + 0*t) if piece == 'c4v' else (lambda t: -t*t/2)
        st = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']
        chi0 = np.array(st['chi_value' if piece == 'c4v' else 'chi_deriv'])
        if stored:
            Y0 = chi0*X.D.n
        else:                                   # new grid: double min-norm solve there; compare with the stored weight on the common nu
            Y0 = X.D.solve(lambda t: np.ones_like(t) if piece == 'c4v' else -0.5*t**2)[1]
            chi0 = np.interp(X.D.nu, (np.arange(len(chi0)) + 0.5)*0.02, chi0, right=0.0)
    else:
        X = XC5(dnu_s, N)
        st = json.load(open(os.path.join(HERE, 'functionals_q5_deshift.json')))
        c0 = np.array(st['c_re']) + 1j*np.array(st['c_im'])
        Phi = lambda t: mp.mpf(1) + 0*t
        if stored:
            Y0 = np.concatenate([(c0*X.D.n).real, (c0*X.D.n).imag])
        else:
            Y0 = X.D.solve(lambda t: np.ones_like(t))[1]
            nu0 = (np.arange(len(c0)) + 0.5)*0.02
            c0 = np.interp(X.D.nu, nu0, c0.real, right=0.0) + 1j*np.interp(X.D.nu, nu0, c0.imag, right=0.0)
    pr(f"\n## {piece} (grid dnu = {dnu_s}, N = {N}, nu_max = {float(X.dnu)*N:.0f}): built extended-precision system (P = {P} bits) [{time.time() - t00:.0f}s]; double residual of the stored weight: "
       f"{np.linalg.norm(X.D.M @ Y0 - np.array([float(x) for x in X.rhs(Phi)]))/np.linalg.norm([float(x) for x in X.rhs(Phi)]):.1e}")
    Yint, hist = refine(X, Y0, Phi, tag=piece)
    Y = [from_int(y) for y in Yint]
    if piece in ('c4v', 'c4d'):
        chi = [Y[m]/X.n[m] for m in range(X.N)]
        dchi = max(abs(float(chi[m]) - chi0[m]) for m in range(X.N))/np.max(np.abs(chi0))
        cvals = chi
    else:
        cvals = [(Y[m] + 1j*Y[X.N + m])/X.n[m] for m in range(X.N)]
        dchi = max(abs(complex(cvals[m]) - c0[m]) for m in range(X.N))/np.max(np.abs(c0))
    nrm = sum(abs(complex(w))*float(n) for w, n in zip(cvals, X.n))*float(X.dnu)
    pr(f"  [{piece}] int |weight| n dnu = {nrm:.4f}; max|weight_xp - weight_stored|/max|weight| = {dchi:.1e}  ({'same grid: same functional to this accuracy' if stored else 'new grid: gauge change vs the stored weight (interpolated)'})")
    ts = [0.5, 1, 2, 3.3, 5, 7.7, 10.1, 12.6, 15.3, 17.9, 20.2, 22.7, 25.1, 27.4, 30.0, 32.5, 35.2, 40.0]
    og = offgrid_W(X, cvals, Phi, ts)
    pr(f"  [{piece}] OFF-grid |W/Phi - 1| (mpmath direct quadrature): " + "; ".join(f"t={t}: {e:.1e}" for t, _, e in og))
    # the stored double weights, same evaluator (for the record)
    if piece in ('c4v', 'c4d'):
        og0 = offgrid_W(X, [mp.mpf(x) for x in chi0], Phi, [8, 10, 12, 14, 16])
    else:
        og0 = offgrid_W(X, [mp.mpc(complex(x)) for x in c0], Phi, [8, 10, 12, 14, 16])
    pr(f"  [{piece}] stored double weights, same evaluator: " + "; ".join(f"t={t}: {e:.1e}" for t, _, e in og0))
    # light-side moments M_n (shifted) as a check that the functional is unchanged
    pickle.dump(dict(piece=piece, P=P, Y=[mp.nstr(y, 70) for y in Y], weights=[mp.nstr(w, 70) for w in cvals], hist=hist, offgrid=og,
                     offgrid_stored=og0, dchi=dchi), open(os.path.join(CACHE, f'stage7_xprec_{piece}{GTAG}.pkl'), 'wb'))
    pr(f"  [{piece}] done [{time.time() - t00:.0f}s]")


def save():
    """Kept weights: C4 value / derivative from the extended grid dnu = 1/25, N = 3072 (nu_max = 122.88; min-norm tol 1e-11 there, refined
    to 384 bits); C5 from dnu = 1/25, N = 2048 (nu_max = 81.92), correction tol 1e-15 (the best C5 run; realised only to t ~ 18-20)."""
    out = dict(description=(
        "Result Q Stage 7: extended-precision de-shift weights (iterative refinement with %d-bit fixed-point residuals; corrections from the "
        "double-precision pseudo-inverse).  Grids nu_k = (k+1/2) dnu with dnu given EXACTLY as a fraction.  C4_value / C4_deriv: Stage-4 "
        "4-shift EE3 system (R_j of functionals_q4_deshift.json) on the extended grid dnu = 1/25, N = 3072 (nu_max = 122.88), min-norm "
        "(tol 1e-11) gauge ON THAT GRID: a gauge change of the stored dnu = 1/50, N = 2048 weights (max relative difference 'dweight', "
        "computed against the stored weight interpolated to the new grid).  C5: Stage-5 3-shift general complex system (c = c_re + i c_im, "
        "functional Re[c Psi_5/Nhat5]) on dnu = 1/25, N = 2048, correction tolerance 1e-15.  'offgrid' lists (t, |W/Phi - 1|) by direct "
        "mpmath quadrature OFF the collocation grid.  C4 pieces: W = Phi to <= 4e-15 for t <= 30.  C5: only to ~1e-6 at t = 15, 8e-4 at "
        "t = 20 (no spectral gap; see RESULT_Q_STAGE7.md).  Weights are decimal strings (64 significant digits); evaluate W(t) for "
        "t > 10 with >= 50 digits." % P))
    files = {'c4v': 'stage7_xprec_c4v_d004_N3072.pkl', 'c4d': 'stage7_xprec_c4d_d004_N3072.pkl', 'c5': 'stage7_xprec_c5_d004_N2048_tol1e-15.pkl'}
    for piece, f in files.items():
        fn = os.path.join(CACHE, f)
        if not os.path.exists(fn):
            continue
        d = pickle.load(open(fn, 'rb'))
        w = [mp.mpmathify(x) for x in d['weights']]
        if piece == 'c5':
            ent = dict(dnu='1/25', N=len(w), c_re=[mp.nstr(mp.re(x), 64) for x in w], c_im=[mp.nstr(mp.im(x), 64) for x in w])
        else:
            ent = dict(dnu='1/25', N=len(w), chi=[mp.nstr(x, 64) for x in w])
        ent['dweight'] = d['dchi']; ent['offgrid'] = [(t, e) for t, _, e in d['offgrid']]; ent['source_cache'] = f
        out[{'c4v': 'C4_value', 'c4d': 'C4_deriv', 'c5': 'C5'}[piece]] = ent
    json.dump(out, open(os.path.join(HERE, 'functionals_q7_xprec.json'), 'w'))
    print("wrote functionals_q7_xprec.json")


if __name__ == '__main__':
    args = sys.argv[1:] or ['c4d', 'c4v', 'c5']
    if args == ['save']:
        save()
    else:
        for a in args:
            run(a)
