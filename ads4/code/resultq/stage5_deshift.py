"""Result Q, Stage 5 (RESULT_Q_STAGE5.md §2): the 3-shift de-shift (Abel inversion) of the raw SPINNING k = 5 functional Psi_{5,nu}[o_H q]
(j' = 1), giving the pure-moment forward functional C5ds := int_0^oo chi(nu) Re[Psi_{5,nu}/Nhat5(nu)] dnu with heavy image
sigma_h P_{-1/2}(cosh beta) (= "C^_5(0)" of the Stage-12/3 model).

Exact image (stage5_hj): Psi_{5,nu} -> N(nu) sum_{j=0,+-2} R_j(nu) P_{l+j}, helicity class ++ ; +- = -(++) exactly (sigma_h).
New feature vs k = 4: the coefficients are real-analytic (conj R_j(nu)/R_0(nu) = R_j(-nu)/R_0(-nu), checked) but NOT shadow-symmetric
(R_{-j}(nu) != c(nu) R_j(-nu)), so no nu-only phase makes them "real" (h_0 real, h_-2 = conj h_2).  A real sum rule is Re[c(nu) Psi_nu]
for any complex c; its image is sum_j a_j P_{l+j} with a_j = (c h_j + conj(c h_{-j}))/2 (a_0 real, a_-2 = conj a_2).
CANONICAL (main) choice: c = chi(nu)/Nhat5(nu), chi REAL, Nhat5 := sum_j h_j (forward normalisation of the raw functional; complex; no
real zero on nu > 0, checked) -> a_j = (hh_j + conj hh_{-j})/2, hh_j := R_j/sum_i R_i (independent of every dropped nu-only factor);
the forward normalisation of Re[Psi/Nhat5] is sum_j a_j = 1 identically.  Mehler-Dirichlet (complex index):
    W[chi](t) = C[chi a_0] + cosh 2t C[chi A2] - sinh 2t S[chi B2],  A2 = Re(hh_2 + hh_-2), B2 = Im(hh_2 - hh_-2)   (all real)
and W = Phi on t >= 0 <=> image = int psi P (Abel injective).  Target Phi = 1 (psi = delta(nu - 0+), image P_{-1/2}).
chi-normalisation (EE / stage-4 way): Y = chi n, n = |a_0| + |a_2| + |a_-2| (hh_{+-2} ~ 1/nu at nu -> 0: bounded operator in Y);
rows divided by cosh 2t.  Midpoint nu-grid, paired t-grid t_k = k pi/(N dnu), dense SVD, min-norm solution.
GENERAL choice (reported for the null-space count): c = (chi_R + i chi_I)/Nhat5, two real unknowns per nu.
Scalar control: the exact rank-2 scalar ('sc2', = stage-4 proxy b = 10; shadow-symmetric, a_j real-structured).
Checks: (V1) residual; (V2) INDEPENDENT image by direct nu-quadrature with mpmath complex-degree Legendre functions; (V3) light-side
moments Re sum_j int chi hh_j (nu - ij)^{2n} = (1, 0, 0) (n = 0, 1, 2); (R) refinement; null space; weight size int|chi|.
Usage: python3 stage5_deshift.py -> cache/stage5_deshift.txt, cache/stage5_deshift_<key>.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import stage5_hj as S5

CACHE = os.path.join(HERE, 'cache')


class Deshift3S:
    def __init__(self, key, dnu, N, general=False, os_=None):
        F = S5.numfuncs(key)
        self.key, self.dnu, self.N, self.general = key, dnu, N, general
        nu = (np.arange(N) + 0.5)*dnu; self.nu = nu
        os_ = (2 if general else 1) if os_ is None else os_         # t-oversampling (square system for the 2-unknown gauge)
        Nt = int(os_*N)
        t = np.arange(Nt)*np.pi/(Nt*dnu); self.t = t
        R = {j: F[j](nu.astype(complex)) for j in F}
        Ns = R[-2] + R[0] + R[2]
        self.Nabs_rel = np.abs(Ns)/sum(np.abs(R[j]) for j in R)
        hh = {j: R[j]/Ns for j in R}
        self.hh = hh
        a0 = hh[0].real
        a2 = 0.5*(hh[2] + np.conj(hh[-2])); am2 = np.conj(a2)
        self.a = {0: a0 + 0j, 2: a2, -2: am2}
        self.n = np.abs(a0) + 2*np.abs(a2)
        self.A2 = (hh[2] + hh[-2]).real; self.B2 = (hh[2] - hh[-2]).imag
        # general (complex-c) columns: a_0 = Re(c hh0), A = Re(c (hh2 + hh-2)), B = Im(c (hh2 - hh-2)), c = chi_R + i chi_I
        Cm = dnu*np.cos(np.outer(t, nu)); Sm = dnu*np.sin(np.outer(t, nu))
        r0 = 1/np.cosh(np.minimum(2*t, 700)); th = np.tanh(2*t)
        self.ch2 = np.cosh(np.minimum(2*t, 700))
        def cols(z0, zA, zB):           # coefficient functions (real) of one real unknown
            return r0[:, None]*Cm*(z0/self.n)[None, :] + Cm*(zA/self.n)[None, :] - th[:, None]*Sm*(zB/self.n)[None, :]
        Sp = hh[2] + hh[-2]; Dm = hh[2] - hh[-2]
        MR = cols(hh[0].real, Sp.real, Dm.imag)
        if general:
            MI = cols(-hh[0].imag, -Sp.imag, Dm.real)
            self.M = np.hstack([MR, MI])
        else:
            self.M = MR
        del Cm, Sm

    def svd(self):
        if not hasattr(self, '_svd'):
            self._svd = np.linalg.svd(self.M, full_matrices=False)
        return self._svd

    def solve(self, Phi, tol=1e-11):
        U, s, Vh = self.svd()
        b = Phi(self.t)/self.ch2
        keep = s > tol*s[0]
        Y = Vh[keep].T @ ((U[:, keep].T @ b)/s[keep])
        res = np.linalg.norm(self.M @ Y - b)/np.linalg.norm(b)
        if self.general:
            c = (Y[:self.N] + 1j*Y[self.N:])/self.n
        else:
            c = Y/self.n + 0j
        return c, Y, res, int(np.sum(~keep)) + max(0, self.M.shape[1] - len(s))

    def W(self, c, t):
        """W(t) = Re sum_j int c hh_j cos((nu - ij) t) dnu by direct quadrature."""
        out = []
        for tt in np.atleast_1d(t):
            cc, ss = np.cos(self.nu*tt), np.sin(self.nu*tt)
            tot = 0
            for j in (0, 2, -2):
                tot = tot + c*self.hh[j]*(cc*np.cosh(j*tt) + 1j*ss*np.sinh(j*tt))
            out.append(self.dnu*np.sum(tot).real)
        return np.array(out)

    def moments(self, c, nmax=2):
        return [self.dnu*np.sum(sum(c*self.hh[j]*(self.nu - 1j*j)**(2*n) for j in (0, 2, -2))).real for n in range(nmax + 1)]


def image_direct(S, c, beta, stride=1):
    """Re sum_j int c hh_j P_{-1/2 + i(nu - ij)}(cosh beta) dnu, mpmath Legendre (independent of the Mehler / DCT machinery)."""
    x = mp.cosh(beta); tot = mp.mpf(0)
    for k in range(0, S.N, stride):
        if abs(c[k]) < 1e-300:
            continue
        nu = S.nu[k]
        acc = 0
        for j in (0, 2, -2):
            acc += complex(S.hh[j][k])*complex(mp.legenp(-0.5 + 1j*(nu - 1j*j), 0, x, type=3))
        tot += (c[k]*acc).real
    return float(tot)*S.dnu*stride


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t00 = time.time()
    one = lambda t: np.ones_like(t)
    TT = np.array([0, 0.5, 1, 2, 4, 6, 8.])
    keys = sys.argv[1:] or ['C5_pp', 'sc2']
    for key in keys:
        pr(f"\n## {key}")
        pr("  (A) restricted gauge: chi REAL on Re[Psi_5/Nhat5] (forward normalisation 1):")
        for (dnu, N) in ((0.04, 1024), (0.02, 2048)):
            t0 = time.time()
            S = Deshift3S(key, dnu, N)
            U, s, Vh = S.svd()
            c, Y, res, nnull = S.solve(one)
            pr(f"    dnu = {dnu}, N = {N}: smallest sv/s_max " + " ".join(f"{v/s[0]:.1e}" for v in s[::-1][:4]) + f"; residual {res:.1e};"
               f" int|chi| = {S.dnu*np.sum(np.abs(c)):.4e}; W(t = 0,.5,1,2,4,6,8) = " + " ".join(f"{v:+.6f}" for v in S.W(c, TT)) +
               f"  [{time.time() - t0:.0f}s]")
        pr("  (B) general gauge: complex c = chi_R + i chi_I on Psi_5/Nhat5 (real sum rule Re[c Psi/Nhat5]); t-grid oversampled x2; min-norm:")
        G = {}
        for (dnu, N) in ((0.04, 1024), (0.02, 2048), (0.02, 3072)):
            t0 = time.time()
            S = Deshift3S(key, dnu, N, general=True)
            U, s, Vh = S.svd()
            c, Y, res, nnull = S.solve(one)
            supp = S.nu[np.abs(c) > 1e-6*np.abs(c).max()].max()
            pr(f"    dnu = {dnu}, N = {N} (nu_max {N*dnu:.0f}): system {S.M.shape}; rank-deficiency (sv < 1e-11 s_max) {nnull}"
               f" (= N - {N - nnull}); gap: sv/s_max at rank edge {s[N*2 - nnull - 1]/s[0]:.1e} | {s[min(2*N - 1, 2*N - nnull)]/s[0]:.1e};"
               f" residual {res:.1e}; int c Nhat-normalised = {S.dnu*np.sum(c.real):+.10f}; int|c| = {S.dnu*np.sum(np.abs(c)):.5f};"
               f" support (|c| > 1e-6 max) nu <= {supp:.1f}; W(t = 0,.5,1,2,4,6,8) = " + " ".join(f"{v:+.10f}" for v in S.W(c, TT)) +
               f"  [{time.time() - t0:.0f}s]")
            G[(dnu, N)] = (S, c)
            del U, s, Vh
        S1, c1 = G[(0.02, 2048)]; S2, c2 = G[(0.02, 3072)]; S0, c0 = G[(0.04, 1024)]
        ci = np.interp(S0.nu, S1.nu, c1.real) + 1j*np.interp(S0.nu, S1.nu, c1.imag)
        pr(f"    (R) nu_max 41 -> 61: max|dc|/max|c| = {np.max(np.abs(c1 - c2[:2048]))/np.max(np.abs(c2)):.1e};  dnu .04 -> .02: "
           f"{np.max(np.abs(c0 - ci))/np.max(np.abs(c1)):.1e}")
        S, c = G[(0.02, 2048)]
        mom = S.moments(c)
        pr(f"    (V3) light-side moments Re sum_j int c hh_j (nu - ij)^(2n), n = 0, 1, 2: " + ", ".join(f"{v:+.3e}" for v in mom) +
           "  (target 1, 0, 0)")
        pr("    profile: nu, |Nhat5|/sum|R|, Re c, Im c, a0 = Re hh0, Re(hh2 + hh-2), Im(hh2 - hh-2)")
        for v in (0.01, 0.3, 0.7, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 40):
            k = min(S.N - 1, int(v/S.dnu))
            pr(f"       {S.nu[k]:7.3f}  {S.Nabs_rel[k]:.3e}  {c[k].real:+.4e}  {c[k].imag:+.4e}  {S.hh[0][k].real:+.4e}  {S.A2[k]:+.4e}  {S.B2[k]:+.4e}")
        pr("    (V2) direct image Re sum_j int c hh_j P_{-1/2+i(nu-ij)}(cosh beta) dnu (mpmath, full grid) vs P_{-1/2}(cosh beta):")
        for beta in (0.0, 0.5, 1.0, 2.0, 4.0, 6.0):
            img = image_direct(S, c, beta, stride=1)
            tgt = float(mp.legenp(-0.5, 0, mp.cosh(beta), type=3).real)
            pr(f"       beta = {beta:3.1f}: image {img:+.12e}  target {tgt:+.12e}  rel.diff {abs(img - tgt)/abs(tgt):.1e}")
        pickle.dump(dict(nu=S.nu, c=c, hh={j: S.hh[j] for j in S.hh}, n=S.n, dnu=S.dnu, key=key),
                    open(os.path.join(CACHE, f'stage5_deshift_{key}.pkl'), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'stage5_deshift.txt'), 'w').write(buf.getvalue())
