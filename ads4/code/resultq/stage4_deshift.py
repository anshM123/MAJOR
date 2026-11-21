"""Result Q, Stage 4 (RESULT_Q_STAGE4.md §2): the 4-shift de-shift (Abel inversion) of the raw k = 4 functional Psi_{4,nu}[EE], giving the
pure-moment FORWARD functionals  C4ds := int_0^oo chi_0(nu) Psi_{4,nu} dnu  (heavy image P_{-1/2}(cosh beta), = "C^_4(0)") and
dC4ds := int chi_2 Psi_{4,nu} dnu (heavy image d/d(nu^2) P_{-1/2+i nu}(cosh beta)|_0 = -D(beta), = "d_{nu^2} C^_4(0)").

Leading-Regge image of Psi_{4,nu} (stage4_hj, exact): sum_{j = +-1, +-3} R_j(nu) P(nu - i j),  P(mu) := P_{-1/2 + i mu}(cosh beta),
(common positive nu-only factor N_com absorbed into chi).  With A_j := R_j + R_{-j}, B_j := -i (R_j - R_{-j}) (both REAL on nu > 0, checked)
and Mehler-Dirichlet for complex index, cos((nu - i j) t) = cos(nu t) cosh(j t) + i sin(nu t) sinh(j t):
    image(beta) = (sqrt2/pi) int_0^beta W(t) (cosh beta - cosh t)^{-1/2} dt,
    W[chi](t) = sum_{j = 1, 3} [ cosh(j t) C[chi A_j](t) - sinh(j t) S[chi B_j](t) ]        (C, S: cos/sin transforms on nu > 0).
The de-shift is the (generalised) Abel equation  W[chi](t) = Phi(t) for all t >= 0, Phi = C[psi]:  psi = delta(nu - 0+)  <->  Phi = 1
(value), psi = d_{nu^2} delta  <->  Phi = -t^2/2 (derivative).  Abel transform injective => image = int psi P exactly.
chi-normalisation as u3_chi_EE (EE_FUNCTIONAL §3): Y := chi n, n := sum_j |R_j| (> 0), a_j = A_j/n, b_j = B_j/n bounded; the forward
normalisation Nhat(nu) := sum_j R_j = A_1 + A_3 (which has REAL zeros) is never divided by.  Rows divided by cosh 3t:
    At[Y] := C[a3 Y] - tanh(3t) S[b3 Y] + (cosh t/cosh 3t) C[a1 Y] - (sinh t/cosh 3t) S[b1 Y] = Phi/cosh 3t.
Discretisation as v4_chi: midpoint nu-grid, paired t-grid t_k = k pi/(N dnu); dense matrix, SVD (null space), min-norm solution.
Checks: (V1) residual; (V2) INDEPENDENT image check by direct nu-quadrature with complex-degree Legendre functions (mpmath) at several beta;
(V3) light-side normalisation int chi Nhat = Phi(0) (R^4: -C4ds|_{R^4} = 2 g4 * Phi(0)); (R) grid refinement; amplification
||chi Nhat||_1 and int |Y| (size of the weight); null-space dimension and support.
Usage: python3 stage4_deshift.py [EE3|sc3 ...]  -> cache/stage4_deshift.txt, cache/stage4_deshift_<key>.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import stage4_hj as HJ

CACHE = os.path.join(HERE, 'cache')


class Deshift4:
    def __init__(self, key, dnu, N):
        F = HJ.numfuncs(key)
        self.key, self.dnu, self.N = key, dnu, N
        nu = (np.arange(N) + 0.5)*dnu; self.nu = nu
        t = np.arange(N)*np.pi/(N*dnu); self.t = t
        Rj = {j: F[j](nu.astype(complex)) for j in F}
        self.Rj = Rj
        A1 = (Rj[1] + Rj[-1]); A3 = (Rj[3] + Rj[-3]); B1 = -1j*(Rj[1] - Rj[-1]); B3 = -1j*(Rj[3] - Rj[-3])
        self.imag_check = max(np.max(np.abs(x.imag)/np.abs(x)) for x in (A1, A3, B1, B3))
        A1, A3, B1, B3 = A1.real, A3.real, B1.real, B3.real
        self.n = sum(np.abs(Rj[j]) for j in Rj)
        self.Nhat = A1 + A3
        self.a1, self.a3, self.b1, self.b3 = A1/self.n, A3/self.n, B1/self.n, B3/self.n
        Cm = dnu*np.cos(np.outer(t, nu)); Sm = dnu*np.sin(np.outer(t, nu))
        r1 = np.exp(t - 3*t)*(1 + np.exp(-2*t))/(1 + np.exp(-6*t))        # cosh t/cosh 3t
        s1 = np.exp(t - 3*t)*(1 - np.exp(-2*t))/(1 + np.exp(-6*t))        # sinh t/cosh 3t
        th = np.tanh(3*t)
        self.M = (Cm*self.a3[None, :] - th[:, None]*Sm*self.b3[None, :] + r1[:, None]*Cm*self.a1[None, :]
                  - s1[:, None]*Sm*self.b1[None, :])
        self.ch3 = np.cosh(np.minimum(3*t, 700))
        del Cm, Sm

    def svd(self):
        if not hasattr(self, '_svd'):
            self._svd = np.linalg.svd(self.M)
        return self._svd

    def solve(self, Phi, tol=1e-11):
        """min-norm least-squares Y (singular values < tol*s_max dropped). Returns chi, Y, residual."""
        U, s, Vh = self.svd()
        b = Phi(self.t)/self.ch3
        keep = s > tol*s[0]
        Y = Vh[keep].T @ ((U[:, keep].T @ b)/s[keep])
        res = np.linalg.norm(self.M @ Y - b)/np.linalg.norm(b)
        return Y/self.n, Y, res

    def W(self, chi, t):
        """W[chi](t) by direct quadrature (any t)."""
        nu = self.nu; d = self.dnu
        A1 = self.a1*self.n; A3 = self.a3*self.n; B1 = self.b1*self.n; B3 = self.b3*self.n
        out = []
        for tt in np.atleast_1d(t):
            c, s = np.cos(nu*tt), np.sin(nu*tt)
            out.append(d*np.sum(chi*(np.cosh(tt)*A1*c - np.sinh(tt)*B1*s + np.cosh(3*tt)*A3*c - np.sinh(3*tt)*B3*s)))
        return np.array(out)


def image_direct(S, chi, beta, numax=None, stride=1):
    """sum_j int chi(nu) R_j(nu) P_{-1/2 + i(nu - i j)}(cosh beta) dnu by direct quadrature, mpmath Legendre functions (independent of
    the Mehler/DCT machinery)."""
    sel = np.arange(0, S.N, stride)
    if numax is not None:
        sel = sel[S.nu[sel] <= numax]
    x = mp.cosh(beta); tot = mp.mpc(0)
    for k in sel:
        nu = S.nu[k]
        if abs(chi[k]) < 1e-300:
            continue
        acc = 0
        for j in (1, -1, 3, -3):
            acc += complex(S.Rj[j][k])*complex(mp.legenp(-0.5 + 1j*(nu - 1j*j), 0, x, type=3))
        tot += chi[k]*acc
    return complex(tot)*S.dnu*stride


def Dbeta(beta):
    """D(beta) = -d_{nu^2} P_{-1/2+i nu}(cosh beta)|_0 = (sqrt2/2pi) int_0^beta t^2 (cosh beta - cosh t)^{-1/2} dt."""
    f = lambda t: t**2/mp.sqrt(mp.cosh(beta) - mp.cosh(t))
    return float(mp.sqrt(2)/(2*mp.pi)*mp.quad(f, [0, beta]))


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    keys = sys.argv[1:] or ['EE3', 'sc3']
    t00 = time.time()
    TGT = {'value': lambda t: np.ones_like(t), 'deriv': lambda t: -0.5*t**2}
    for key in keys:
        pr(f"\n## {key}: 4-shift de-shift of Psi_4 (R_j from stage4_hj)")
        res_by_grid = {}
        for (dnu, N) in ((0.04, 1024), (0.02, 2048), (0.02, 3072)):
            t0 = time.time()
            S = Deshift4(key, dnu, N)
            U, s, Vh = S.svd()
            nsmall = int(np.sum(s < 1e-9*s[0]))
            pr(f"  grid dnu = {dnu}, N = {N} (nu_max = {N*dnu:.0f}): Im/|.| of A_j, B_j = {S.imag_check:.1e}; "
               f"smallest sv/s_max: " + " ".join(f"{v/s[0]:.1e}" for v in s[::-1][:6]) + f"; # < 1e-9: {nsmall}  [{time.time() - t0:.0f}s]")
            out = {}
            for nm, Phi in TGT.items():
                chi, Y, res = S.solve(Phi)
                w0 = S.dnu*np.sum(chi*S.Nhat)
                l1 = S.dnu*np.sum(np.abs(chi*S.Nhat)); l1Y = S.dnu*np.sum(np.abs(Y))
                tt = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
                Wd = S.W(chi, tt)
                supp = S.nu[np.abs(Y) > 1e-6*np.abs(Y).max()].max()
                out[nm] = dict(chi=chi, Y=Y, res=res, w0=w0, l1=l1, l1Y=l1Y, nu=S.nu)
                pr(f"    {nm}: residual {res:.1e}; int chi Nhat = {w0:+.10f} (target {Phi(np.array([0.0]))[0]:+.1f}); "
                   f"||chi Nhat||_1 = {l1:.3f}, int|Y| = {l1Y:.3e}; |Y| > 1e-6 max up to nu = {supp:.1f}; "
                   f"W(t) direct at t = 0,.5,1,2,3: " + " ".join(f"{v:+.8f}" for v in Wd))
            res_by_grid[(dnu, N)] = (S, out)
        # refinement
        (S1, o1), (S2, o2) = res_by_grid[(0.02, 2048)], res_by_grid[(0.02, 3072)]
        for nm in TGT:
            d = np.max(np.abs(o1[nm]['chi'] - o2[nm]['chi'][:2048]))/np.max(np.abs(o2[nm]['chi']))
            pr(f"  (R) {nm}: nu_max 41 -> 61: max|d chi|/max|chi| = {d:.1e}")
        S0, o0 = res_by_grid[(0.04, 1024)]
        for nm in TGT:
            ci = np.interp(S0.nu, S1.nu, o1[nm]['chi'])
            d = np.max(np.abs(o0[nm]['chi'] - ci))/np.max(np.abs(ci))
            pr(f"  (R) {nm}: dnu 0.04 -> 0.02: max|d chi|/max|chi| = {d:.1e}")
        # profile
        S, out = res_by_grid[(0.02, 2048)]
        pr("  profile (dnu = 0.02): nu, Nhat(nu)/|Nhat(0)|, chi_value*Nhat, chi_deriv*Nhat")
        N0 = abs(S.Nhat[0])
        for v in (0.01, 0.5, 1.0, 1.0685, 1.5, 2, 3, 4.546, 5, 6, 8.169, 10, 13.094, 15, 20, 25, 30):
            k = min(S.N - 1, int(v/S.dnu))
            pr(f"     {S.nu[k]:7.3f}  {S.Nhat[k]/N0:+.4e}  {out['value']['chi'][k]*S.Nhat[k]:+.4e}  {out['deriv']['chi'][k]*S.Nhat[k]:+.4e}")
        # (V2) independent image check
        pr("  (V2) direct image (mpmath complex-degree Legendre, nu <= 30, stride 2) vs target:")
        for beta in (0.0, 0.3, 1.0, 2.0, 4.0, 6.0):
            for nm in TGT:
                tgt = float(mp.legenp(-0.5, 0, mp.cosh(beta), type=3).real) if nm == 'value' else -Dbeta(beta)
                if beta == 0.0 and nm == 'deriv':
                    tgt = 0.0
                img = image_direct(S, out[nm]['chi'], beta, numax=30, stride=2)
                sc = abs(tgt) if tgt != 0 else 1.0
                pr(f"     beta = {beta:4.1f} {nm:5s}: image {img.real:+.8e} (Im {img.imag:+.1e})  target {tgt:+.8e}  rel.diff {abs(img - tgt)/sc:.1e}")
        # null vectors
        U, s, Vh = S.svd()
        nn = int(np.sum(s < 1e-9*s[0]))
        for kk in range(nn):
            v = Vh[-1 - kk]
            a = np.abs(v)
            pr(f"  null vector {kk + 1}: sv/s_max {s[-1 - kk]/s[0]:.1e}; peak nu = {S.nu[np.argmax(a)]:.2f}; |Y_h| > 1e-3 max up to nu = "
               f"{S.nu[a > 1e-3*a.max()].max():.1f}")
        pickle.dump(dict(nu=S.nu, chi_value=out['value']['chi'], chi_deriv=out['deriv']['chi'], Y_value=out['value']['Y'],
                         Y_deriv=out['deriv']['Y'], Nhat=S.Nhat, n=S.n, sv=s, dnu=S.dnu),
                    open(os.path.join(CACHE, f'stage4_deshift_{key}.pkl'), 'wb'))
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'stage4_deshift.txt'), 'w').write(buf.getvalue())
