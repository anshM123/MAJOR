"""Result Q, Stage 4 (RESULT_Q_STAGE4.md §3): the k = 5 forward piece.  Stage 3 called raw C^_5(0) "fine (J0 < 3/2)" -- true for
CONVERGENCE, but its leading-Regge heavy image is not a pure moment: the k = 5 kernel is -16 eta^2 r^2 x (k = 3 kernel) (Lemma B), the
k = 3 image is a pure P_l (D, Stage 12 §1.4), so the raw k = 5 image is sum_{j = 0, +-2} h_j P_{l+j}; at nu = 0 it grows like
beta P_{3/2}(cosh beta) ~ eta^{3/2} log eta, with the helicity sign sigma_h of B5.  Stage 3/12's heavy models used the pure moment
sigma_h u^3 P0 instead, i.e. implicitly a DE-SHIFTED C^_5.
(a) raw image, scalar proxy h_j = a^(2)_j(l) E_b(nu - i j) (the proxy is exact for the scalar rank-3 case, stage4_hj (C2)); b = 8, 10
    bracket the unknown radial weight of the j' = 1 kernel; the spinning (j' = 1, o_H q) coefficients are NOT computed [I].
(b) the leading-Regge model of stage12_forward_regge with the raw C^_5 image in place of P0: where does the Stage-0/3 functional go negative?
(c) the 3-shift de-shift of C^_5 (proxy coefficients): W[chi](t) = C[chi h0] + cosh 2t C[chi A2] - sinh 2t S[chi B2] = 1.
Output: cache/stage4_c5.txt.  ~ minutes.
"""
import os, sys, io, json, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
from scipy.special import loggamma
import stage12_forward_regge as FR

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def a2(l):
    """eta^2 P_l = sum_j a_j P_{l+j} (Legendre recurrence twice)."""
    c = {0: 1}
    for _ in range(2):
        n = {}
        for j, v in c.items():
            m = l + j
            n[j + 1] = n.get(j + 1, 0) + v*(m + 1)/(2*m + 1)
            n[j - 1] = n.get(j - 1, 0) + v*m/(2*m + 1)
        c = n
    return c


def hj_proxy(nu, b):
    """h_j(nu), j = -2, 0, 2 (scalar proxy; common nu-only factor removed by E_b(nu))."""
    l = -0.5 + 1j*np.asarray(nu, complex)
    A = a2(l)
    lE = lambda mu: 2*(loggamma((b - 0.5 - 1j*mu)/2) + loggamma((b - 0.5 + 1j*mu)/2))
    return {j: A[j]*np.exp(lE(nu - 1j*j) - lE(nu + 0j)) for j in (-2, 0, 2)}


def raw_image(beta, b, nu=1e-7):
    h = hj_proxy(np.array([nu]), b)
    x = mp.cosh(beta)
    num = sum(complex(h[j][0])*complex(mp.legenp(-0.5 + 1j*(nu - 1j*j), 0, x, type=3)) for j in h)
    den = sum(complex(h[j][0]) for j in h)
    return (num/den).real


def load_func(name):
    for fn in ('functionals_q.json', 'functionals_q3.json'):
        D = json.load(open(os.path.join(HERE, fn)))
        if name in D:
            F = D[name]
            return F, np.array(F['psi2'], float), np.array(F['psi3'], float) if len(F['psi3']) else np.zeros(1)
    raise KeyError(name)


def lr_scan(Lam, name, H5fun, betas, us):
    """stage12_forward_regge part (b) model; B5 image sigma u^3 x H5(beta) (H5 = P0 reproduces stage 12)."""
    F, c2, c3 = load_func(name)
    a, c = F['a'], F['c']
    C4 = {+1: 1.0, -1: 41.0}
    omax = Lam*betas.max() + 10
    om = np.linspace(0, omax, int(omax*8) + 2)
    T = FR.cos_transforms(c2, c3, om)
    worst = (np.inf, None); first_neg = None
    for beta in betas:
        t, W = FR.mehler_nodes(beta, 1600 if Lam*beta > 2000 else 800)
        tr = [np.interp(Lam*t, om, T[i]) for i in range(3)]
        I20 = (W*tr[0]).sum(); I22 = (W*tr[1]).sum(); I3 = (W*tr[2]).sum()
        P0 = W.sum(); D = 0.5*(W*t**2).sum()
        H5 = H5fun(beta, P0)
        for h in (+1, -1):
            jmin = 0 if h > 0 else 4
            phys = Lam*beta >= (2*jmin + 1)*np.sqrt(us)
            B2 = 2*I20 - us*I22
            B3 = h*us*I3
            FW = a*(2*us**2*Lam**2*D - C4[h]*us**3*P0 + h*us**3*H5) - 2*c*us**2*P0
            FF = B2 + B3 + FW
            Nrm = np.abs(B2) + np.abs(B3) + a*(2*us**2*Lam**2*D + 42*us**3*P0 + us**3*abs(H5)) + 2*c*us**2*P0
            rat = np.where(phys, FF/Nrm, np.inf)
            i = int(np.argmin(rat))
            if Lam*beta < 20:                       # b < 20: outside the leading-Regge model's validity (Stage 12 §1.5)
                continue
            if rat[i] < worst[0]:
                worst = (rat[i], (beta, us[i], h))
            if first_neg is None and rat[i] < 0:
                first_neg = (beta, us[i], h)
    return worst, first_neg


def deshift3(b, dnu=0.02, N=2048):
    nu = (np.arange(N) + 0.5)*dnu
    t = np.arange(N)*np.pi/(N*dnu)
    h = hj_proxy(nu, b)
    R0 = h[0]; A2 = h[2] + h[-2]; B2 = -1j*(h[2] - h[-2])
    imag = max(np.max(np.abs(x.imag)/np.abs(x)) for x in (R0, A2, B2))
    R0, A2, B2 = R0.real, A2.real, B2.real
    n = np.abs(h[0]) + np.abs(h[2]) + np.abs(h[-2])
    Cm = dnu*np.cos(np.outer(t, nu)); Sm = dnu*np.sin(np.outer(t, nu))
    r0 = 1/np.cosh(np.minimum(2*t, 700))
    M = Cm*(A2/n)[None, :] - np.tanh(2*t)[:, None]*Sm*(B2/n)[None, :] + r0[:, None]*Cm*(R0/n)[None, :]
    U, s, Vh = np.linalg.svd(M)
    bvec = r0.copy()
    keep = s > 1e-11*s[0]
    Y = Vh[keep].T @ ((U[:, keep].T @ bvec)/s[keep])
    res = np.linalg.norm(M @ Y - bvec)/np.linalg.norm(bvec)
    chi = Y/n
    Nhat = R0 + A2
    Wt = []
    for tt in (0, 1, 2, 4, 6):
        cc, ss = np.cos(nu*tt), np.sin(nu*tt)
        Wt.append(dnu*np.sum(chi*(R0*cc + np.cosh(2*tt)*A2*cc - np.sinh(2*tt)*B2*ss)))
    zeros = nu[1:][np.sign(Nhat[1:]) != np.sign(Nhat[:-1])]
    return dict(imag=imag, sv=s[::-1][:5]/s[0], res=res, W=Wt, l1=dnu*np.sum(np.abs(chi*Nhat)), zeros=zeros, nnull=int(np.sum(~keep)))


if __name__ == '__main__':
    t0 = time.time()
    pr("# stage4_c5: raw vs de-shifted k = 5 forward piece (scalar proxy for the eta^2 mixing)")
    pr("(a) raw C^_5(0) leading-Regge image H5(beta) (normalised H5(0) = 1), proxy radial weight b:")
    pr("    beta    P0=P_{-1/2}     H5 (b=8)        H5 (b=10)       H5/(beta P_{3/2}) b=8, b=10")
    for beta in (0.01, 0.5, 1, 2, 4, 6, 8, 12):
        P0 = float(mp.legenp(-0.5, 0, mp.cosh(beta), type=3).real); P32 = float(mp.legenp(1.5, 0, mp.cosh(beta), type=3).real)
        h8 = raw_image(beta, 8); h10 = raw_image(beta, 10)
        pr(f"   {beta:5.2f}  {P0:.6e}  {h8:+.6e}  {h10:+.6e}   {h8/(beta*P32):+.4f} {h10/(beta*P32):+.4f}")
    # (b)
    pr("(b) leading-Regge model (stage12_forward_regge (b), c = 1/2): min of FF/Nrm over physical states with bhat = Lam beta >= 20, beta <= 20, u in [1e-4, 1];")
    pr("    'pure' = de-shifted C^_5 (B5 -> sigma u^3 P0, = Stage 12/3 model); 'raw' = raw C^_5 (B5 -> sigma u^3 H5, proxy b = 8)")
    betas = np.concatenate([np.geomspace(2e-3, 0.05, 15), np.linspace(0.06, 20, 160)])
    us = np.concatenate([np.geomspace(1e-4, 0.05, 15), np.linspace(0.06, 1.0, 40)])
    Hcache = {}
    def H5raw(beta, P0):
        k = round(float(beta), 10)
        if k not in Hcache:
            Hcache[k] = raw_image(beta, 8)
        return Hcache[k]
    for Lam in (100., 300., 1000.):
        for name in ('Q0-full-robust', 'Q3-full-delta0.2'):
            wp, fp = lr_scan(Lam, name, lambda beta, P0: P0, betas, us)
            wr, fr = lr_scan(Lam, name, H5raw, betas, us)
            def J_of(beta, u):
                m = Lam/np.sqrt(u)
                return 0.5*(m*np.sqrt(2*(np.cosh(beta) - 1)) - 1)
            fneg = (f"first negative at beta = {fr[0]:.2f}, u = {fr[1]:.3f} ({'++' if fr[2] > 0 else '+-'}), m = {Lam/np.sqrt(fr[1]):.0f}, J ~ {J_of(fr[0], fr[1]):.3g}"
                    if fr else "none")
            pr(f"   Lam = {Lam:6.0f} {name:18s}: pure min {wp[0]:+.4f} (beta {wp[1][0]:.2f});  raw min {wr[0]:+.4f} (beta {wr[1][0]:.2f}, u {wr[1][1]:.3f}, "
               f"{'++' if wr[1][2] > 0 else '+-'}); {fneg}")
    # (c)
    pr("(c) 3-shift de-shift of C^_5 (proxy coefficients), target Phi = 1 (pure P0 image):")
    for b in (8, 10):
        R = deshift3(b)
        pr(f"   b = {b}: Im/|.| {R['imag']:.1e}; smallest sv/s_max " + " ".join(f"{v:.1e}" for v in R['sv']) +
           f"; null dim {R['nnull']}; residual {R['res']:.1e}; W(t = 0,1,2,4,6) = " + " ".join(f"{v:.8f}" for v in R['W']) +
           f"; ||chi Nhat||_1 = {R['l1']:.3f}; real zeros of the forward normalisation (nu < 41): {np.round(R['zeros'], 3).tolist()}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage4_c5.txt'), 'w').write(buf.getvalue())
