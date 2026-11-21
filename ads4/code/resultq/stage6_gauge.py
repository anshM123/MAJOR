"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §2): gauge robustness of the two de-shifted forward pieces.

Gauge families (every member: same Abel equation W = Phi, hence the same leading-Regge image and the same light-side moments):
  C5 (3-shift, complex c, stage5_deshift.Deshift3S general, Y = (Y_R, Y_I), c = (Y_R + i Y_I)/n):
     G0 min-norm (Stage 5 reference), G1 (1+(nu/5)^2), G2 Im x10, G3 Re x10   (Stage 5's four)
     G4 smooth: min ||Y'||^2 + 1e-2 ||Y||^2 (first differences in nu, per block)
     G5 amp-L2: min int (|c| w)^2,  w(nu) = sum_j |hh_j(nu)| (the amplification density of a per-harmonic error)
     G6 amp-L1: min int |c| w  (LP over the null space, |c| by a 16-gon; minimal worst-case amplification)
  C4 (4-shift, real chi, stage4_deshift.Deshift4, Y = chi n, n = sum_j |R_j|; value Phi = 1 and derivative Phi = -t^2/2):
     H0 min-norm tol 1e-11 (Stage 4 reference; = amp-L2 since the amplification density of chi is n|chi| = |Y|), H1 tol 1e-6 (drops the
     ill-conditioned direction), H2 (1+(nu/5)^2), H3 smooth, H4 amp-L1 (LP).
Per gauge: residual, W(t) check, shifted moments M_n = Re sum_j int c hh_j (nu - ij)^{2n} (n <= 4; target delta_{n0} resp. delta_{n1}),
amplifications A0 = int |c| w, A2 = int |c| sum_j |hh_j| |nu - ij|^2, distance to the reference, and the NON-COVARIANT probes
   mu_n  = Re int c (sum_j hh_j) nu^{2n}            ("no-shift" model NS1: the per-harmonic argument nu - ij replaced by nu)
   kap_n = Re int c sum_j (j/2) hh_j nu^{2n}         ("j-odd" model NS2)
which are the only gauge-dependent numbers entering the stress models of stage6_minima.py.  Covariant test: for a flat scalar-skeleton
g(mu^2) = (2-y)(1-y)^-3 P_J(1-2y), y = mu^2/m^2 (entire up to |mu| = m) at the binding states, Re sum_j int c hh_j g((nu-ij)^2) must equal
g(0) (value; C5) resp. g'(0) (derivative) in EVERY gauge.
Output: cache/stage6_gauge.txt, cache/stage6_gauge.pkl
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import linprog
from scipy.special import eval_legendre
import stage5_deshift as DS3
import stage4_deshift as DS4

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def constrained_min(M, b, L, tol=1e-11):
    """min ||L Y||_2 s.t. M Y = b (least squares, s.v. < tol s_max dropped); L square invertible."""
    Li = np.linalg.inv(L)
    U, s, Vh = np.linalg.svd(M @ Li, full_matrices=False)
    keep = s > tol*s[0]
    Z = Vh[keep].T @ ((U[:, keep].T @ b)/s[keep])
    return Li @ Z


def diff_op(N, blocks, eps=1e-2, dnu=1.0):
    D = np.zeros((N*blocks, N*blocks))
    for bl in range(blocks):
        o = bl*N
        for k in range(N):                      # L Y = sqrt(eps) Y_k + (Y_k - Y_{k+1}): inverse decays like (1+sqrt eps)^-k
            D[o + k, o + k] = 1.0 + np.sqrt(eps)
            if k + 1 < N:
                D[o + k, o + k + 1] = -1.0
    return D


def l1_null(M, b, wvec, blocks, N, ngon=16, tol=1e-11):
    """min sum_k w_k |Y_k| (complex if blocks = 2: Y_k = Y_R,k + i Y_I,k, |.| by an outer ngon) s.t. M Y = b, via Y = Y0 + V0 z."""
    U, s, Vh = np.linalg.svd(M, full_matrices=True)
    r = int(np.sum(s > tol*s[0]))
    Y0 = Vh[:r].T @ ((U[:, :r].T @ b)/s[:r])
    V0 = Vh[r:].T
    nz = V0.shape[1]
    th = np.arange(ngon)*2*np.pi/ngon
    rows, rhs = [], []
    if blocks == 2:
        for t in th:
            # w_k (cos t (Y0R + V0R z) + sin t (Y0I + V0I z)) <= T_k
            A = wvec[:, None]*(np.cos(t)*V0[:N] + np.sin(t)*V0[N:])
            rows.append(np.hstack([A, -np.eye(N)])); rhs.append(-wvec*(np.cos(t)*Y0[:N] + np.sin(t)*Y0[N:]))
    else:
        for sg in (1, -1):
            A = sg*wvec[:, None]*V0
            rows.append(np.hstack([A, -np.eye(N)])); rhs.append(-sg*wvec*Y0)
    Aub = np.vstack(rows); bub = np.concatenate(rhs)
    obj = np.concatenate([np.zeros(nz), np.ones(N)/np.cos(np.pi/ngon if blocks == 2 else 0)])
    res = linprog(obj, A_ub=Aub, b_ub=bub, bounds=[(None, None)]*nz + [(0, None)]*N, method='highs')
    if res.status != 0:
        return None, res.message, nz
    return Y0 + V0 @ res.x[:nz], 'ok', nz


def gflat(J, m, nmax=4):
    """Taylor coefficients g_n (in mu^2) of (2-y)(1-y)^-3 P_J(1-2y), y = mu^2/m^2, n <= nmax."""
    PJ = np.array([(-1)**k*np.exp(sum(np.log((J - i)*(J + 1 + i)/((i + 1)**2)) for i in range(k))) if k <= J else 0. for k in range(nmax + 1)])
    inv = np.array([(n + 1)*(n + 2)/2 for n in range(nmax + 1)])
    a = np.convolve(np.convolve([2., -1.], inv)[:nmax + 1], PJ)[:nmax + 1]
    return a/m**(2*np.arange(nmax + 1))


def gfun(J, m, mu2):
    y = mu2/m**2
    return (2 - y)*(1 - y)**-3*eval_legendre(J, 1 - 2*y)


def c5_family(dnu, N):
    S = DS3.Deshift3S('C5_pp', dnu, N, general=True)
    b = np.ones_like(S.t)/S.ch2
    w = sum(np.abs(S.hh[j]) for j in (0, 2, -2))
    one = np.ones(S.N)
    G = {}
    for k, d in (('G0 min-norm', np.ones(2*N)), ('G1 (1+(nu/5)^2)', np.concatenate([1 + (S.nu/5)**2]*2)),
                 ('G2 Im x10', np.concatenate([one, 10*one])), ('G3 Re x10', np.concatenate([10*one, one])),
                 ('G5 amp-L2', np.concatenate([w/S.n]*2))):
        G[k] = constrained_min(S.M, b, np.diag(d))
    G['G4 smooth'] = constrained_min(S.M, b, diff_op(N, 2))
    t0 = time.time()
    Y, msg, nz = l1_null(S.M, b, w/S.n, 2, N, ngon=16 if N <= 512 else 8)      # 8-gon at N = 1024 (memory <= 3 GB)
    if Y is not None:
        G['G6 amp-L1'] = Y
    pr(f"    [C5 dnu={dnu} N={N}: amp-L1 LP over {nz} null directions: {msg} ({time.time() - t0:.0f}s)]")
    out = {}
    for k, Y in G.items():
        c = (Y[:N] + 1j*Y[N:])/S.n
        out[k] = dict(c=c, res=np.linalg.norm(S.M @ Y - b)/np.linalg.norm(b))
    return S, w, out


def c5_report(S, w, out, ref='G0 min-norm'):
    c0 = out[ref]['c']; d = S.dnu
    for k, o in out.items():
        c = o['c']
        Wt = S.W(c, np.array([0., 1., 3., 6.]))
        M = [d*np.sum(sum(c*S.hh[j]*(S.nu - 1j*j)**(2*n) for j in (0, 2, -2))).real for n in range(5)]
        A0 = d*np.sum(np.abs(c)*w); A2 = d*np.sum(np.abs(c)*sum(np.abs(S.hh[j])*np.abs(S.nu - 1j*j)**2 for j in (0, 2, -2)))
        mu = [d*np.sum((c*sum(S.hh[j] for j in (0, 2, -2))).real*S.nu**(2*n)) for n in range(4)]
        kap = [d*np.sum((c*sum(j/2*S.hh[j] for j in (0, 2, -2))).real*S.nu**(2*n)) for n in range(4)]
        dist = np.sum(np.abs(c - c0))/np.sum(np.abs(c0)); amp_d = d*np.sum(np.abs(c - c0)*w)
        o.update(M=M, A0=A0, A2=A2, mu=mu, kap=kap, dist=dist, ampd=amp_d, intc=d*np.sum(np.abs(c)))
        pr(f"    {k:16s} res {o['res']:.0e} W(0,1,3,6)-1 " + " ".join(f"{v - 1:+.0e}" for v in Wt) +
           f" | M1..4 " + " ".join(f"{v:+.0e}" for v in M[1:]) + f" | int|c| {o['intc']:.3f} A0 {A0:.3f} A2 {A2:.1f}"
           f" | dist {dist:.2f} int|dc|w {amp_d:.3f} | mu1,2 {mu[1]:+.3f} {mu[2]:+.2f} | kap0,1,2 {kap[0]:+.3f} {kap[1]:+.3f} {kap[2]:+.2f}")


def c4_family(dnu, N, Phi, deriv):
    S = DS4.Deshift4('EE3', dnu, N)
    b = Phi(S.t)/S.ch3
    G = {}
    G['H0 min-norm'] = constrained_min(S.M, b, np.eye(N))
    G['H1 tol 1e-6'] = constrained_min(S.M, b, np.eye(N), tol=1e-6)
    G['H2 (1+(nu/5)^2)'] = constrained_min(S.M, b, np.diag(1 + (S.nu/5)**2))
    G['H3 smooth'] = constrained_min(S.M, b, diff_op(N, 1))
    t0 = time.time()
    Y, msg, nz = l1_null(S.M, b, np.ones(N), 1, N)
    if Y is not None:
        G['H4 amp-L1'] = Y
    pr(f"    [C4 {'deriv' if deriv else 'value'} dnu={dnu} N={N}: amp-L1 LP over {nz} null directions: {msg} ({time.time() - t0:.0f}s)]")
    out = {}
    for k, Y in G.items():
        out[k] = dict(chi=Y/S.n, Y=Y, res=np.linalg.norm(S.M @ Y - b)/np.linalg.norm(b))
    return S, out


def c4_report(S, out, Phi, ref='H0 min-norm'):
    d = S.dnu; Y0 = out[ref]['Y']
    for k, o in out.items():
        chi = o['chi']
        tt = np.array([0., 1., 3., 6.])
        Wt = S.W(chi, tt)
        M = [d*np.sum(chi*sum(S.Rj[j]*(S.nu - 1j*j)**(2*n) for j in (1, -1, 3, -3))).real for n in range(5)]
        A0 = d*np.sum(np.abs(o['Y'])); A2 = d*np.sum(np.abs(chi)*sum(np.abs(S.Rj[j])*np.abs(S.nu - 1j*j)**2 for j in S.Rj))
        mu = [d*np.sum(chi*S.Nhat*S.nu**(2*n)) for n in range(4)]
        kap = [d*np.sum((chi*sum(j/2*S.Rj[j] for j in S.Rj)).real*S.nu**(2*n)) for n in range(4)]
        dist = np.sum(np.abs(o['Y'] - Y0))/np.sum(np.abs(Y0))
        o.update(M=M, A0=A0, A2=A2, mu=mu, kap=kap, dist=dist, ampd=d*np.sum(np.abs(o['Y'] - Y0)))
        pr(f"    {k:16s} res {o['res']:.0e} W(0,1,3,6)-Phi " + " ".join(f"{v - Phi(np.array([x]))[0]:+.0e}" for v, x in zip(Wt, tt)) +
           f" | M0..4 " + " ".join(f"{v:+.1e}" for v in M) + f" | A0 {A0:.3f} A2 {A2:.1f} | dist {dist:.2f} int|dY| {o['ampd']:.3f}"
           f" | mu1,2 {mu[1]:+.3f} {mu[2]:+.2f} | kap0,1,2 {kap[0]:+.3f} {kap[1]:+.3f} {kap[2]:+.2f}")


def covariant_test(S, fam, kind, states):
    """Re sum_j int c hh_j g((nu - ij)^2): must be g(0) (value/C5) or g'(0) (deriv) in every gauge."""
    for (lab, J, m) in states:
        gn = gflat(J, m, 2)
        vals = []
        for k, o in fam.items():
            if o['res'] > 1e-8:
                continue
            if kind == 'C5':
                tot = sum(o['c']*S.hh[j]*gfun(J, m, (S.nu - 1j*j)**2) for j in (0, 2, -2))
            else:
                tot = o['chi']*sum(S.Rj[j]*gfun(J, m, (S.nu - 1j*j)**2) for j in S.Rj)
            vals.append(S.dnu*np.sum(tot).real)
        tgt = gn[1] if kind == 'C4d' else gn[0]
        vals = np.array(vals)
        pr(f"      covariant test {kind} {lab} (J={J}, m={m:.1f}): target {tgt:+.10e}; max rel dev over gauges {np.max(np.abs(vals/tgt - 1)):.1e};"
           f" spread {np.ptp(vals)/abs(tgt):.1e}")


if __name__ == '__main__':
    t00 = time.time()
    pr("# stage6_gauge: gauge families of the de-shifted forward pieces (C5 3-shift complex; C4 4-shift value/derivative)")
    L100 = np.sqrt(98*99.)
    states = [('Dg100 +- J5 1.2Lam', 5, 1.2*L100), ('Dg100 +- J4 gap', 4, np.sqrt(98*(99 + 8.))), ('Dg100 +- J8 1.1Lam', 8, 1.1*L100),
              ('Dg100 J24 1.3Lam', 24, 1.3*L100), ('Dg100 J32 3Lam', 32, 3*L100)]
    res = {}
    for (dnu, N) in ((0.08, 512), (0.04, 1024)):
        pr(f"\n## C5 (spinning 3-shift), dnu = {dnu}, N = {N}")
        S, w, fam = c5_family(dnu, N)
        c5_report(S, w, fam)
        covariant_test(S, fam, 'C5', states)
        res[('C5', dnu)] = {k: {kk: vv for kk, vv in o.items()} for k, o in fam.items()}
        res[('C5', dnu, 'nu')] = S.nu; res[('C5', dnu, 'w')] = w
        del S
    for (dnu, N) in ((0.04, 1024), (0.02, 2048)):
        for deriv in (False, True):
            Phi = (lambda t: -t**2/2) if deriv else (lambda t: np.ones_like(t))
            pr(f"\n## C4 (EE 4-shift) {'derivative' if deriv else 'value'}, dnu = {dnu}, N = {N}")
            S, fam = c4_family(dnu, N, Phi, deriv)
            c4_report(S, fam, Phi)
            covariant_test(S, fam, 'C4d' if deriv else 'C4v', states)
            res[('C4d' if deriv else 'C4v', dnu)] = fam; res[('C4', dnu, 'nu')] = S.nu
            del S
    pickle.dump(res, open(os.path.join(CACHE, 'stage6_gauge.pkl'), 'wb'))
    pr(f"[{time.time() - t00:.0f}s]")
    open(os.path.join(CACHE, 'stage6_gauge.txt'), 'w').write(buf.getvalue())
