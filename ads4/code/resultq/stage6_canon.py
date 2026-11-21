"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §2.2): a canonical ('shift-blind', SB) gauge for both de-shifted forward pieces.

Extra condition.  The leading-Regge image is exactly covariant (sum_j R_j(nu) P(nu - ij), stages 4/5), and every gauge maps any covariant image to
the same pure moment.  The first structure that can distinguish gauges is a beyond-LR term that is NOT a function of (nu - ij): at O(1/m^2) the
simplest one is the shift-blind quadratic  dF_j(nu) ∝ nu^2 - (nu - ij)^2  (the finite-J constants read with the unshifted argument).  Its image is
proportional to  s := mu_1 - M_1,  mu_1 = Re sum_j int c h_j nu^2 (unshifted), M_1 = Re sum_j int c h_j (nu - ij)^2 (shifted; = 0 value/C5, = 1
derivative).  SB gauge := min ||Y||_2 subject to the Abel equation AND s = 0 (one extra linear condition).  For C4 the condition is met through
the near-null direction (s.v. ~2e-11; Stage 4's 'ill-conditioned' one), hence tol 1e-9 for the C4 null space; for C5 through its large null space.
Checks: residual, W(t), shifted moments, covariant test (stage6_gauge.covariant_test), amplification A0, refinement (C4 dnu .04 -> .02,
C5 .08 -> .04), and the second shift-blind probe s_b (j^2-weighted, NOT imposed) for comparison with the other gauges.
Writes functionals_q6_canonical.json (SB weights, C4 value/derivative on dnu = 0.02 N = 2048, C5 on dnu = 0.04 N = 1024) and
cache/stage6_canon.{txt,pkl}.
"""
import os, sys, io, json, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import stage5_deshift as DS3
import stage4_deshift as DS4
import stage6_gauge as G6

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def sb_solve(M, b, r, t, tol):
    """min ||Y|| s.t. M Y = b (s.v. < tol s_max treated as null) and r.Y = t."""
    U, s, Vh = np.linalg.svd(M, full_matrices=True)
    k = int(np.sum(s > tol*s[0]))
    Y0 = Vh[:k].T @ ((U[:, :k].T @ b)/s[:k])
    V0 = Vh[k:].T
    g = V0.T @ r
    z = g*(t - r @ Y0)/(g @ g)
    return Y0 + V0 @ z, Y0, V0.shape[1]


def c4_sb(dnu, N, deriv):
    S = DS4.Deshift4('EE3', dnu, N)
    Phi = (lambda t: -t**2/2) if deriv else (lambda t: np.ones_like(t))
    b = Phi(S.t)/S.ch3
    r = dnu*S.Nhat*S.nu**2/S.n                      # mu_1 as a function of Y (chi = Y/n)
    Y, Y0, nz = sb_solve(S.M, b, r, 1.0 if deriv else 0.0, 1e-9)
    return S, Y, Y0, nz, Phi, b


def c5_sb(dnu, N):
    S = DS3.Deshift3S('C5_pp', dnu, N, general=True)
    b = np.ones_like(S.t)/S.ch2
    r = np.concatenate([dnu*S.nu**2/S.n, np.zeros(N)])     # mu_1 = int Re(c) nu^2 (sum_j hh_j = 1); c = (Y_R + i Y_I)/n
    Y, Y0, nz = sb_solve(S.M, b, r, 0.0, 1e-11)
    return S, Y, Y0, nz, b


def probes4(S, chi):
    d = S.dnu
    M = [d*np.sum(chi*sum(S.Rj[j]*(S.nu - 1j*j)**(2*n) for j in S.Rj)).real for n in range(4)]
    mu1 = d*np.sum(chi*S.Nhat*S.nu**2)
    sb = d*np.sum(chi*sum((j*j/9.)*S.Rj[j] for j in S.Rj)).real*0 + d*np.sum((chi*sum((j*j/9.)*S.Rj[j] for j in S.Rj)).real*S.nu**2)
    return M, mu1, sb, d*np.sum(np.abs(chi*S.n)), d*np.sum(np.abs(chi)*sum(np.abs(S.Rj[j])*np.abs(S.nu - 1j*j)**2 for j in S.Rj))


def probes5(S, c):
    d = S.dnu
    M = [d*np.sum(sum(c*S.hh[j]*(S.nu - 1j*j)**(2*n) for j in (0, 2, -2))).real for n in range(4)]
    mu1 = d*np.sum(c.real*S.nu**2)
    sb = d*np.sum((c*sum((j*j/4.)*S.hh[j] for j in (0, 2, -2))).real*S.nu**2)
    w = sum(np.abs(S.hh[j]) for j in (0, 2, -2))
    return M, mu1, sb, d*np.sum(np.abs(c)*w), d*np.sum(np.abs(c)*sum(np.abs(S.hh[j])*np.abs(S.nu - 1j*j)**2 for j in (0, 2, -2)))


if __name__ == '__main__':
    t0 = time.time()
    L100 = np.sqrt(98*99.)
    states = [('Dg100 +- J5 1.2Lam', 5, 1.2*L100), ('Dg100 +- J4 gap', 4, np.sqrt(98*(99 + 8.))), ('Dg100 J24 1.3Lam', 24, 1.3*L100)]
    out = {}
    pr("# stage6_canon: shift-blind (SB) gauge = min-norm + [mu_1 = M_1]")
    for deriv in (False, True):
        lab = 'C4d' if deriv else 'C4v'
        sols = {}
        for (dnu, N) in ((0.04, 1024), (0.02, 2048)):
            S, Y, Y0, nz, Phi, b = c4_sb(dnu, N, deriv)
            chi = Y/S.n; chi0 = Y0/S.n
            res = np.linalg.norm(S.M @ Y - b)/np.linalg.norm(b)
            Wt = S.W(chi, np.array([0., 1., 3., 6., 8.]))
            for tag, cc in (('min-norm(tol 1e-9)', chi0), ('SB', chi)):
                M, mu1, sb, A0, A2 = probes4(S, cc)
                pr(f"  {lab} dnu={dnu} N={N} {tag:18s}: null dim {nz}; M0..3 " + " ".join(f"{v:+.2e}" for v in M) +
                   f"; mu1 {mu1:+.4f}; s_a = mu1 - M1 {mu1 - M[1]:+.2e}; s_b {sb:+.4f}; A0 {A0:.3f}; A2 {A2:.1f}")
            pr(f"      SB residual {res:.1e}; W(0,1,3,6,8) - Phi = " + " ".join(f"{v - Phi(np.array([x]))[0]:+.1e}" for v, x in zip(Wt, (0, 1, 3, 6, 8))))
            fam = {'SB': dict(chi=chi, res=res)}
            G6.covariant_test(S, fam, lab, states)
            sols[dnu] = (S.nu, chi, S)
        nu1, c1, S1 = sols[0.04]; nu2, c2, S2 = sols[0.02]
        ci = np.interp(nu1, nu2, c2*S2.n)
        pr(f"    refinement dnu .04 -> .02: int|dY|/int|Y| = {np.sum(np.abs(c1*S1.n - ci))/np.sum(np.abs(c1*S1.n)):.2e}")
        out[lab] = dict(nu=nu2, chi=c2, nu_04=nu1, chi_04=c1)
    sols = {}
    for (dnu, N) in ((0.08, 512), (0.04, 1024)):
        S, Y, Y0, nz, b = c5_sb(dnu, N)
        c = (Y[:N] + 1j*Y[N:])/S.n; c0 = (Y0[:N] + 1j*Y0[N:])/S.n
        res = np.linalg.norm(S.M @ Y - b)/np.linalg.norm(b)
        Wt = S.W(c, np.array([0., 1., 3., 6., 8.]))
        for tag, cc in (('min-norm', c0), ('SB', c)):
            M, mu1, sb, A0, A2 = probes5(S, cc)
            pr(f"  C5 dnu={dnu} N={N} {tag:18s}: null dim {nz}; M0..3 " + " ".join(f"{v:+.2e}" for v in M) +
               f"; mu1 {mu1:+.4f}; s_a {mu1 - M[1]:+.2e}; s_b {sb:+.4f}; A0 {A0:.3f}; A2 {A2:.1f}; int|c| {dnu*np.sum(np.abs(cc)):.3f}")
        pr(f"      SB residual {res:.1e}; W(0,1,3,6,8) - 1 = " + " ".join(f"{v - 1:+.1e}" for v in Wt))
        G6.covariant_test(S, {'SB': dict(c=c, res=res)}, 'C5', states)
        sols[dnu] = (S.nu, c)
    nu1, c1 = sols[0.08]; nu2, c2 = sols[0.04]
    ci = np.interp(nu1, nu2, c2.real) + 1j*np.interp(nu1, nu2, c2.imag)
    pr(f"    refinement dnu .08 -> .04: int|dc|/int|c| = {np.sum(np.abs(c1 - ci))/np.sum(np.abs(c1)):.2e}")
    out['C5'] = dict(nu=nu2, c=c2)
    pickle.dump(out, open(os.path.join(CACHE, 'stage6_canon.pkl'), 'wb'))
    J = {'description': "Result Q Stage 6 shift-blind (SB) canonical gauge for the de-shifted forward pieces: min-norm subject to the Abel "
                        "equation W = Phi and mu_1 = M_1 (unshifted nu^2-moment = shifted one). C4: int chi(nu) Psi_{4,nu}[EE] dnu "
                        "(stage4 conventions, R_j from functionals_q4_deshift.json; value Phi = 1, derivative Phi = -t^2/2). C5: int Re[c(nu) "
                        "Psi_{5,nu}/Nhat5(nu)] dnu (stage5 conventions, functionals_q5_deshift.json).",
         'C4_value': {'dnu': 0.02, 'nu': out['C4v']['nu'].tolist(), 'chi': out['C4v']['chi'].tolist()},
         'C4_deriv': {'dnu': 0.02, 'nu': out['C4d']['nu'].tolist(), 'chi': out['C4d']['chi'].tolist()},
         'C5': {'dnu': 0.04, 'nu': out['C5']['nu'].tolist(), 'c_re': out['C5']['c'].real.tolist(), 'c_im': out['C5']['c'].imag.tolist()}}
    json.dump(J, open(os.path.join(HERE, 'functionals_q6_canonical.json'), 'w'))
    pr(f"  wrote functionals_q6_canonical.json  [{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage6_canon.txt'), 'w').write(buf.getvalue())
