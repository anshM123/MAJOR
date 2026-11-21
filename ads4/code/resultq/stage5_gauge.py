"""Result Q, Stage 5 (RESULT_Q_STAGE5.md §5): gauge (null-space) dependence of the spinning 3-shift de-shift weight c(nu) of C5ds, and the
saved functional file.  All gauges below have the SAME leading-Regge heavy image sigma_h P_{-1/2}(cosh beta) and the same light-side moments;
they differ only beyond leading Regge.  Gauges (weighted min-norm, || D Y ||_2 minimal subject to the Abel equation, Y = (Y_R, Y_I)):
  G0: D = 1 (min-norm, the reference);  G1: D = (1 + (nu/5)^2) (penalise large nu);  G2: Im block weighted x 10 (as real as possible);
  G3: Re block weighted x 10 (as imaginary as possible).
Reports int|c|, the Im fraction int|Im c|/int|c|, support, and int|c - c_G0| / int|c_G0|.
Also writes functionals_q5_deshift.json: exact R_j (k = 5, ++ class; +- = -(++)), exact k = 3 R_0, raw-image constants k1, k2, k3, and the
G0 weight c(nu) (Re, Im) on the dnu = 0.02, N = 2048 grid (functional: int Re[c(nu) Psi_{5,nu}/Nhat5(nu)] dnu, Nhat5 = sum_j R_j).
Output: cache/stage5_gauge.txt, functionals_q5_deshift.json
"""
import os, sys, io, json, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import stage5_hj as S5
import stage5_deshift as DS

buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def weighted_solve(S, d):
    """min ||d * Y|| s.t. M Y = b  ->  Z = d Y, min-norm for M diag(1/d) Z = b."""
    b = np.ones_like(S.t)/S.ch2
    Md = S.M/d[None, :]
    U, s, Vh = np.linalg.svd(Md, full_matrices=False)
    keep = s > 1e-11*s[0]
    Z = Vh[keep].T @ ((U[:, keep].T @ b)/s[keep])
    Y = Z/d
    res = np.linalg.norm(S.M @ Y - b)/np.linalg.norm(b)
    return (Y[:S.N] + 1j*Y[S.N:])/S.n, res


if __name__ == '__main__':
    t0 = time.time()
    S = DS.Deshift3S('C5_pp', 0.04, 1024, general=True)
    one = np.ones(S.N)
    G = {'G0 min-norm': np.concatenate([one, one]),
         'G1 (1+(nu/5)^2)': np.concatenate([1 + (S.nu/5)**2]*2),
         'G2 Im x10': np.concatenate([one, 10*one]),
         'G3 Re x10': np.concatenate([10*one, one])}
    sol = {}
    pr("# stage5_gauge: spinning C5 de-shift, dnu = 0.04, N = 1024 (general complex gauge, t oversampled x2)")
    for k, d in G.items():
        c, res = weighted_solve(S, d)
        sol[k] = c
        W = S.W(c, np.array([0., 2., 4., 6.])); mom = S.moments(c)
        pr(f"  {k:17s}: residual {res:.1e}; W(0,2,4,6) = " + " ".join(f"{v:+.9f}" for v in W) + f"; moments " + ", ".join(f"{v:+.1e}" for v in mom) +
           f"; int|c| = {S.dnu*np.abs(c).sum():.4f}; int|Im c|/int|c| = {np.abs(c.imag).sum()/np.abs(c).sum():.3f}; "
           f"support (|c| > 1e-4 max) nu <= {S.nu[np.abs(c) > 1e-4*np.abs(c).max()].max():.1f}")
    c0 = sol['G0 min-norm']
    for k, c in sol.items():
        pr(f"  {k:17s}: int|c - c_G0|/int|c_G0| = {np.abs(c - c0).sum()/np.abs(c0).sum():.3f}")
    # saved functional (G0, dnu = 0.02, N = 2048 from stage5_deshift)
    P = pickle.load(open(os.path.join(HERE, 'cache', 'stage5_deshift_C5_pp.pkl'), 'rb'))
    D = S5.load()
    out = {'description': "Result Q Stage 5: spinning (j'=1, o_H q) k=5 de-shift. Raw leading-Regge image of Psi_{5,nu} (++ class) = "
                          "N(nu) sum_{j in {-2,0,2}} R_j(nu) P_{-1/2+i(nu-ij)}(cosh beta); +- class = -(++) (sigma_h). Functional "
                          "C5ds = int_0^nu_max Re[c(nu) Psi_{5,nu}/Nhat5(nu)] dnu, Nhat5 = sum_j R_j; heavy image sigma_h P_{-1/2}(cosh beta); "
                          "min-norm gauge (stage5_deshift.py general gauge, t-grid oversampled x2).",
           'R_j_k5_pp': {str(j): str(v) for j, v in D['C5_pp'].items()},
           'R_0_k3_pp': str(D['C3_pp'][0]),
           'raw_image_nu0': {'k1_dmuP_3/2': '-192491/1028120', 'k2_P_3/2': '1252043/1028120', 'k3_P_-1/2': '-223923/1028120'},
           'dnu': float(P['dnu']), 'nu': P['nu'].tolist(), 'c_re': P['c'].real.tolist(), 'c_im': P['c'].imag.tolist()}
    json.dump(out, open(os.path.join(HERE, 'functionals_q5_deshift.json'), 'w'))
    pr(f"  wrote functionals_q5_deshift.json (N = {len(P['nu'])})  [{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'stage5_gauge.txt'), 'w').write(buf.getvalue())
