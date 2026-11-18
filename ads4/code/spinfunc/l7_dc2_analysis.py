"""
l7 analysis: from the l7_dc2_scan pickles, per even J
 (i)  LO check (item 3): at order 0, the 2x2 action over the x grid lies in a 2D span of two rank-1 matrices R_++ = a a^T, R_+- = b b^T
      (third singular value), and in the helicity basis T = [a b]:  d_++/Y_scal ~ (1-x^2)^2,  d_+-/d_++ ~ d^J_44/((1-x^2)^4 P_J)
      (stripped flat densities, doc 9/17).  Helicities are identified by these shapes (no external input).
 (ii) delta_h(x; m) = d_h(x) / (C_h W_h(x) Y_scal(x)) - 1, W_++ = (1-x^2)^2, W_+- = (1-x^2)^{-2} d44/P_J, C_h fixed by delta_h(x_ref) = 0:
      the graviton-specific relative correction (graviton action vs the validated scalar k = 2 action, same nu, Delta, J), up to a per-state,
      per-helicity constant.  Fitted across Delta at fixed physical x: delta = e1/m + e2/m^2 + e3/m^3 + e4/m^4  (dc2 := e2 + ...).
 (iii) helicity mixing: M~_{+,-}/sqrt(M~_{++} M~_{+-}) at order 2 (T from the order-0 decomposition at the same Delta).
Usage: python3 l7_dc2_analysis.py [xref]  -> cache/l7_dc2_analysis.txt, cache/l7_dc2.pkl
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import eval_jacobi, eval_legendre


def d44(J, c):
    return ((1 + c)/2)**4*eval_jacobi(J - 4, 0, 8, c)


def W(h, J, x):
    c = 1 - 2*x*x
    return (1 - x*x)**2 if h == '++' else d44(J, c)/((1 - x*x)**2*eval_legendre(J, c))


def Dh(h, J, x):
    c = 1 - 2*x*x
    return eval_legendre(J, c) if h == '++' else d44(J, c)


def vec(M):
    return np.array([M[0, 0], M[0, 1], M[1, 1]])


def rank1(Ms):
    V = np.array([vec(M) for M in Ms])
    _, s, Vh = np.linalg.svd(V)
    mat = lambda v: np.array([[v[0], v[1]], [v[1], v[2]]])
    A1, A2 = mat(Vh[0].conj()), mat(Vh[1].conj())
    a = np.linalg.det(A2); b = A1[0, 0]*A2[1, 1] + A1[1, 1]*A2[0, 0] - 2*A1[0, 1]*A2[0, 1]; c = np.linalg.det(A1)
    R = [A1 + t*A2 for t in np.roots([a, b, c])]
    vs = []
    for Ri in R:                                           # R = v v^T
        i = int(np.argmax(np.abs(np.diag(Ri)))); v = Ri[:, i]/np.sqrt(Ri[i, i]); vs.append(v)
    return s/s[0], vs


def analyse(J, res, xref):
    Dls = sorted(set(k[0] for k in res)); xs = np.array(sorted(set(k[1] for k in res)))
    out = {}
    for Dl in Dls:
        if not all((Dl, x) in res for x in xs):
            continue
        rr = [res[(Dl, x)] for x in xs]
        m = rr[0]['m']
        pref = np.array([[np.exp(r['logC'] - r['logK'] - r['sl'])/r['sS'][o] for o in range(3)] for r in rr])   # (nx, 3)
        if J >= 4:
            sv, vs = rank1([r['M'][:, :, 0] for r in rr])
            T = np.array(vs).T                                         # columns a, b
            Ti = np.linalg.inv(T)
            Mt = np.array([[Ti @ r['M'][:, :, o] @ Ti.T for o in range(3)] for r in rr])   # (nx, 3, 2, 2)
            d = np.array([[Mt[i, o, 0, 0], Mt[i, o, 1, 1]] for i in range(len(xs)) for o in range(3)]).reshape(len(xs), 3, 2)
            # identify ++: d/Y_scal / (1-x^2)^2 flattest at order 0
            flat = []
            for j in range(2):
                q = d[:, 0, j]*pref[:, 0]/W('++', J, xs)
                flat.append(np.std(np.abs(q/q[np.argmin(abs(xs - xref))])))
            ipp = int(np.argmin(flat)); ipm = 1 - ipp
            hel = {'++': d[:, :, ipp], '+-': d[:, :, ipm]}
            mix = np.array([Mt[i, 2, 0, 1]/np.sqrt(Mt[i, 2, 0, 0]*Mt[i, 2, 1, 1]) for i in range(len(xs))])
        else:
            sv = None; hel = {'++': np.array([[r['M'][0, 0, o] for o in range(3)] for r in rr])}; mix = None
        iref = int(np.argmin(abs(xs - xref)))
        rec = dict(m=m, sv=sv, mix=mix, J=J)
        for h, dd in hel.items():
            q = dd*pref/W(h, J, xs)[:, None]                            # (nx, 3)
            rec[h] = q/q[iref][None, :]                                 # normalised: 1 + delta_h
        out[Dl] = rec
    return xs, out


def fit_dc2(xs, out, h, nterm=4):
    """log(1+delta) at order 2, fitted across Delta at fixed x: sum_k e_k(x)/m^k; then per k over x:
    e_k(x) = alpha_k x^-k [nu-only factor, removable by the functional's nu-only normalisation N^_k(nu)] + b_k + g_k x^2 + d_k x^4.
    Regular part r_k(x) = g_k x^2 + d_k x^4 (relative to x = 0; the per-state constant b_k is irrelevant).  Points with |D_h| < 0.1 get
    weight 0.1 (near zeros of the leading density the relative correction is ill-defined)."""
    Dls = sorted(out); ms = np.array([out[D]['m'] for D in Dls])
    L = np.array([np.log(out[D][h][:, 2]) for D in Dls])               # (nD, nx)
    A = np.array([ms**(-k) for k in range(1, nterm + 1)]).T
    e = np.linalg.lstsq(A, L, rcond=None)[0]                            # (nterm, nx)
    A0 = np.array([ms**(-k) for k in range(0, nterm + 1)]).T
    e0 = np.linalg.lstsq(A0, L, rcond=None)[0][0]
    J = out[Dls[0]]['J']
    w = np.maximum(np.minimum(np.abs(Dh(h, J, xs))/0.1, 1.0), 0.1)
    reg = []
    for k in range(1, 3):
        B = np.array([xs**(-k), np.ones_like(xs), xs**2, xs**4]).T
        cf = np.linalg.lstsq(B*w[:, None], e[k - 1]*w, rcond=None)[0]
        resid = np.abs(B @ cf - e[k - 1])
        reg.append(dict(alpha=cf[0], b=cf[1], g=cf[2], d=cf[3], r=cf[2]*xs**2 + cf[3]*xs**4, resid=resid))
    return ms, e, e0, reg


def main(xref=0.05):
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr(f"# l7: graviton C^_2 vs scalar k = 2 crystal actions; x_ref = {xref} (delta_h(x_ref) := 0 per state and helicity)")
    store = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'l7_J*.pkl')), key=lambda f: int(f.split('_J')[1][:-4])):
        J = int(fn.split('_J')[1][:-4])
        res = pickle.load(open(fn, 'rb'))
        xs, out = analyse(J, res, xref)
        if len(out) < 4:
            pr(f"  J={J}: only {len(out)} complete Delta rows, skipped"); continue
        Dls = sorted(out)
        pr(f"  J={J}: Delta = {Dls}, x = {list(xs)}")
        if J >= 4:
            pr("    LO (order 0) 2D-span test, 3rd singular value / 1st: " + "  ".join(f"{D:.0f}: {out[D]['sv'][2]:.1e}" for D in Dls))
        for h in (('++', '+-') if J >= 4 else ('++',)):
            lo = np.array([np.abs(out[D][h][:, 0] - 1).max() for D in Dls])
            pr(f"    {h}: LO shape test max_x |q_0(x)/q_0(x_ref) - 1| vs Delta: " + " ".join(f"{v:.1e}" for v in lo))
            ms, e, e0, reg = fit_dc2(xs, out, h)
            ph = np.array([np.abs(np.angle(out[D][h][:, 2])).max() for D in Dls])
            pr(f"       order-2 phase max_x |arg(1+delta)| vs Delta: " + " ".join(f"{v:.1e}" for v in ph) + f";  LO limit |e0| max {np.abs(e0).max():.1e}")
            for k in (1, 2):
                R = reg[k - 1]
                pr(f"       e{k}(x) = alpha x^-{k} + b + g x^2 + d x^4: alpha = {R['alpha']:.3f}, regular r{k}(x) = g x^2 + d x^4 with "
                   f"g = {R['g']:.2f}, d = {R['d']:.2f}; max fit resid {R['resid'].max():.2f}")
            r1, r2 = reg[0]['r'], reg[1]['r']
            dc2 = (r2 + 0.5*r1**2).real
            pr("       regular dc2(x) [= r2 + r1^2/2, rel. to x = 0]: " + "  ".join(f"{x:.2f}:{v:+.1f}" for x, v in zip(xs, dc2)))
            pr("       regular r1(x) [genuine x-dependent 1/m]:          " + "  ".join(f"{x:.2f}:{v.real:+.2f}{v.imag:+.2f}i" for x, v in zip(xs, r1)))
            store[(J, h)] = dict(xs=xs, e=e, reg=reg, dc2=dc2, r1=r1, ms=ms)
        if J >= 4:
            mixs = np.array([np.abs(out[D]['mix']).max() for D in Dls])
            pr("    helicity mixing max_x |M~_{+,-}|/sqrt|M~_++ M~_+-| (order 2) vs Delta: " + " ".join(f"{v:.1e}" for v in mixs)
               + "  ;  times m: " + " ".join(f"{v*out[D]['m']:.2f}" for v, D in zip(mixs, Dls)))
    pickle.dump(store, open(os.path.join(HERE, 'cache', 'l7_dc2_xref%.2f.pkl' % xref), 'wb'))
    open(os.path.join(HERE, 'cache', 'l7_dc2_analysis_xref%.2f.txt' % xref), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.05)
