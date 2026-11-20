"""
sm_offset (SMOOTH_BOUND.md §3): the 'smooth offset' of the large-J near-gap rows (EE_DG1000_LARGEJ §3a: +200...+500 c2 at J = 48).
Per row (K2 or direct-state, EE1000_DSTATE) and key (C2/C3H, ++/+-): least-squares fit of Gamma - kappa-image (sm_curvfix, the
leading-Regge curvature in partial-wave phase) on the oscillating class {x^2k D_h, x^2k x D_h'/(2J+1)} (k <= 1) plus a constant r,
on x windows [XLO, 0.45] and [XLO, 0.7]; r with its least-squares error (residual-scaled).  Control: the same fit on the EXACT scalar
Gamma_s = m^2 (ex_2/E_2 - P_J) at the same nodes (sm_exact record with SM_XPROBE, J = 48 only), which has no smooth part.
Also printed: the zero of the heavy reference's LO shape, x0 = j_{0,1} m_ref/((2J+1) m) (the K2 reference's cumulant ratio is singular
there; the Hampel filter drops the node next to it: J = 40 x = 0.575, J = 48 x = 0.475).
Usage: W2_DMIN=1000 python3 sm_offset.py "J:D,J:D" [XLO]   (EE1000_DSTATE=1 for the direct-state rows)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ['SM_CURV'] = 'none'
import numpy as np
import ee1000_chain as CH
import sm_curvfix as CF
Q = CH.Q
XLO = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
DREF = 2e4


def fitr(xs, y, h, J, sel, K=1, refpole=None):
    D = CF.Dsh(h, J, xs); hh = 1e-6
    Dp = (CF.Dsh(h, J, xs + hh) - CF.Dsh(h, J, xs - hh))/(2*hh)*xs/(2*J + 1)
    Pw = [xs**(2*k) for k in range(K + 1)]
    extra = []
    if refpole is not None:                 # c(x) with the reference's 1/D_ref(x_ref) structure (oscillating class, non-polynomial c)
        Dr = CF.Dsh(h, J, xs*refpole)
        extra = [D/Dr, Dp/Dr]
    A = np.array([p*D for p in Pw] + [p*Dp for p in Pw] + extra + [np.ones_like(xs)]).T[sel]
    c, *_ = np.linalg.lstsq(A, y[sel], rcond=None)
    res = y[sel] - A @ c
    dof = max(len(res) - A.shape[1], 1)
    cov = np.linalg.inv(A.T @ A)*(res @ res)/dof
    return c[-1], np.sqrt(cov[-1, -1]), np.sqrt(np.mean(res**2)), np.sqrt(np.mean(y[sel]**2))


if __name__ == '__main__':
    tag = 'direct' if os.environ.get('EE1000_DSTATE') == '1' else 'K2'
    for spec in sys.argv[1].split(','):
        J, D0 = int(spec.split(':')[0]), float(spec.split(':')[1])
        g = Q.even_grid(J)
        if D0 not in g:
            print(f"J={J} D={D0}: no row ({sorted(g)})"); continue
        r = g[D0]; xs, m = r['xs'], r['m']
        mr = np.sqrt((DREF - J - 2)*(DREF + J - 1))
        x0 = 2.404826*mr/((2*J + 1)*m)
        print(f"J={J} Delta={D0:.0f} [{tag}] m={m:.1f} J^2/m={J*J/m:.2f} nodes {len(xs)}; reference LO zero at x0 = {x0:.3f}; "
              f"dropped near x0: {[x for x in np.arange(0.05, 0.81, 0.025) if abs(x - x0) < 0.02 and np.min(np.abs(xs - x)) > 1e-6]}")
        for kd in ('C2', 'C3H'):
            for h in ('++', '+-'):
                y = np.asarray(r[(kd, h)], float) - CF.curv_image(kd, h, J, m, xs, Q.A.env, Q.A.SIG)
                out = []
                for xh in (0.45, 0.7):
                    sel = (xs >= XLO - 1e-9) & (xs <= xh + 1e-9)
                    if sel.sum() < 8:
                        continue
                    rr, er, res, rms = fitr(xs, y, h, J, sel)
                    out.append(f"x in [{XLO},{xh}] (n={sel.sum()}): r = {rr:+6.0f} +- {er:4.0f} (resid {res:5.0f} of {rms:5.0f})")
                    if xh < x0 - 0.01:
                        rr, er, res, rms = fitr(xs, y, h, J, sel, refpole=m/mr)
                        out.append(f"+1/D_ref: r = {rr:+6.0f} +- {er:4.0f} (resid {res:5.0f})")
                print(f"   {kd:3s} {h}: " + "   ".join(out))
        # exact scalar control at the same nodes (if an sm_exact record with x-probes exists)
        pk = os.path.join(HERE, 'cache', 'sm_exact_1000.pkl')
        db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
        recs = [v for k, v in db.items() if k[0] == J and len(k) > 3]
        if recs:
            from scipy.special import j0
            import remk
            rec = recs[0]; n = rec['npn']; xp = rec['x'][n + 19:]; mm = rec['st']['m']
            b = 2*np.arcsinh((J + .5)/mm)
            P = remk.Pders_vec(xp*mm, np.cosh(b), nmax=0)[0]
            E2 = (2 - xp**2)/(1 - xp**2)**2
            Gs = mm*mm*(rec['rho'][2][n + 19:] + 2*(P - j0((2*J + 1)*xp))/E2)
            y = Gs - CF.curv_image('C2', '++', J, mm, xp, Q.A.env, Q.A.SIG)
            for xh in (0.45, 0.7):
                sel = (xp >= XLO - 1e-9) & (xp <= xh + 1e-9)
                rr, er, res, rms = fitr(xp, y, '++', J, sel)
                print(f"   exact scalar control (Delta={rec['st']['Delta']:.0f}, ++, x-probes, n={sel.sum()}): x <= {xh}: r = {rr:+6.1f} +- {er:4.1f} "
                      f"(resid {res:5.1f} of {rms:5.0f})")
