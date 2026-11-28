"""Flip-channel bound:  c * sum_n kappa_n rhs_n = sum_s xi_s G(s),  |xi| <= sqrt(rho0 rho4),
|G(s)| <= 2 t sqrt(W0 W4) on heavy states  =>  |c| <= t * sum_s(rho0 W0 + rho4 W4) <= 3 t  (endpoint, D=0, + light share).
Light states (2<twist<g, parity-even even J>=4): coefficient  ell = max |G|/(2 sqrt(W0 W4)) multiplies the light endpoint share.
LP: minimise t over kappa, with sum kappa rhs = 1 and optional cap ell <= Theta."""
import sys, time, numpy as np
from scipy.optimize import linprog
from flip import flip_weights
from vec import Wv


def grid_states(g, Jlin, Jmax, nJlog, ts):
    Js = np.array(sorted(set(list(range(4, Jlin + 1, 2)) + [2 * (int(x) // 2) for x in np.geomspace(Jlin + 2, Jmax, nJlog)])), float)
    JJ, DD = np.meshgrid(Js, ts, indexing='ij')
    return JJ.ravel(), (JJ + DD).ravel()


def solve(g, N, theta=None, Jlin=None, Jmax=None, nJlog=150, nt=60):
    Jlin = Jlin or 12 * g; Jmax = Jmax or 5000 * g
    ts = np.concatenate([[0.0], np.geomspace(1e-3, 1e5, nt)]) * 1.0 + g
    J, D = grid_states(g, Jlin, Jmax, nJlog, ts)
    Fw, rhs, lab = flip_weights(J, D, N)
    sw = 2 * np.sqrt(Wv(np.zeros_like(J), D) * Wv(4 * np.ones_like(J), D))
    A = Fw / sw[:, None]
    rows = [A, -A]; b = [np.zeros(len(J)), np.zeros(len(J))]
    tcol = [-np.ones(len(J)), -np.ones(len(J))]
    if theta is not None:
        lh = np.concatenate([[2.0 + 1e-9, 2.01, 2.05, 2.2, 2.5, 3, 4, 5, 6, 8], np.geomspace(10, g - 1e-6, 40)])
        Jl, Dl = grid_states(g, Jlin, Jmax, nJlog, lh)
        Fl, _, _ = flip_weights(Jl, Dl, N)
        swl = 2 * np.sqrt(Wv(np.zeros_like(Jl), Dl) * Wv(4 * np.ones_like(Jl), Dl))
        Al = Fl / swl[:, None]
        rows += [Al, -Al]; b += [theta * np.ones(len(Jl))] * 2; tcol += [np.zeros(len(Jl))] * 2
    Aub = np.hstack([np.vstack(rows), np.concatenate(tcol)[:, None]])
    bub = np.concatenate(b)
    nv = Aub.shape[1]
    sc = np.abs(Aub).max(axis=0); sc[sc == 0] = 1
    Aeq = np.concatenate([rhs, [0.0]])[None, :]
    c = np.zeros(nv); c[-1] = 1.0
    res = linprog(c / sc, A_ub=Aub / sc, b_ub=bub, A_eq=Aeq / sc, b_eq=[1.0], bounds=[(None, None)] * (nv - 1) + [(0, None)], method='highs')
    if res.status != 0:
        return None, res
    x = res.x / sc
    return dict(t=x[-1], kappa=x[:-1], lab=lab), res


if __name__ == '__main__':
    g = int(sys.argv[1]); th = None if len(sys.argv) < 4 else float(eval(sys.argv[3]))
    for N in [int(x) for x in sys.argv[2].split(',')]:
        t0 = time.time()
        o, res = solve(g, N, th)
        if o is None: print(g, N, 'fail', res.message); continue
        tau = 12 * o['t']
        print(f"g={g} N={N} theta={th}: t={o['t']:.4e}  tau<=12t={tau:.4e}  tau*g={tau*g:.3f} tau*g^2={tau*g*g:.2f}  [{time.time()-t0:.0f}s]", flush=True)
