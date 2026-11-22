"""
lsd_hairpin (LIGHTSIDE_DIST.md): EIN's whole-function swapped form on a contour family C_tau that interpolates between the circle C
(tau = 0, EIN's E+) and a contour that leaves the endpoint t = 0 TANGENTIALLY to the double-commutator cut (tau = 1, the hairpin
boundary values t -+ i0 near the Regge corner), with everything else unchanged.  Read-only imports: ein_dbar, ein_run, og_free.

    t(theta) = (1 - cos theta)/2 - (i/2) kappa(theta) sin theta,   kappa(theta) = (1 - tau) + tau sin^2(theta/2)
    (kappa = 1 at theta = pi for every tau; near theta -> 0: Im t ~ -theta^3/8 vs Re t ~ theta^2/4 when tau = 1).

The sigma part is continued along the physical route (real spacelike axis t = inf -> 1, then along C_tau), exactly as EIN's
continue_ellipse, with mpmath jets at every node (EIN §2.5 pitfall 2).

    (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_hairpin.py <tag> <k> <tau> <n>)   -> cache/lsd_hp_<tag>_k<k>_tau<tau>_n<n>.pkl
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ein_dbar as E
import og_free as F
import lsd_pos as P

CLS = ('1234', '1243', '1324', 's', 't', 'u')


def contour_tau(n, tau):
    x, wx = np.polynomial.legendre.leggauss(n)
    psi = (x + 1)*np.pi/2; wpsi = wx*np.pi/2
    th = np.pi*(1 - np.cos(psi)); dth = np.pi*np.sin(psi)
    kap = (1 - tau) + tau*np.sin(th/2)**2
    dkap = tau*np.sin(th/2)*np.cos(th/2)
    t = (1 - np.cos(th))/2 - 0.5j*kap*np.sin(th)
    dt = (np.sin(th)/2 - 0.5j*(dkap*np.sin(th) + kap*np.cos(th)))*dth
    return t, wpsi*dt, th


def _tt(th, tau):
    kap = (1 - tau) + tau*np.sin(th/2)**2
    return (1 - np.cos(th))/2 - 0.5j*kap*np.sin(th)


def continue_tau(th1, th3, y, tau, eps=1e-3, n_real=600, n_line=3000, bump=0.05):
    """physical route of ein_dbar.continue_ellipse, with the contour C_tau."""
    th1, th3 = np.broadcast_arrays(np.asarray(th1, float), np.asarray(th3, float))
    U0, V0 = E.UV_of(eps, eps, y)
    st = E.principal_state(np.full(th1.shape, U0, complex), np.full(th1.shape, V0, complex))
    one = np.ones(th1.shape, complex)
    for s in np.linspace(0, 1, n_real + 1)[1:]:
        r = eps + (1 - eps)*s + 1j*bump*np.sin(np.pi*s)
        st = E.step(st, *E.UV_of(r*one, r*one, y))
    P1 = E._theta_path(th1, n_line); P3 = E._theta_path(th3, n_line)
    for th in P1:
        st = E.step(st, *E.UV_of(1/_tt(th, tau), one, y))
    for th in P3:
        st = E.step(st, *E.UV_of(1/_tt(th1, tau), 1/_tt(th, tau), y))
    return st


def sigma_grid_tau(n, y, tau, pieces, chunk=10):
    t, w, th = contour_tau(n, tau)
    out = {p: np.empty((n, n), complex) for p in pieces}
    for i0 in range(0, n, chunk):
        st = continue_tau(th[i0:i0 + chunk, None], th[None, :], y, tau)
        T1, T3 = np.broadcast_arrays(t[i0:i0 + chunk, None], t[None, :])
        D, _ = E.phi_derivs_mixed(st, T1, T3, y, 6)
        Xs = E.Xgauge(T1, T3, y)
        for p in pieces:
            out[p][i0:i0 + chunk] = E.evaluate(pieces[p], Xs, D)
    return out


def Etable_tau(ks, phi1, phi3, p2, p4, n, tau, S=None):
    y = 2*np.sin((phi1 - phi3)/2)**2
    t, w, th = contour_tau(n, tau)
    PC = E.pieces()
    if S is None:
        S = sigma_grid_tau(n, y, tau, PC)
    z1 = F.zvec(phi1); z3 = F.zvec(phi3)
    G = F.correlator(F.links(t[:, None], t[None, :], z1, z3, p2, p4, y13=y), CLS)
    T1b, T3b = np.broadcast_arrays(t[:, None], t[None, :])
    out = {}
    for k in ks:
        h = T1b**5*T3b**5*F.fk(k, T1b, T3b, z1, z3, y13=y)*(w[:, None]*w[None, :])
        out[k] = {(p, c): complex(np.sum(h*S[p]*G[c])) for p in S for c in CLS}
    return out, S


if __name__ == '__main__':
    tag, n = sys.argv[1], int(sys.argv[4])
    ks = tuple(int(x) for x in sys.argv[2].split(','))
    tau = float(sys.argv[3])
    a, b, p2, p4 = P.config(tag)
    t0 = time.time()
    res = {}
    S = None
    # sigma grid depends on y only (same for (a,b) and (b,a))
    res['ab'], S = Etable_tau(ks, a, b, p2, p4, n, tau)
    res['ba'], _ = Etable_tau(ks, b, a, p2, p4, n, tau, S=S)
    print(f'{tag} k={ks} tau={tau} n={n}: {time.time() - t0:.0f}s', flush=True)
    for k in ks:
        sg = 1 if k % 2 == 0 else -1
        for p in ('c0', 'c2', 'c3', 'c4', 'xs', 'xt', 'xu'):
            A = sum(res['ab'][k][(p, c)] for c in CLS); B = sum(res['ba'][k][(p, c)] for c in CLS)
            print(f'  k={k} {p}: E(a,b) {A.real:+.6e}{A.imag:+.1e}i  E(b,a) {B.real:+.6e}  sym {(A + sg*B).real:+.4e}', flush=True)
    with open(os.path.join(HERE, 'cache', f'lsd_hp_{tag}_k{sys.argv[2].replace(",", "")}_tau{tau}_n{n}.pkl'), 'wb') as fh:
        pickle.dump(res, fh)
