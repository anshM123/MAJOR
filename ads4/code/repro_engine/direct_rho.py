"""repro_engine/direct_rho.py -- DIRECT evaluation (no 1/Delta expansion) of the EE k = 4 image I(nu; Delta) of the J = 0 state on the
real section D: seven-dimensional tensor-product PRODUCT-INTEGRATION rule (quad1d.py) applied to the exact point-wise integrand
(integrand.py).  The nu-dependence  A^{i nu/2} = 4^{i nu/2} prod atoms^{i nu E1} (1 + eps)  is split into the exactly-integrated atom
powers (in the 1-D weights) and the tiny remainder eps (|eps| < 1e-6), kept to second order.

Symmetries used (both exact, checked point-wise): s -> -s averaging (t = s^2 variable), and the transverse reflection
(s, phi, alpha) -> (-s, -phi, -alpha) (alpha restricted to (0, pi)).  With them the integrand divided by its real weight
X1^2.25 X2^1.75 Y1^0.25 Y2^-0.25 e^{-S} w^1.75 cos(phi)^0.5 sin(alpha)^-0.5 is analytic in X1, X2, Y1, Y2, t, phi, alpha.

Usage: python3 direct_rho.py Delta nX1,nX2,nY1,nY2,nt,nphi,nal [f64|f80] [numax]
       -> cache/direct_D{Delta}_n{...}_{prec}.npz  (nus, S0, S1, S2, shift)
"""
import os, sys, time
import numpy as np
import integrand as IG
import quad1d as Q

HERE = os.path.dirname(os.path.abspath(__file__))
DCOEF = float(os.environ['DCOEF']) if os.environ.get('DCOEF') else None     # order-0 diagnostics: <TTO> coefficients at DCOEF
NUS = np.round((np.arange(1025) + 0.5)*0.04, 4)

# (kind, real exponent, E1 of the nu-atom) per variable;  variables: X1, X2, Y1, Y2 (X-form), t = s^2, phi, alpha
SPEC = [('lag', 2.25, -0.5), ('lag', 1.75, +0.5), ('lag', 0.25, +0.5), ('lag', -0.25, -0.5),
        ('tjac', 1.75, -0.5), ('cosw', 0.5, +1.0), ('sinw', -0.5, -1.0)]


def rules(ns):
    return [Q.Rule(k, n, p, e) for (k, p, e), n in zip(SPEC, ns)]


def run(Delta, ns, prec='f64', numax=12.0, verbose=True):
    typ = np.float64 if prec == 'f64' else np.longdouble
    nus = NUS[NUS < numax]
    R = rules(ns)
    W = [r.weights(nus) for r in R]                     # (n_d, nnu) complex
    X1, X2, Y1, Y2 = (R[i].x for i in range(4))
    g = np.meshgrid(X1, X2, Y1, Y2, indexing='ij')
    g = [v.ravel() for v in g]
    shp = tuple(ns[:4])
    S = [np.zeros(len(nus), complex) for _ in range(3)]
    shift = None
    t0 = time.time()
    logW_rad = 2.25*np.log(g[0]) - g[0] + 1.75*np.log(g[1]) - g[1] + 0.25*np.log(g[2]) - g[2] - 0.25*np.log(g[3]) - g[3]
    lat_rad = -0.5*np.log(g[0]) + 0.5*np.log(g[1]) + 0.5*np.log(g[2]) - 0.5*np.log(g[3])      # sum E1 log atom (radial)
    epsmax = 0.0
    nb = 0
    for it, t in enumerate(R[4].x):
        w = 1.0 - t
        for ip, ph in enumerate(R[5].x):
            cf = np.cos(ph)
            for ia, al in enumerate(R[6].x):
                sa = np.sin(al)
                logW = logW_rad + 1.75*np.log(w) + 0.5*np.log(cf) - 0.5*np.log(sa)
                lat = lat_rad - 0.5*np.log(w) + np.log(cf) - np.log(sa)
                Hs = [0.0, 0.0, 0.0]
                for sg in (1.0, -1.0):
                    args = [v.astype(typ) for v in g] + [np.full(g[0].shape, sg*np.sqrt(t), dtype=typ),
                                                         np.full(g[0].shape, ph, dtype=typ), np.full(g[0].shape, al, dtype=typ)]
                    if prec != 'f64':
                        # node values in extended precision (sqrt(t), angles re-derived from the float64 nodes)
                        pass
                    P = IG.point(*args, typ(Delta), Dcoef=DCOEF)
                    lh = np.asarray(P['lG'] - logW, dtype=np.float64)
                    if shift is None:
                        shift = float(np.max(lh))
                    h = np.asarray(P['sG'], dtype=np.float64)*np.exp(lh - shift)
                    eps = np.asarray(P['lA'] - 2*lat - np.log(typ(4.0)), dtype=np.float64)
                    epsmax = max(epsmax, float(np.max(np.abs(eps))))
                    Hs[0] = Hs[0] + h; Hs[1] = Hs[1] + h*eps; Hs[2] = Hs[2] + h*eps*eps
                ang = W[4][it]*W[5][ip]*W[6][ia]                                  # (nnu,)
                for m in range(3):
                    Hm = Hs[m].reshape(shp)
                    T = np.tensordot(Hm, W[0], axes=([0], [0]))                   # (n2, n3, n4, nnu)
                    T = np.einsum('bcdv,bv->cdv', T, W[1])
                    T = np.einsum('cdv,cv->dv', T, W[2])
                    T = np.einsum('dv,dv->v', T, W[3])
                    S[m] += ang*T
                nb += 1
        if verbose:
            print(f'   D={Delta:g} {prec} t-node {it + 1}/{ns[4]}: {time.time() - t0:.0f}s, max|eps| {epsmax:.1e}', flush=True)
    fn = os.path.join(HERE, 'cache', 'direct_D%g_n%s_%s%s.npz' % (Delta, '-'.join(map(str, ns)), prec,
                                                                  '' if DCOEF is None else '_coef%g' % DCOEF))
    np.savez(fn, nus=nus, S0=S[0], S1=S[1], S2=S[2], shift=shift, epsmax=epsmax, ns=np.array(ns), Delta=Delta)
    return fn


def load(fn):
    d = np.load(fn)
    nus = d['nus']
    z = 0.5j*nus
    I = d['S0'] + z*d['S1'] + 0.5*z*z*d['S2']
    return nus, I, d


if __name__ == '__main__':
    Delta = float(sys.argv[1]); ns = [int(v) for v in sys.argv[2].split(',')]
    prec = sys.argv[3] if len(sys.argv) > 3 else 'f64'
    numax = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0
    t0 = time.time()
    fn = run(Delta, ns, prec, numax)
    print(f'direct_rho D={Delta:g} n={ns} {prec}: {time.time() - t0:.0f}s -> {fn}', flush=True)
