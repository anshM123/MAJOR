"""
pcyc_dquad (PCYC_CLOSE.md section 3.1): the integrals over the PHYSICAL chain D' (= CLS-D causal domain, pcyc_region) of the scalar
and of the EE (k = 4, J = 0) CLS-D integrands, EXACT in Delta (no 1/Delta expansion, the cut at the light cones of x5 = e included),
by a seven-dimensional tensor-product PRODUCT-INTEGRATION rule in the coordinates of pcyc_dlib (own code):

  * t_x, t_y : generalised Gauss-Laguerre, weight t^p exp(-lam t)  (p = 6 scalar / 8 EE;  lam = 1 - 1/Delta / 1 - 5/Delta: the measured
               decay rates of the integrands toward the B zero, so that the remainders are smooth and bounded);
  * r_x, r_y : r^c1 (1-r)^c2 x [r^e1 (1-r)^e2]^{i nu}   (scalar (2.25, 1.75)/(2.25, 1.75); EE (2.25, 1.75)/(0.25, -0.25));
  * t = s^2  : (1-t)^1.75 t^-1/2 x (1-t)^{-i nu/2}, both signs of s evaluated (the s-even part is analytic in all variables);
  * phi      : cos^h x cos^{i nu}  (h = 4.5 / 0.5);   alpha in (0, pi): sin^-1/2 x sin^{-i nu} (transverse parity doubles it).
  The nu-dependence  A^{i nu/2} = 4^{i nu/2} prod atoms^{...} exp(i nu eps/2)  is split into the exactly integrated atom powers (in the
  1-D product weights, computed by fine tanh-sinh quadrature) and exp(i nu eps/2), evaluated EXACTLY at the nodes (no expansion).
  Nodes: Gauss nodes of the real weight; weights omega_k(nu) = int W atom^{i nu} l_k (interpolatory: exact for polynomial remainders).
Outputs I_S(nu), I_E(nu) (constants dropped: same constants on the crystal side cancel in <Q> = I_E/I_S).
Usage: python3 pcyc_dquad.py Delta kind n1,...,n7 [nus] -> cache/dq_D<Delta>_<kind>_<grid>.npz       (kind S or E)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import roots_genlaguerre, roots_jacobi, gammaln
import pcyc_dlib as P
np.seterr(all='ignore')

NUS_DEF = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0)


def bary(x):
    n = len(x); w = np.ones(n)
    for j in range(n):
        w[j] = 1.0/np.prod(x[j] - np.delete(x, j))
    return w


def lagmat(xn, t):
    wb = bary(xn)
    D = t[None, :] - xn[:, None]
    hit = np.abs(D) < 1e-300
    D = np.where(hit, 1.0, D)
    T = wb[:, None]/D
    L = T/np.sum(T, axis=0, keepdims=True)
    for i in np.where(hit.any(axis=0))[0]:
        L[:, i] = hit[:, i].astype(float)
    return L


def ts01(h=1/256., tmax=4.3):
    tau = np.arange(-tmax, tmax + h/2, h)
    zz = 0.5*np.pi*np.sinh(tau)
    t = 1/(1 + np.exp(-2*zz)); u = 1/(1 + np.exp(2*zz))
    w = h*np.pi*np.cosh(tau)*t*u
    k = (t > 0) & (u > 0) & (w > 0)
    return t[k], u[k], w[k]


class Rule:
    """1-D product rule.  kind: 'lag' (p, lam), 'jac' (c1, c2, e1, e2), 'tjac' (g, e), 'cos' (h, e), 'sin' (k, e)."""

    def __init__(self, kind, n, *par):
        self.kind, self.n, self.par = kind, n, par
        if kind == 'lag':
            p, lam = par
            x, w = roots_genlaguerre(n, p)
            self.x = x/lam; self.w0 = w/lam**(p + 1); self.real = True
            self.logW = lambda t: p*np.log(t) - lam*t
            return
        self.real = False
        t, u, w = ts01()
        if kind == 'jac':
            c1, c2, e1, e2 = par
            y, _ = roots_jacobi(n, c2, c1); x = (y + 1)/2
            fx = t; dw = w*t**c1*u**c2; la = e1*np.log(t) + e2*np.log(u)
            self.logW = lambda r: c1*np.log(r) + c2*np.log1p(-r)
        elif kind == 'tjac':
            g, e = par
            y, _ = roots_jacobi(n, g, -0.5); x = (y + 1)/2
            fx = t; dw = w*u**g/np.sqrt(t); la = e*np.log(u)
            self.logW = lambda tt: g*np.log1p(-tt)          # t^-1/2 = 2 ds/dt is the Jacobian of t = s^2, not part of the integrand
        elif kind == 'cos':
            h, e = par
            y, _ = roots_jacobi(n, h, h); x = y*np.pi/2
            fx = np.pi*(t - 0.5); c = np.sin(np.pi*np.minimum(t, u))
            dw = np.pi*w*c**h; la = e*np.log(c)
            self.logW = lambda ph: h*np.log(np.cos(ph))
        elif kind == 'sin':
            k, e = par
            y, _ = roots_jacobi(n, k, k); x = (y + 1)*np.pi/2
            fx = np.pi*t; c = np.sin(np.pi*np.minimum(t, u))
            dw = np.pi*w*c**k; la = e*np.log(c)
            self.logW = lambda al: k*np.log(np.sin(al))
        else:
            raise ValueError(kind)
        self.x = np.asarray(x, float)
        self._L = lagmat(self.x, fx)*dw[None, :]
        self._la = la

    def weights(self, nus):
        nus = np.asarray(nus, float)
        if self.real:
            return np.repeat(self.w0[:, None], len(nus), axis=1).astype(complex)
        return self._L @ np.exp(1j*np.outer(self._la, nus))


SPEC = {'S': [('lag', 6.0, None), ('jac', 2.25, 1.75, -0.5, 0.5), ('tjac', 1.75, -0.5), ('lag', 6.0, None),
              ('jac', 2.25, 1.75, 0.5, -0.5), ('cos', 4.5, 1.0), ('sin', -0.5, -1.0)],
        'E': [('lag', 8.0, None), ('jac', 2.25, 1.75, -0.5, 0.5), ('tjac', 1.75, -0.5), ('lag', 8.0, None),
              ('jac', 0.25, -0.25, 0.5, -0.5), ('cos', 0.5, 1.0), ('sin', -0.5, -1.0)]}
LAMC = {'S': 1.0, 'E': 5.0}


def rules(kind, ns, Delta):
    out = []
    for sp, n in zip(SPEC[kind], ns):
        if sp[0] == 'lag':
            out.append(Rule('lag', n, sp[1], 1.0 - LAMC[kind]/Delta))
        else:
            out.append(Rule(sp[0], n, *sp[1:]))
    return out


def run(Delta, kind, ns, nus=NUS_DEF, verbose=True):
    nus = np.asarray(nus, float)
    R = rules(kind, ns, Delta)
    W = [r.weights(nus) for r in R]
    tx, rx, ts, ty, ry, ph, al = (r.x for r in R)
    G = np.meshgrid(tx, rx, ty, ry, indexing='ij')
    G = [g.ravel() for g in G]
    shp = (ns[0], ns[1], ns[3], ns[4])
    lw4 = R[0].logW(G[0]) + R[1].logW(G[1]) + R[3].logW(G[2]) + R[4].logW(G[3])
    I = np.zeros(len(nus), complex)
    shift = None; epsmax = 0.0; nbad = 0
    t0 = time.time()
    for it, tt in enumerate(ts):
        for ip, pp in enumerate(ph):
            for ia, aa in enumerate(al):
                H = np.zeros((len(nus), G[0].size), complex)
                lw = lw4 + R[2].logW(tt) + R[5].logW(pp) + R[6].logW(aa)
                for sg in (1.0, -1.0):
                    o = P.point(G[0], G[1], np.full(G[0].size, sg*np.sqrt(tt)), G[2], G[3], np.full(G[0].size, pp),
                                np.full(G[0].size, aa), Delta, want_Q=(kind == 'E'))
                    lh = o['L'] - lw
                    if shift is None:
                        shift = float(np.max(lh[np.isfinite(lh)]))
                    h = np.exp(lh - shift)
                    if kind == 'E':
                        h = h*o['Q']
                    bad = ~np.isfinite(h)
                    nbad += int(bad.sum()); h = np.where(bad, 0.0, h)
                    epsmax = max(epsmax, float(np.max(np.abs(o['eps'][~bad])))) if (~bad).any() else epsmax
                    H += h[None, :]*np.exp(0.5j*np.outer(nus, np.where(bad, 0.0, o['eps'])))
                ang = W[2][it]*W[5][ip]*W[6][ia]
                Hm = H.reshape((len(nus),) + shp)
                T = np.einsum('vabcd,av,bv,cv,dv->v', Hm, W[0], W[1], W[3], W[4], optimize=True)
                I += ang*T
        if verbose:
            print(f'   D={Delta:g} {kind} t-node {it + 1}/{len(ts)}: {time.time() - t0:.0f}s, max|eps| {epsmax:.1e}, nbad {nbad}', flush=True)
    tag = '' if tuple(np.round(nus, 6)) == NUS_DEF else '_nu' + '-'.join('%g' % v for v in nus)
    fn = os.path.join(HERE, 'cache', 'dq_D%g_%s_%s%s.npz' % (Delta, kind, '-'.join(map(str, ns)), tag))
    np.savez(fn, nus=nus, I=I, shift=shift, epsmax=epsmax, nbad=nbad, ns=np.array(ns), Delta=Delta)
    return fn, nus, I*np.exp(shift)


def load(Delta, kind, ns):
    d = np.load(os.path.join(HERE, 'cache', 'dq_D%g_%s_%s.npz' % (Delta, kind, '-'.join(map(str, ns)))))
    return d['nus'], d['I']*np.exp(float(d['shift'])), d


if __name__ == '__main__':
    Delta = float(sys.argv[1]); kind = sys.argv[2]; ns = [int(v) for v in sys.argv[3].split(',')]
    nus = [float(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else NUS_DEF
    t0 = time.time()
    fn, nus, I = run(Delta, kind, ns, nus)
    print(f'pcyc_dquad D={Delta:g} {kind} n={ns}: {time.time() - t0:.0f}s -> {fn}')
    for v, x in zip(nus, I):
        print(f'   nu={v:5.2f}  I = {x.real:+.10e} {x.imag:+.10e}i')
