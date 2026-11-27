"""md_dquad.py -- DIRECT per-nu image of the literal CLS-D k = 2 EE functional (C^'_EE, j' = 0, f = v'-u') on a J = 0 block of
dimension Delta, as the absolutely convergent integral over the physical chain D' (CLS-D's causal domain 4>5>3>0>1>2, Lorentzian
shadow representation of the block; PCYC_CLOSE section 1), EXACT in Delta (no 1/Delta expansion; the cut at the light cones of x5
included).  This is a route independent of the detector-ladder formula (fp-target)/(fp-H): no light-ray matrix elements, no CDH
polynomials, no rest-frame reduction.

Integrand:  Phi_EE2 = a * Q2 * Phi_scal,   Q2 = u'v' prefh/(m12^2 m34^2) = Q_k4 / R,   R = ((v'-u')^2 - 1)/(4u'v'),
  Phi_scal = the scalar k = 2, J = 0 CLS-D integrand (rd_core.scalar_factors; pcyc_dlib), Q_k4 = rd_core.Qratio (= repro_engine EE k=4
  integrand / Phi_scal, RD_CLOSE section 1.1 [P + DV]).  Check of the k = 2 factor: w0_ee_common._gs (the project's EE k = 2 builder)
  differs from the scalar build by (v'-u') m12^-5 m34^-7/2 m35^-3/2 m45^-3/2 V042^-4 pref / [fsub m12^-3 m34^-3/2 m35^-3/2 m45^-3/2]
  = u'v' prefh/(m12^2 m34^2).  (md_check_k2.py verifies this point-wise against w0's factor list.)
Quadrature: the 7-D product rule of pcyc_dquad (imported read-only; its run() is NOT called, nothing is written to pcyc/cache) with
weights matched to the measured boundary exponents of Phi_EE2 (md_probe.py): t^6 e^{-(1-5/Delta)t} (t_x, t_y), r_x: (2.25, 1.75),
r_y: (0.25, -0.25), s: (1-t)^1.75, phi: cos^0.5, alpha: sin^-0.5.  Converges for Delta > 5.
Usage: python3 md_dquad.py Delta kind n1,...,n7 [nus]   kind in {S, E2}  -> out/dq_D<Delta>_<kind>_<grid>.npz
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import pcyc_dquad as PQ
import rd_core as C
np.seterr(all='ignore')

NUS_DEF = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0)
SPEC = {'S': [('lag', 6.0, None), ('jac', 2.25, 1.75, -0.5, 0.5), ('tjac', 1.75, -0.5), ('lag', 6.0, None),
              ('jac', 2.25, 1.75, 0.5, -0.5), ('cos', 4.5, 1.0), ('sin', -0.5, -1.0)],
        'E2': [('lag', 6.0, None), ('jac', 2.25, 1.75, -0.5, 0.5), ('tjac', 1.75, -0.5), ('lag', 6.0, None),
               ('jac', 0.25, -0.25, 0.5, -0.5), ('cos', 0.5, 1.0), ('sin', -0.5, -1.0)]}
LAMC = {'S': 1.0, 'E2': 5.0}


def Q2(g, Delta):
    Q, parts = C.Qratio(g, Delta, parts=True)
    up = g['m12']*g['m34']/(g['m13']*g['m24']); vp = g['m14']*g['m23']/(g['m13']*g['m24'])
    R = (parts['vmu']**2 - 1.0)/(4*up*vp)
    return np.real(Q/R)


def rules(kind, ns, Delta):
    out = []
    for sp, n in zip(SPEC[kind], ns):
        if sp[0] == 'lag':
            out.append(PQ.Rule('lag', n, sp[1], 1.0 - LAMC[kind]/Delta))
        else:
            out.append(PQ.Rule(sp[0], n, *sp[1:]))
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
    I = np.zeros(len(nus), complex); Iabs = np.zeros(len(nus))
    shift = None; epsmax = 0.0; nbad = 0
    t0 = time.time()
    for it, tt in enumerate(ts):
        for ip, pp in enumerate(ph):
            for ia, aa in enumerate(al):
                H = np.zeros((len(nus), G[0].size), complex)
                lw = lw4 + R[2].logW(tt) + R[5].logW(pp) + R[6].logW(aa)
                for sg in (1.0, -1.0):
                    o = P.point(G[0], G[1], np.full(G[0].size, sg*np.sqrt(tt)), G[2], G[3], np.full(G[0].size, pp),
                                np.full(G[0].size, aa), Delta, want_Q=False)
                    lh = o['L'] - lw
                    if shift is None:
                        shift = float(np.max(lh[np.isfinite(lh)]))
                    h = np.exp(lh - shift)
                    if kind == 'E2':
                        h = h*Q2(o['g'], Delta)
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
    fn = os.path.join(HERE, 'out', 'dq_D%g_%s_%s.npz' % (Delta, kind, '-'.join(map(str, ns))))
    np.savez(fn, nus=nus, I=I, shift=shift, epsmax=epsmax, nbad=nbad, ns=np.array(ns), Delta=Delta)
    return fn, nus, I*np.exp(shift)


def load(Delta, kind, ns):
    d = np.load(os.path.join(HERE, 'out', 'dq_D%g_%s_%s.npz' % (Delta, kind, '-'.join(map(str, ns)))))
    return d['nus'], d['I']*np.exp(float(d['shift'])), d


if __name__ == '__main__':
    Delta = float(sys.argv[1]); kind = sys.argv[2]; ns = [int(v) for v in sys.argv[3].split(',')]
    nus = [float(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else NUS_DEF
    t0 = time.time()
    fn, nus, I = run(Delta, kind, ns, nus)
    print(f'md_dquad D={Delta:g} {kind} n={ns}: {time.time() - t0:.0f}s -> {fn}')
    for v, x in zip(nus, I):
        print(f'   nu={v:5.2f}  I = {x.real:+.10e} {x.imag:+.10e}i')
