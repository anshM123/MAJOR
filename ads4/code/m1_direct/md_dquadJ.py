"""md_dquadJ.py -- direct D' per-nu image of the literal CLS-D k = 2 EE functional on a spin-J block (md_spinJ integrand), all
conserved parity-even structure pairs (a: <TTO_Delta>, c: <TTO~_{3-Delta}>) -> matrix I[a, c](nu).  Product rule as md_dquad with
t-weights t^6 e^{-lam t}, lam = 1 - (5 + J)/Delta (measured decay toward the B zero, md_test_spinJ); requires Delta - J > 5.
Usage: python3 md_dquadJ.py Delta J n1,...,n7 [nus] -> out/dqJ_D<Delta>_J<J>_<grid>.npz
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import pcyc_dquad as PQ
import md_spinJ as SJ
from md_dquad import NUS_DEF
np.seterr(all='ignore')

SPEC = [('lag', 6.0), ('jac', 2.25, 1.75, -0.5, 0.5), ('tjac', 1.75, -0.5), ('lag', 6.0),
        ('jac', 0.25, -0.25, 0.5, -0.5), ('cos', 0.5, 1.0), ('sin', -0.5, -1.0)]


def rules(ns, Delta, J):
    lam = 1.0 - (5.0 + J)/Delta
    assert lam > 0, 'D\' integral diverges for Delta - J <= 5'
    out = []
    for sp, n in zip(SPEC, ns):
        out.append(PQ.Rule('lag', n, sp[1], lam) if sp[0] == 'lag' else PQ.Rule(sp[0], n, *sp[1:]))
    return out


def run(Delta, J, ns, nus=NUS_DEF, verbose=True):
    nus = np.asarray(nus, float)
    SJ.conserved(Delta, J); SJ.conserved(3.0 - Delta, J)
    R = rules(ns, Delta, J)
    W = [r.weights(nus) for r in R]
    tx, rx, ts, ty, ry, ph, al = (r.x for r in R)
    G = np.meshgrid(tx, rx, ty, ry, indexing='ij')
    G = [g.ravel() for g in G]
    shp = (ns[0], ns[1], ns[3], ns[4])
    lw4 = R[0].logW(G[0]) + R[1].logW(G[1]) + R[3].logW(G[2]) + R[4].logW(G[3])
    I = None
    shift = None; nbad = 0
    t0 = time.time()
    for it, tt in enumerate(ts):
        for ip, pp in enumerate(ph):
            for ia, aa in enumerate(al):
                lw = lw4 + R[2].logW(tt) + R[5].logW(pp) + R[6].logW(aa)
                Hs = None
                for sg in (1.0, -1.0):
                    o = P.point(G[0], G[1], np.full(G[0].size, sg*np.sqrt(tt)), G[2], G[3], np.full(G[0].size, pp),
                                np.full(G[0].size, aa), Delta)
                    lh = o['L'] - lw
                    if shift is None:
                        shift = float(np.max(lh[np.isfinite(lh)]))
                    rat = SJ.ratio(o['g'], Delta, J)                                    # (na, nc, N)
                    h = np.exp(lh - shift)[None, None]*rat
                    bad = ~np.isfinite(h)
                    nbad += int(bad.sum()); h = np.where(bad, 0.0, h)
                    ph_ = np.exp(0.5j*np.outer(nus, np.where(np.isfinite(o['eps']), o['eps'], 0.0)))   # (nnu, N)
                    contrib = h[:, :, None, :]*ph_[None, None]
                    Hs = contrib if Hs is None else Hs + contrib
                ang = W[2][it]*W[5][ip]*W[6][ia]
                Hm = Hs.reshape(Hs.shape[:3] + shp)
                T = np.einsum('acvwxyz,wv,xv,yv,zv->acv', Hm, W[0], W[1], W[3], W[4], optimize=True)
                I = ang*T if I is None else I + ang*T
        if verbose:
            print(f'   D={Delta:g} J={J} t-node {it + 1}/{len(ts)}: {time.time() - t0:.0f}s, nbad {nbad}', flush=True)
    fn = os.path.join(HERE, 'out', 'dqJ_D%g_J%d_%s.npz' % (Delta, J, '-'.join(map(str, ns))))
    np.savez(fn, nus=nus, I=I, shift=shift, nbad=nbad, ns=np.array(ns), Delta=Delta, J=J)
    return fn, nus, I*np.exp(shift)


if __name__ == '__main__':
    Delta = float(sys.argv[1]); J = int(sys.argv[2]); ns = [int(v) for v in sys.argv[3].split(',')]
    nus = [float(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else NUS_DEF
    t0 = time.time()
    fn, nus, I = run(Delta, J, ns, nus)
    print(f'md_dquadJ D={Delta:g} J={J} n={ns}: {time.time() - t0:.0f}s -> {fn}')
    for k, v in enumerate(nus):
        print(f'   nu={v:5.2f}  ' + '  '.join(f'{x.real:+.8e}{x.imag:+.1e}i' for x in I[:, :, k].ravel()))
