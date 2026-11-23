r"""
pcyc_region (PCYC_CLOSE.md section 1): the physical integration domain of the CLS-D functional (arXiv:2311.04271, eqs.
scalarfunctional_plus_final, block_shadowrep0, gaugefixing) is the causal region 4 > 5 > 3 > 0 > 1 > 2.  In the gauge
x3 = -x1 = y, x4 = -x2 = x/(-x^2), x5 = e = (1,1,0) this script checks, on random REAL configurations, with an independent
Minkowski implementation (t, x^1, x^2 coordinates; u = t + x^1, v = t - x^1):

  (a) [phys]  <=>  [D']   where  phys = all causal relations of 4>5>3>0>1>2 (future-timelike / future-null as required),
                          D' = {x, y future timelike; 1 - s_x + q_x > 0, u_x < 1; 1 - s_y + q_y > 0, u_y < 1} x (x0 circle);
  (b) on D' every base f_k of a non-integer power of the scalar integrand (rd_core.scalar_factors) has ONE sign, A > 0, B > 0, and the
      EE/scalar ratio Q (rd_core.Qratio) is real: the physical (absolute-valued) integrand is therefore a z-independent constant
      times the analytic integrand on all of D' (one connected chain);
  (c) on D \ D' (the uncut real section): which relation fails, and the sign of B there.
Usage: python3 pcyc_region.py [N] -> cache/pcyc_region.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import rd_core as C
np.seterr(all='ignore')


def mink(p):
    """lightcone (u, v, perp) -> (t, x1, x2)."""
    return np.array([(p[0] + p[1])/2, (p[0] - p[1])/2, p[2]])


def sq(P):
    return -P[0]**2 + P[1]**2 + P[2]**2


def ftl(P, tol=0.0):            # future timelike
    return (sq(P) < -tol) & (P[0] > 0)


def fnull(P, tol=1e-9):          # future null (or within tol)
    return (np.abs(sq(P)) <= tol*(P[0]**2 + P[1]**2 + P[2]**2 + 1e-300)) & (P[0] > 0)


def sample(N, rng, wide=True):
    """random real x, y future timelike (wide ranges, so that D' and D \\ D' are both populated), x0 on the circle."""
    lo, hi = (1e-3, 30.0) if wide else (1e-3, 1.0)
    ux, vx, uy, vy = (np.exp(rng.uniform(np.log(lo), np.log(hi), N)) for _ in range(4))
    s = rng.uniform(-1, 1, N); sf = rng.uniform(-1, 1, N); psi = rng.uniform(-np.pi, np.pi, N)
    xp = s*np.sqrt(ux*vx); yp = sf*np.sqrt(uy*vy)
    cf = np.sqrt(1 - sf**2)
    xt0 = cf*np.cos(psi); y0 = np.sqrt(uy*vy)*np.sin(psi)
    z = np.array([xt0, ux, vx, uy, vy, xp, yp], dtype=complex)
    return z, y0.astype(complex), dict(s=s, sf=sf, psi=psi)


def predicates(g):
    x, y, x0 = (np.real(np.array(g[k])) for k in ('x', 'y', 'x0'))
    x4 = np.real(np.array(g['x4'])); x2 = -x4; x1 = -y; x3 = y
    e = np.array([np.ones_like(x[0]), np.ones_like(x[0]), np.zeros_like(x[0])])
    M = {k: mink(v) for k, v in dict(x=x, y=y, x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, e=e).items()}
    rel = {
        'x ftl': ftl(M['x']), 'y ftl': ftl(M['y']),
        '4>5': ftl(M['x4'] - M['e']), '5>3': ftl(M['e'] - M['x3']),
        '3>0 (null)': fnull(M['x3'] - M['x0']), '0>1 (null)': fnull(M['x0'] - M['x1']),
        '1>2': ftl(M['x1'] - M['x2']), '4>3': ftl(M['x4'] - M['x3']),
    }
    phys = np.all(np.array(list(rel.values())), axis=0)
    return phys, rel


def dprime(z):
    ux, vx, uy, vy, xp, yp = (np.real(z[i]) for i in range(1, 7))
    sx, qx = ux + vx, ux*vx - xp**2
    sy, qy = uy + vy, uy*vy - yp**2
    return (1 - sx + qx > 0) & (ux < 1) & (1 - sy + qy > 0) & (uy < 1)


def main(N=200000):
    rng = np.random.default_rng(20261123)
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    pr('# pcyc_region: physical domain of CLS-D (4>5>3>0>1>2, gauge x3=-x1=y, x4=-x2=x/(-x^2), x5=e) vs the cut real cycle D\'')
    tot = dict(both=0, phys_only=0, dp_only=0, neither=0)
    fails = {}
    signs = {}
    Bneg = 0; nDDc = 0; Qimag = 0.0; Amin = np.inf; Bmin = np.inf; nD = 0
    fac = [k for k, e in C.scalar_factors(4.0, 300.0)]
    expo = dict(C.scalar_factors(4.0, 300.0))
    for chunk in range(N//20000):
        z, y0, aux = sample(20000, rng, wide=(chunk % 2 == 0))
        g = C.geometry([z[i] for i in range(7)], y0val=y0)
        phys, rel = predicates(g)
        dp = dprime(z)
        tot['both'] += int((phys & dp).sum()); tot['phys_only'] += int((phys & ~dp).sum())
        tot['dp_only'] += int((~phys & dp).sum()); tot['neither'] += int((~phys & ~dp).sum())
        for k, v in rel.items():
            fails[k] = fails.get(k, 0) + int((~v & ~dp).sum())
        # signs of every factor base on D'
        for k in fac + ['negxsq']:
            v = np.asarray(g[k])[dp]
            if np.max(np.abs(v.imag)) > 1e-9*np.max(np.abs(v) + 1e-300):
                signs.setdefault(k, set()).add('complex')
            sg = set(np.unique(np.sign(v.real)).astype(int).tolist())
            signs.setdefault(k, set()).update(sg)
        A = (g['Nn']**2*g['c13']/(g['m24']*g['m40']*g['m20']*g['c10']*g['c30']))[dp]
        B = (g['m35']*g['m45']*g['m12']/(g['m15']*g['m25']*g['m34']))
        Amin = min(Amin, float(np.min(A.real))); Bmin = min(Bmin, float(np.min(B[dp].real)))
        Bneg += int((B[~dp].real <= 0).sum()); nDDc += int((~dp).sum())
        nD += int(dp.sum())
        if chunk < 2:
            idx = np.where(dp)[0][:300]
            gs = C.geometry([z[i][idx] for i in range(7)], y0val=y0[idx])
            Q = C.Qratio(gs, 300.0)
            Qimag = max(Qimag, float(np.max(np.abs(Q.imag)/np.abs(Q))))
    pr(f'  samples {N} (half wide u, v in [1e-3, 30], half in [1e-3, 1]); x, y future timelike by construction, x0 on the circle')
    pr(f'  confusion:  phys & D\' {tot["both"]}   phys only {tot["phys_only"]}   D\' only {tot["dp_only"]}   neither {tot["neither"]}')
    pr('  (a) phys <=> D\': ' + ('EXACT on all samples' if tot['phys_only'] == 0 and tot['dp_only'] == 0 else 'MISMATCH'))
    pr('  (c) on D \\ D\', relation failing (counts; several can fail):')
    for k, v in fails.items():
        pr(f'        {k:12s} {v}')
    pr(f'      B <= 0 on {Bneg} of {nDDc} points of D \\ D\' (B changes sign across the light cones of e: the uncut D leaves the physical chain)')
    pr('  (b) signs on D\' of every base f_k (exponent e_k of the scalar integrand at nu = 4, Delta = 300; non-integer or nu-dependent marked *):')
    for k in fac + ['negxsq']:
        e = expo.get(k, 0)
        mark = '*' if (abs(np.imag(e)) > 0 or abs(np.real(e) - round(np.real(e))) > 1e-12) else ' '
        pr(f'        {k:7s} e = {complex(e).real:+8.3f}{complex(e).imag:+.2f}i {mark}  signs {sorted(signs.get(k, []), key=str)}')
    pr(f'      min A on D\' = {Amin:.3e} (> 0);  min B on D\' = {Bmin:.3e} (> 0);  max |Im Q|/|Q| on 600 D\' points = {Qimag:.1e}')
    pr(f'  points of D\' sampled: {nD}')
    open(os.path.join(HERE, 'cache', 'pcyc_region.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200000)
