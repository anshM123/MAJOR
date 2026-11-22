"""
rd_divisors (RD_CLOSE.md section 1.2): local exponents of the scalar integrand and of Q = Phi_EE/(a Phi_scal) on every component of
the singular divisor that could carry a pole of Q.

For each divisor f = 0 (a building block of rd_core.geometry): pick a random complex raw point, move one coordinate by Newton onto
f = 0 (generic point of the divisor), then step off along a random complex direction, z = p + t n, t = 10^-2 ... 10^-7, and fit
   e_scal = d sum_k Re(e_k) log|f_k| / d log|f|        (scalar; nu-independent part; = Re of the exponent)
   e_Q    = d log|Q| / d log|f|                       (Q evaluated in 60-digit arithmetic, rd_core.Qratio_mp)
Integer e_Q < 0 = pole of Q of that order; the divisor is then admissible for the EE only if it is a branch divisor of the scalar
(non-integer e_scal) or a zero of B^{Delta/2} strong enough to absorb the pole.
Usage: python3 rd_divisors.py -> cache/rd_divisors.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import rd_core as C
np.seterr(all='ignore')

DIVS = ['FP1', 'FP2', 'FP3', 'm15', 'm25', 'm35', 'm45', 'm20', 'm40', 'Nn', 'c10', 'c30', 'c13', 'mB', 'mA', 'negxsq', 'vmuN']
NOTE = {'FP1': 'm12 = m34 = 0 (m12 = -FP1/x^2)', 'FP2': 'm14 = m23 = 0 (scalar simple pole)', 'FP3': '1 - x^2 y^2 = 0',
        'm15': 'pole of B (y side)', 'm25': 'pole of B (x side)', 'm35': 'zero of B (y side)', 'm45': 'zero of B (x side)',
        'm20': 'x2 - x0 null', 'm40': 'x4 - x0 null', 'Nn': 'V_{0,42} = 0 (A divisor)', 'c10': 'A divisor', 'c30': 'A divisor',
        'c13': 'A divisor', 'mB': 'y^2 = 0 (u_y v_y - y_p^2)', 'mA': 'y0 = 0 (meas)', 'negxsq': 'x^2 = 0', 'vmuN': "v' - u' = 0 (zero of f)"}


def fval(z, y0ref, key):
    g = C.geom_plain(z, y0ref)
    return g[key], g['y0']


def onto(key, z, y0, coord, it=60):
    for _ in range(it):
        f, y0 = fval(z[:, None], np.array([y0]).ravel(), key)
        f = np.ravel(f); y0 = complex(np.ravel(y0)[0])
        h = 1e-7*max(1.0, abs(z[coord]))
        zp = z.copy(); zp[coord] += h
        fp, _ = fval(zp[:, None], np.array([y0]), key)
        d = (np.ravel(fp)[0] - f[0])/h
        step = f[0]/d
        z = z.copy(); z[coord] -= step
        if abs(step) < 1e-15*max(1.0, abs(z[coord])):
            break
    f, y0 = fval(z[:, None], np.array([y0]).ravel(), key)
    return z, complex(np.ravel(y0)[0]), abs(np.ravel(f)[0])


def run(Delta=300.0, ntries=3, seed=5):
    rng = np.random.default_rng(seed)
    F = C.Flow(1.0, Delta, np.ones(7))
    out = [f'# rd_divisors: local exponents at generic points of each divisor (Delta = {Delta:g}; scalar part nu-independent)',
           '#   divisor                         e_scal (fits)                 e_Q (fits, 60-digit Q)          verdict']
    for key in DIVS:
        es, eq = [], []
        for tr in range(ntries):
            for attempt in range(40):
                z = (0.3*np.exp(rng.normal(size=7)*0.7))*np.exp(1j*rng.uniform(-np.pi, np.pi, 7))
                z[0] = 0.5*(rng.normal() + 1j*rng.normal())
                y0 = complex(np.ravel(C.geom_plain(z[:, None])['y0'])[0])
                coord = rng.integers(1, 7)
                zz, yy, res = onto(key, z, y0, coord)
                if res < 1e-11 and np.all(np.isfinite(zz)):
                    break
            n = rng.normal(size=7) + 1j*rng.normal(size=7); n /= np.linalg.norm(n)
            ts = 10.0**-np.arange(2, 8)
            Ls, Lq, Lf = [], [], []
            yref = yy
            for t in ts:
                zt = zz + t*n
                g = C.geom_plain(zt[:, None], np.array([yref]))
                yref = complex(np.ravel(g['y0'])[0])
                vals = {k: v for k, v in g.items() if not isinstance(v, list)}
                Ls.append(float(np.ravel(F.logabs_parts(vals))[0])); Lf.append(float(np.log(np.abs(np.ravel(g[key])[0]))))
                Lq.append(np.log(abs(C.Qratio_mp(zt, yref, Delta, dps=60))))
            Ls, Lq, Lf = map(np.array, (Ls, Lq, Lf))
            es.append(np.polyfit(Lf[2:], Ls[2:], 1)[0]); eq.append(np.polyfit(Lf[2:], Lq[2:], 1)[0])
        es = np.array(es); eq = np.array(eq)
        eqr = np.round(np.median(eq))
        isint = abs(np.median(es) - np.round(np.median(es))) < 0.02
        if eqr >= 0:
            ver = 'Q regular (no pole)'
        elif not isint:
            ver = f'Q pole of order {-int(eqr)} on a scalar BRANCH divisor (in S)'
        else:
            ver = f'Q pole of order {-int(eqr)} on an INTEGER-exponent divisor of the scalar'
        out.append(f'  {key:7s} {NOTE[key]:30s} {np.array2string(es, precision=3):28s} {np.array2string(eq, precision=3):28s} {ver}')
    return out


if __name__ == '__main__':
    L = run()
    L.append('# note: m15/m25 carry -Delta/2 and m35/m45 +Delta/2 from B^{Delta/2}; these lamF exponents are included in e_scal.')
    print('\n'.join(L))
    open(os.path.join(HERE, 'cache', 'rd_divisors.txt'), 'w').write('\n'.join(L) + '\n')
