"""
bf_validate (BELOW_FOLD.md §3): scalar validation of the real-section method (bf_core, theta = 0) against the exact Mellin skeleton
(n3_exact.X = g C(nu+i) + g~ C(nu-i), records via bf_exact): state / reference (Delta_ref = 2e4, same J) ratio on common QMC points,
residual m^2 eps = m^2 [log(R_D/R_exact) - const] (const fitted on nu < 1.4), +- replica error.  Also the EE k = 4 precision diagnostic
and the zero of N^_4 (1.068482) of the Regge-limit EE k = 4 integral.
Usage: python3 bf_validate.py  -> cache/bf_validate.txt  (about 3 min on 2 cores)
"""
import os, sys, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.special import eval_legendre
import bf_core as B, bf_exact as BE, n3_exact as NX

lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


def scalar_check(Dst, J, m2=15, R=8, numax=26, sampler='gamma'):
    rs, rr = BE.get(Dst, J), BE.get(2e4, J)
    nus = rs['nus'][rs['nus'] < numax]
    Rex = np.array([NX.X(nu, rs)[0]/NX.X(nu, rr)[0] for nu in nus])
    out = []; SH = {}
    for seed in range(R):
        S = B.Sample(m2, seed)
        V = {}
        for Dl in (Dst, 2e4):
            tot = 0
            for n0 in (1, -1):
                P = B.scalar_parts(S, Dl, 0, n0)
                if J:
                    Z, _, _ = S.raw(1/Dl, n0)
                    P['G'] = P['G'] + np.log(eval_legendre(J, B.t_ZW(Z).real) + 0j)
                v, sh = B.combine(P, nus, Dl, S.logp, shift=SH.get((Dl, n0))); SH[(Dl, n0)] = sh
                tot = tot + v
            V[Dl] = tot
        out.append(V[Dst]/V[2e4])
    out = np.array(out); Rm = out.mean(0); dR = out.std(0, ddof=1)/np.sqrt(R)
    L = np.log(Rm/Rex); c = np.mean(L[nus < 1.4]); m = rs['m']
    return nus, m*m*(L - c), m*m*dR/np.abs(Rm), m


if __name__ == '__main__':
    t0 = time.time()
    pr('# bf_validate: real-section (future-cone) crystal integral vs exact Mellin skeleton; m^2 eps (c2-type units, relative x m^2)')
    for (D, J) in ((100.0, 0), (100.0, 2), (100.0, 4), (300.0, 0), (336.0, 24)):
        nus, e, d, m = scalar_check(D, J)
        pr(f'## Delta = {D}, J = {J} (m = {m:.2f}), 2^15 x 8 points, |Delta - 2e4| on common points')
        for i, v in enumerate(nus):
            pr(f'   nu = {v:6.2f}:  Re {e[i].real:+11.2f}  Im {e[i].imag:+11.2f}  +- {d[i]:9.2f}   ({abs(e[i].real)/max(d[i], 1e-9):.1f} sigma)')
    pr(f'# time {time.time() - t0:.0f}s')
    open(os.path.join(HERE, 'cache', 'bf_validate.txt'), 'w').write('\n'.join(lines) + '\n')
