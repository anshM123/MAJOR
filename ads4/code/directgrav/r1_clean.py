"""
r1: re-extraction of the directgrav waves (A: m/Lam=1.04, J=4; C: 1.3, J=6) with the t19 cleaning of SPINNING_TTTT_EXT.md sect. 7
(fourpt/tttt/t19_k2.py imported READ-ONLY): mask cut 40 (production 25 was calibrated for l = 1), per-node noise floor from the imaginary
parts of the real-by-symmetry modes, power-law row tails at tiny t, light-cone wedge extrapolation; variants -> model spread, constant
tails + cut 32 -> envelope.  Kernels K*, A (odd), S_ff (info), leading mode; waves ss (++), oo (+-); k = 2..7; nu = the functional's Gauss
nodes with nu <= NUMAX (and 0, 1, 3).
Usage: python3 r1_clean.py cache/a1_gw_<tag>.pkl cache/gw_<tag>.npz [NUMAX]  -> cache/r1_<tag>.pkl (incremental)
"""
import os, sys, time, pickle, functools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
for p in (TT, os.path.join(TT, '..'), os.path.join(TT, '..', '..', 'e2e')):
    sys.path.insert(0, p)
import numpy as np
import t19_k2 as T19
from t13_analysis_l2 import Modes
T19.design = functools.lru_cache(None)(T19.design)          # same function, cached across kernels/waves
a1fn, wfn = sys.argv[1], sys.argv[2]
NUMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 10.7
a1 = pickle.load(open(a1fn, 'rb'))
NUS = [nu for nu in a1['NUS'] if nu <= NUMAX]
mo = Modes(wfn); mo.fn = wfn
tag = os.path.basename(wfn).replace('.npz', '')
out = os.path.join(HERE, 'cache', 'r1_%s.pkl' % tag)
db = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
c = T19.crit(mo)
print(f"# r1 {tag}: Delta={mo.Delta:.6f} J={mo.J}; nodes kept at cut 25/32/40/50: " +
      "/".join(str(int(((c > cc) & mo.g['ok']).sum())) for cc in (25, 32, 40, 50)) + f"; {len(NUS)} nu", flush=True)
t0 = time.time()
WAVES = [('same', 'same')] + ([('opp', 'opp')] if mo.J >= 4 else [])
for nu in NUS:
    for k in (2, 3, 4, 5, 6, 7):
        for wave in WAVES:
            for kn in ('K*', 'A', 'Sff', 'lead'):
                key = (k, nu, wave, kn)
                if key in db:
                    continue
                if kn == 'Sff':
                    r = T19.cleaned(mo, k, nu, wave, T19.KS[kn]); db[key] = dict(R=r['R'], raw=r['raw'], spread=np.nan, env=(np.nan, np.nan))
                else:
                    res, env, Rn, sp = T19.analyse([mo], [1.0], k, nu, wave, kn)
                    db[key] = dict(R=Rn[0], raw=res['nominal'][0]['raw'], spread=sp[0], env=env[0])
    pickle.dump(db, open(out + '.tmp', 'wb')); os.replace(out + '.tmp', out)
    print(f"  nu={nu:.3f} done ({time.time() - t0:.0f} s)", flush=True)
print("#  k   nu    | per wave: K* [spread] {env} | A [spread] | Sff | lead   (raw K*)")
for k in (2, 3, 4, 5, 6, 7):
    for nu in NUS:
        cells = []
        for wave in WAVES:
            g = lambda kn: db[(k, nu, wave, kn)]
            cells.append(f"{wave[0][0]}{wave[1][0]}: K* {g('K*')['R']:+.4f} [{g('K*')['spread']:.1e}] {{{g('K*')['env'][0]:.3f}..{g('K*')['env'][1]:.3f}}}"
                         f" A {g('A')['R']:+.4f} [{g('A')['spread']:.1e}] Sff {g('Sff')['R']:+.4f} lead {g('lead')['R']:+.4f} (raw K* {g('K*')['raw']:+.3e})")
        print(f"  {k}  {nu:6.3f} | " + " | ".join(cells), flush=True)
print(f"# total {time.time() - t0:.0f} s")
