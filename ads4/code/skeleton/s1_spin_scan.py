"""
s1: exact two-pair <JJJJ> moments over nu at fixed (Delta, J): R_{kern,h}(k, nu) = exact/leading, with the tail estimate of the
dropped small-u region.  Usage: python3 s1_spin_scan.py Delta J [grid] [numax] [dnu]   (dnu < 0: nu = x m on the x-grid 0, |dnu|, 2|dnu|, ... up to
numax (in nu); default numax = min(0.65 m, 16), x-step 0.05; k = 2, 3, 4, 6)
Writes cache/s1_<Delta>_<J>_<grid>.pkl  (dict (k, nu) -> ratios) and prints a table.
"""
import sys, os, time, pickle
from common import CACHE, info
import numpy as np, mpmath as mp
import spin

if __name__ == "__main__":
    Delta, J = float(sys.argv[1]), int(sys.argv[2])
    grid = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != 'none' else None
    m, beta, eta = info(Delta, J); m = float(m)
    numax = float(sys.argv[4]) if len(sys.argv) > 4 else min(0.65*m, 16.0)
    dnu = float(sys.argv[5]) if len(sys.argv) > 5 else -0.05
    nus = np.arange(0.0, numax + 1e-9, dnu) if dnu > 0 else np.array([x*m for x in np.arange(0.0, 0.70001, -dnu) if x*m <= numax + 1e-9])
    t0 = time.time()
    mom = spin.load_case(Delta, J, grid=grid)
    print(f"# s1: Delta={Delta} J={J} m={m:.3f} beta={float(beta):.3f} grid={grid or 'std(fourpt cache)'}  load {time.time()-t0:.0f} s", flush=True)
    fn = os.path.join(CACHE, f"s1_D{Delta:.6f}_J{J}_{grid or 'std'}.pkl")
    res = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    for k in (2, 3, 4, 6):
        print(f"#  k={k}:  nu   x     R_even_same  R_odd_same   R_even_opp   R_odd_opp    (tail est. max)", flush=True)
        for nu in nus:
            key = (k, round(float(nu), 6))
            if key not in res:
                res[key] = spin.ratios(mom, k, float(nu))
                pickle.dump(res, open(fn, 'wb'))
            r = res[key]
            tl = max(v[1] for v in r.values())
            print(f"    {nu:5.2f} {nu/m:.3f}  " + "  ".join(f"{r[(kern, h)][0]:11.6f}" for h in ('same', 'opp') for kern in ('even', 'odd'))
                  + f"   {tl:.1e}", flush=True)
