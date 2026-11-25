"""
fullpos2/f2_polyext.py (FULLPOS2.md §5): route D. Re-runs the FULLPOS per-J every-mass certificate (fullpos/fp_poly.py, imported
read-only; nothing in fullpos/ is written) on new spins J0..J1 with the patched K = 1000 moments (moments_1000_K1000_an.json,
acb.integral with the holomorphy flag).  Output: fullpos2/cache/f2_polyext_<J0>-<J1>.txt.   [certified] given the moment balls.
"""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
from flint import arb, ctx
import fp_common as FC
import fp_poly as FP

def main():
    Dg = 1000; K = 1000
    J0, J1 = int(sys.argv[1]), int(sys.argv[2]); step = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    d = json.load(open(os.path.join(FC.CACHE, f'moments_{Dg}_K{K}_an.json')))
    ctx.prec = max(ctx.prec, d['prec'])
    Z = [arb(s) for s in d['Z']]; Y = [arb(s) for s in d['Y']]
    lam = (Dg - 2)*(Dg - 1)
    Q = FP.Q_exact(K, lam)
    A = FP.Action(Dg, Z, Y, Q)
    fn = os.path.join(HERE, 'cache', f'f2_polyext_{J0}-{J1}_s{step}.txt')
    with open(fn, 'a') as fo:
        for J in range(J0, J1 + 1, step):
            for b in (0, 4):
                if (b == 0 and J % 2 == 1) or (b == 4 and J < 4) or J - b > K:
                    continue
                t0 = time.time()
                P = A.P_poly(J, b)
                wmax = arb(lam)/((J + Dg - 2)*(J + Dg - 1))
                ok, lo, nc = FP.certify_pos(P, arb(0), wmax)
                p1 = P(wmax)
                line = (f'J={J:4d} b={b} n={J-b:4d}  {"PASS" if ok else "FAIL"} cells={nc:3d}  P(wmax)={p1.str(6)}  '
                        f'min Bernstein lower={lo:.3e}  {time.time()-t0:.1f}s')
                print(line, flush=True); fo.write(line + '\n'); fo.flush()

if __name__ == '__main__':
    main()
