"""
z1 (EE_DG300_SMALLNU.md): location of the fold nu*(Delta, J) of the symmetric crystal saddle of the scalar flow function Phi_nu
(n3_saddles / n3_flow, imported read-only) at large J (24, 28, 32) for the Delta_g = 300 near-gap rows and the heavy references.
Symmetric saddle continued from nu = 20 (or higher) downwards in steps of 0.25 (Newton + sigma-symmetrisation); fold = where the
soft 4x4 determinant -> 0 / the symmetric Newton stops converging; then a symmetry-breaking kick tests for the A+- pair.
Usage: OMP_NUM_THREADS=2 python3 z1_fold.py Delta J [nu_top] [step]   -> appends to cache/z1_fold.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import n3_saddles as SD
import n4_thimble as N4
import l1_crystal as Cr


def scan(Delta, J, nu_top=24.0, step=0.25, nu_min=3.0, pr=print):
    w = np.zeros(7, complex); cprev = None
    last = None
    nu = nu_top
    rows = []
    while nu >= nu_min - 1e-9:
        cn, _ = Cr.saddle(nu, Delta)
        if cprev is not None:
            w = N4.wshift(w, cprev, cn)
        try:
            w2, g, H, P0, c = SD.newton(w.copy(), float(nu), Delta, J, sym=True, it=80)
        except Exception as ex:
            pr(f"   nu={nu:.3f}: newton exception {ex}"); break
        gm = float(np.abs(g).max())
        ds = complex(np.linalg.det(-H[1:5, 1:5]))
        ev = np.sort(np.abs(np.linalg.eigvals(-H)))
        rows.append((nu, gm, ds, ev[0]))
        pr(f"   Delta={Delta:g} J={J} nu={nu:7.3f}: |grad| {gm:.1e}  det_soft {ds.real:+.3e}{ds.imag:+.3e}i  min|eig| {ev[0]:.3f}  "
           f"Phi0 {P0.real:+.3f}")
        if gm > 1e-6 or not np.all(np.isfinite(w2)):
            break
        w = w2; cprev = cn; last = nu
        nu -= step
    return rows, last


def kick(Delta, J, nu, pr=print):
    """does an asymmetric pair exist at nu (below the fold)?  n4_thimble.track's first step only."""
    try:
        recs = N4.track(Delta, J, [nu], kinds=())
        a = {wh: N4.asym(recs[wh][0]['wc']) for wh in ('A+', 'A-')}
        pr(f"   kick at nu={nu}: A+- exist, asym {a}, Phi0 {recs['A+'][0]['Phi0']:.4f} / {recs['A-'][0]['Phi0']:.4f}")
        return True
    except Exception as ex:
        pr(f"   kick at nu={nu}: {ex}")
        return False


if __name__ == '__main__':
    Dl, J = float(sys.argv[1]), int(sys.argv[2])
    top = float(sys.argv[3]) if len(sys.argv) > 3 else 24.0
    st = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
    out = open(os.path.join(HERE, 'cache', 'z1_fold.txt'), 'a')
    def pr(s):
        print(s, flush=True); out.write(s + '\n'); out.flush()
    t = time.time()
    pr(f"## fold scan Delta={Dl:g} J={J} from nu={top} step {st}")
    rows, last = scan(Dl, J, top, st, pr=pr)
    pr(f"## Delta={Dl:g} J={J}: last symmetric saddle at nu={last}  ({time.time() - t:.0f}s)")
    if last is not None:
        for d in (0.05, 0.3, 1.0):
            if kick(Dl, J, last - d, pr=pr):
                break
