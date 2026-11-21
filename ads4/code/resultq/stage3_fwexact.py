"""Result Q, Stage 3 (RESULT_Q_STAGE3.md §3): exact scalar-skeleton FORWARD moments on near-gap / small-bhat heavy states.

For a d = 3 scalar block (Delta, J) with identical external scalars Df (default 3 = Delta_T), compute with the exact CHMRS Mellin
representation (ads4/code/exactskel/exactmom.py, read-only import, [DV] there to 1e-11 for k = 4, 1e-10 for odd k):

    ex_k(nu) := m^{2k} C_{k,nu}[Delta, J]       k = 4, 5   (k = 5: Gamma-ratio continuation in k, exactmom docstring; a DEFINITION [M])

at small nu (nu = 0 itself is a removable 0/0 of the residue formula: the poles T = A +- i nu merge), and extrapolate to nu = 0 by a
cubic in nu^2 through nu = 0.1, 0.2, 0.3, 0.4; an independent fit through nu = 0.15, 0.25, 0.35, 0.45 gives the error estimate.
Outputs per state:  ex4(0), d ex4/d nu^2 |_0, ex5(0) (+ fit errors, max cancellation ratio).

States: Lam = sqrt((Dg-2)(Dg-1)); twist >= Dg; m^2 = alpha (alpha + 2J + 1), alpha = Delta - J - 2;  m/Lam on MFL (>= gap), J in JS.
Usage:  OMP_NUM_THREADS=1 python3 stage3_fwexact.py Dg [iworker nworkers]   -> cache/stage3_fw_<Dg>_w<i>.pkl (incremental, resumable)
"""
import os, sys, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import exactmom as E

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
DF = float(os.environ.get('S3_DF', '3'))
N0 = int(os.environ.get('S3_N0', '200'))
JS = [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24]
MFL = [1.03, 1.06, 1.1, 1.15, 1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]
NU_A = np.array([0.1, 0.2, 0.3, 0.4]); NU_B = np.array([0.15, 0.25, 0.35, 0.45])


def lam(Dg):
    return np.sqrt((Dg - 2.0) * (Dg - 1.0))


def states(Dg):
    L = lam(Dg); out = []
    for J in JS:
        mg = np.sqrt((Dg - 2.0) * (Dg - 1.0 + 2 * J))
        ms = [mg] + [f * L for f in MFL if f * L > mg * 1.0005]
        if Dg >= 300 and J >= 16:          # near the gap these have bhat >= 25: covered by the leading-Regge model (stage 12)
            if J in (18, 22):
                continue
            ms = [x for x in ms if x >= 1.5 * L * 0.9999]
        for m in ms:
            al = 0.5 * (-(2 * J + 1) + np.sqrt((2 * J + 1) ** 2 + 4 * m * m))
            out.append((float(al + J + 2), J, float(m / L)))
    return out


def fit0(nus, v):
    """cubic in nu^2 through 4 points: returns (c0, c1) = value and d/d(nu^2) at nu = 0."""
    V = np.vander(nus ** 2, 4, increasing=True)
    c = np.linalg.solve(V, v)
    return c[0], c[1]


def one(Delta, J):
    hc = E.HeavyC(Delta, J, DF, N0=N0)
    rec = {}
    canc = 0.
    for k in (4, 5):
        va, ca, _ = hc.value(k, NU_A, diag=True)
        vb, cb, _ = hc.value(k, NU_B, diag=True)
        a0, a1 = fit0(NU_A, va); b0, b1 = fit0(NU_B, vb)
        rec[k] = dict(v0=a0, d0=a1, ev0=abs(a0 - b0), ed0=abs(a1 - b1), va=va, vb=vb)
        canc = max(canc, ca.max(), cb.max())
    rec['canc'] = canc
    return rec


def main():
    Dg = int(sys.argv[1]); iw = int(sys.argv[2]) if len(sys.argv) > 2 else 0; nw = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    fn = os.path.join(CACHE, 'stage3_fw_%d_Df%g_w%d.pkl' % (Dg, DF, iw))
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    st = states(Dg)
    # balance by cost ~ J^2: deal states round-robin after sorting by J
    st = sorted(st, key=lambda s: (s[1], s[2]))[iw::nw]
    for (Delta, J, mfl) in st:
        key = (round(Delta, 6), J)
        if key in db:
            continue
        t = time.time()
        r = one(Delta, J); r['mfl'] = mfl; r['Dg'] = Dg
        db[key] = r
        pickle.dump(db, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
        print(f"Dg={Dg} J={J:2d} m/Lam={mfl:.4f} Delta={Delta:.3f}: ex4(0)={r[4]['v0']:.8f} (+-{r[4]['ev0']:.1e}) "
              f"d4={r[4]['d0']:.4e} (+-{r[4]['ed0']:.1e}) ex5(0)={r[5]['v0']:.8f} canc {r['canc']:.1e}  {time.time() - t:.1f}s", flush=True)


if __name__ == '__main__':
    main()
