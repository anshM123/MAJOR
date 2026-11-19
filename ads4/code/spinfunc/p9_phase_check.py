"""
p9 (PARITY_ODD_ROUTEA.md sect. 12): origin of the e^{-i pi (Delta - Delta_ref)} state/reference phase, and the closed form of N.

  caches   : (i)  Im logK = pi (Delta + J) exactly, for every shadow entry (p2_shadow.pkl: odd sector; l9_shadow_k2.pkl: even sector);
             (ii) Im(logC_state - logC_ref) = 0 mod 2 pi at equal nu (all p3 / p8_k2 Laplace records, even and odd J).
  even J   : direct (n1_run.point, RQMC) state/reference ratio at fixed nu for Delta = 100, 100.5, 101, 104, 105 (ref 2e4), in a
             canonical real structure basis (M' = C^T M C, C = pinv(sol(Delta)) sol(100); the O~ basis cancels in M = E inv(b)^T),
             with and without the correction logK -> logK - i pi (Delta + J).   Read-only use of n1_*/l5_*; nothing cached.
  N        : N(Delta, J) = tan^2(pi (Delta - J)/2) at the p5 crystal rows and the PARITY_ODD_HEAVY sect. 4 Regge points.
Usage: python3 p9_phase_check.py caches | even J kind x | N        (outputs appended to cache/p9_phase_check.txt)
"""
import os, sys, glob, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

OUT = os.path.join(HERE, 'cache', 'p9_phase_check.txt')


def pr(s=""):
    print(s, flush=True)
    open(OUT, 'a').write(s + "\n")


def caches():
    for fn in ('p2_shadow.pkl', 'l9_shadow_k2.pkl'):
        db = pickle.load(open(os.path.join(HERE, 'cache', fn), 'rb'))
        w = max(abs(v[1].imag - np.pi*(k[0] + k[1])) for k, v in db.items())
        fr = sorted(set(round(k[0] % 2, 4) for k in db))
        pr(f"(i) {fn}: {len(db)} entries, max |Im logK - pi(Delta+J)| = {w:.1e}; Delta mod 2 values present: {fr[:8]}")
    worst, n = 0.0, 0
    for fn in glob.glob(os.path.join(HERE, 'cache', 'p3_pts_C*_J[4-9].pkl')) + glob.glob(os.path.join(HERE, 'cache', 'p8_k2_*.pkl')):
        db = pickle.load(open(fn, 'rb'))
        refs = {round(k[1], 6): v for k, v in db.items() if k[0] == 20000.0 and k[2] == 1.0}
        for k, v in db.items():
            if k[0] < 2e4 and k[2] == 1.0 and round(k[1], 6) in refs:
                d = (v['logC'] - refs[round(k[1], 6)]['logC']).imag
                worst = max(worst, abs((d + np.pi) % (2*np.pi) - np.pi)); n += 1
    pr(f"(ii) Laplace records: {n} state/reference pairs at equal nu, max |Im(logC_s - logC_r) mod 2pi| = {worst:.1e}")


def even(J=4, kind='C2', x=0.5):
    os.environ.setdefault('N1_M', '12'); os.environ.setdefault('N1_MREF', '12')
    import n1_run as R
    import l5_shadow_matrix as SM
    import l5_spin as S
    nu = x*R.mphys(100.0, J)
    Ds = [100.0, 100.5, 101.0, 104.0, 105.0, 20000.0]
    out = {}
    for D in Ds:
        t = time.time()
        sh = SM.shadow_matrix(D, J, 0.5, ntrip=6, seed=7)
        r = R.point(kind, J, D, nu)
        out[D] = dict(logK=sh['logK'], H=r['E'] @ np.linalg.inv(sh['b']).T, logN=r['logN'])
        pr(f"  even J={J} {kind} Delta={D} nu={nu:.4f}: Im logK/pi = {sh['logK'].imag/np.pi:.9f}, max|Im b|/|b| = "
           f"{np.abs(sh['b'].imag).max()/np.abs(sh['b']).max():.1e}, RQMC rel err {np.abs(r['err']).max()/np.abs(r['E']).max():.1e} ({time.time()-t:.0f}s)")
    r0 = S.conserved(100.0, J)[1]
    def can(D, M):
        C = np.linalg.pinv(S.conserved(D, J)[1]) @ r0
        return C.T @ M @ C
    ref = out[20000.0]; Mr = can(20000.0, ref['H'])
    pr(f"  arg/pi of the state/reference ratio (canonical real basis; elements 11 12 21 22; per-row real constant dropped):")
    for D in Ds[:-1]:
        s = out[D]
        dl = (s['logN'] - s['logK']) - (ref['logN'] - ref['logK'])
        v = np.exp(1j*dl.imag)*can(D, s['H'])/Mr
        f = lambda w: " ".join(f"{np.angle(z)/np.pi:+.4f}" for z in w.ravel())
        pr(f"    Delta={D:7.2f}: as computed {f(v)} | logK - i pi(Delta+J): {f(v*np.exp(1j*np.pi*(D - 20000.0)))}")


def Ncheck():
    import mpmath as mp
    rows = [(1600.37, 5, '2.3175632436 (C2, C3H)'), (3200.37, 5, '2.3175632436'), (400.37, 5, '2.31756'), (1600.37, 7, '2.31756'),
            (3200.37, 7, '2.31756'), (1600.37, 19, '2.31756'), (800.37, 19, '2.31755'), (400.37, 4, '0.43149'), (1600.37, 4, '0.43149'),
            (3200.37, 4, '0.43149'), (43.697092, 19, '3.76609 (Regge moments)'), (78.320346, 35, '3.30023 (Regge moments)')]
    for D, J, meas in rows:
        pr(f"  N: Delta={D} J={J}: tan^2(pi(Delta-J)/2) = {mp.nstr(mp.tan(mp.pi*(D - J)/2)**2, 11)}   measured {meas}")


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'caches'
    pr(f"# p9_phase_check {mode} {' '.join(sys.argv[2:])}")
    if mode == 'caches':
        caches()
    elif mode == 'even':
        even(int(sys.argv[2]), sys.argv[3], float(sys.argv[4]))
    else:
        Ncheck()
