"""
w4 (EE_NEARGAP.md): tables from w2's per-state records (cache/w2_recs_<tag>.pkl):
  (a) per J: worst state (corrected min lambda/N_s), uncorrected, signed worst diagonal dF_hh/N_h, off-diagonal ratio, and the
      "diagonal-only" eigenvalue (off-diagonal correction switched off) -- to isolate the graviton diagonal correction;
  (b) per (J, m/Lam) grid: corrected | uncorrected | signed min diag/N_h;
  (c) margin factor t_d of the DIAGONAL correction alone per state (scale on dF_hh at which lambda_min crosses 0, off-diagonal kept),
      and the combined n.
Usage: python3 w4_tables.py <tag> [Dg ...] -> stdout + cache/w4_tables_<tag>.txt
"""
import os, sys, io, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


def lam(base, Ns, dF, t=1.0, td=None):
    D = np.array(dF, float).copy()
    if td is not None:
        D[np.diag_indices_from(D)] *= td
    Fm = np.diag(base) + t*D; n = 1/np.sqrt(Ns)
    return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]


def tdiag(r):
    """smallest positive factor on the diagonal correction (off-diagonal and base fixed) at which lambda_min < 0 (inf if none <= 1e3)."""
    f = lambda t: lam(r['base'], r['Ns'], r['dF'], 1.0, t) > 0
    if not f(0.0):
        return 0.0
    lo, hi = 0.0, 1.0
    while f(hi) and hi < 1e3:
        lo, hi = hi, 2*hi
    if hi >= 1e3:
        return np.inf
    for _ in range(40):
        mid = 0.5*(lo + hi); lo, hi = (mid, hi) if f(mid) else (lo, mid)
    return lo


def main(tag, Dgs, mode='ee'):
    st = pickle.load(open(os.path.join(HERE, 'cache', 'w2_recs_%s.pkl' % tag), 'rb'))[mode]
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    for Dg in Dgs:
        recs = st['recs'][Dg]
        cov = [r for r in recs if r['cov']]
        pr(f"\n#### Delta_g = {Dg} [{tag}, {mode}]: {len(recs)} states, {len(cov)} with crystal corrections")
        grid = sorted(set(r['mfl'] for r in recs))
        nval = next((n0 for n0 in grid if all(r['lam'] > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
        w = min(recs, key=lambda r: r['lam'])
        pr(f"  n = {nval}; min lam/N_s = {w['lam']:+.4f} at (m/Lam, J, Df, hs) = ({w['mfl']}, {w['J']}, {w['DF']}, {w['hs']})")
        pr("  per J: worst corrected lam/N_s [m/Lam] | uncorrected min | worst signed diag dF_hh/N_h [m/Lam] | max off | lam (diag only) | "
           "min t_d (diag factor to failure) [m/Lam]")
        for J in sorted(set(r['J'] for r in cov)):
            rj = [r for r in cov if r['J'] == J]
            wj = min(rj, key=lambda r: r['lam'])
            dmin = min(rj, key=lambda r: min(r['diag']))
            don = min(lam(r['base'], r['Ns'], np.diag(np.diag(r['dF']))) for r in rj)
            tds = [(tdiag(r), r['mfl']) for r in rj]
            td = min(tds)
            pr(f"   J={J:2d}: {wj['lam']:+.4f} [{wj['mfl']:.2f}{'h' if wj['DF'] == 'hyb' else ''}] | {min(r['lam0'] for r in rj):+.4f} | "
               f"{min(dmin['diag']):+.4f} [{dmin['mfl']:.2f}] | {max(r['off'] for r in rj):.4f} | {don:+.4f} | {td[0]:.2f} [{td[1]:.2f}]")
        Js = sorted(set(r['J'] for r in cov))
        pr("  grid (m/Lam x J): corrected lam/N_s (signed min diag/N_h); '.' = no state")
        mf = sorted(set(r['mfl'] for r in cov))
        pr("   J \\ m/Lam " + " ".join(f"{m:>14.2f}" for m in mf))
        for J in Js:
            row = []
            for m in mf:
                rr = [r for r in cov if r['J'] == J and abs(r['mfl'] - m) < 1e-9]
                if not rr:
                    row.append(f"{'.':>14s}"); continue
                r = min(rr, key=lambda r: r['lam'])
                row.append(f"{r['lam']:+.3f}({min(r['diag']):+.3f})")
            pr(f"   J={J:2d}       " + " ".join(row))
    open(os.path.join(HERE, 'cache', 'w4_tables_%s.txt' % tag), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main(sys.argv[1], [int(a) for a in sys.argv[2:]] or [300, 1000], os.environ.get('W4_MODE', 'ee'))
