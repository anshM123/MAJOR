"""
gapscale/gs_state.py (GAP_SCALING.md §2): per-state Omega^EE action in GAP UNITS (J, u = Lam^2/m^2) as a function of Delta_g.

Decomposition used by every per-state budget (Supplement §S7; spinfunc/ee1000_common.hyb, read-only import):
    base(J, u; Lam) = B2 + B3            (flat NF02 partial waves: depend on (J, u, h) ONLY -- exactly Lam-independent)
                    + cv2 + cv3          (leading-Regge curvature: int psi_k [P_{-1/2+i Lam p}(cosh beta) - J0(b p)] x factor)
    lambda0 = base/N_s,  N_s = |B2 + cv2| + |B3 + cv3|
and the sign-free Lemma C bound (Supplement §S7, [P]):  |P_{-1/2+i nu}(cosh beta) - J0(nu B)| <= (1 - sech(beta/2)) + min(2, nu(B - beta)),
B = 2 sinh(beta/2) = (2J+1)/m.  At fixed (J, u):  b = (2J+1) sqrt(u) is fixed, beta = 2 asinh(b/(2 Lam)) decreases with Lam, and
   1 - sech(beta/2)                    decreases with Lam,
   nu (B - beta) = Lam p (B - beta) = b p (1 - asinh(s)/s),  s = b/(2 Lam),   decreases with Lam  (asinh(s)/s is decreasing in s),
so the Lemma-C curvature bound LC(J, u; Lam) is pointwise MONOTONE DECREASING in Lam [P].
This script tabulates, at fixed (J, u, h): flat value/N_flat, exact curvature (and Lam^2 x curvature, to exhibit the 1/Lam^2 law),
the LC bound, the sign-free lower bound  lb = (B2 + B3 - LC)/(N_flat + LC)  and lambda0, at Dg = 300, 1000, 2000, 5000, 20000.
Usage: python3 gs_state.py -> cache/state.txt, cache/state.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E
Q = E.Q

OUT = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def lc_bound(st, J, mfl, Lm, h):
    """sign-free Lemma C bounds on |cv2|, |cv3| at the state (same quadrature nodes as hyb)."""
    m = mfl*Lm
    B = (2*J + 1)/m
    beta = 2*np.arcsinh(B/2)
    nu = st['pn']*Lm
    lc = (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, nu*(B - beta))
    b2 = (np.abs(st['p2'])*2*lc).sum()
    b3 = (np.abs(st['p3'])*st['u']*lc).sum()
    return b2, b3


def one(Dg, J, mfl, h, npn):
    Lm = E.lam_of(Dg)
    st = E.hyb(Dg, J, mfl, h, npn)
    b2, b3 = lc_bound(st, J, mfl, Lm, h)
    flat = st['B2'] + st['B3']; Nf = abs(st['B2']) + abs(st['B3'])
    cv = st['cv2'] + st['cv3']
    LC = b2 + b3
    return dict(Lam=float(Lm), flat=float(flat), Nf=float(Nf), cv=float(cv), LC=float(LC),
                flatN=float(flat/Nf), lam0=float(st['base']/st['Ns']), lb=float((flat - LC)/(Nf + LC)),
                cvN=float(cv/Nf), LCN=float(LC/Nf))


def main():
    t0 = time.time()
    DGS = (300, 1000, 2000, 5000, 20000)
    # (J, m/Lam): fixed GAP-UNIT states.  m/Lam = 1.0 is the Dg -> oo limit of every twist edge.
    states = [(5, 1.0049925), (5, 1.0), (4, 1.0), (6, 1.0), (7, 1.0), (8, 1.0), (12, 1.0), (16, 1.0), (21, 1.0), (26, 1.026), (32, 1.0315),
              (48, 1.05), (12, 1.3), (32, 1.3), (64, 1.3), (5, 2.0), (32, 2.0), (100, 2.0)]
    res = {}
    pr("# gs_state: per-state Omega^EE hybrid pieces in gap units (J, m/Lam) vs Dg.  [flat: exact; curvature: exact LR; LC: Lemma C bound]")
    pr("# columns per Dg: lambda0 | cv/N_f | Lam^2 cv/N_f | LC/N_f | sign-free lb = (flat-LC)/(N_f+LC)")
    for (J, mf) in states:
        hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
        for h in hs:
            npn = 256 if J < 100 else 768
            row = []
            for Dg in DGS:
                if mf < E.mgap(Dg, J) - 1e-12 and mf != 1.0:
                    row.append(None); continue
                try:
                    r = one(Dg, J, mf, h, npn)
                except Exception as ex:          # noqa
                    row.append(None); pr(f"   ! {Dg} {J} {mf} {h}: {ex}"); continue
                row.append(r)
            res[f"{J}_{mf}_{h}"] = row
            fl = next(r for r in row if r is not None)
            pr(f"\nJ={J:3d} m/Lam={mf:<9} {h}: flat/N_f = {fl['flatN']:+.6f} (Lam-independent; m >= Lam twist-allowed at Dg >= "
               f"{'oo' if mf == 1.0 else '%.0f' % (2*J/(mf**2 - 1) + 1)})")
            for Dg, r in zip(DGS, row):
                if r is None:
                    pr(f"   Dg={Dg:6d}: below the twist edge"); continue
                pr(f"   Dg={Dg:6d}: lam0 {r['lam0']:+.6f} | cv/N {r['cvN']:+.3e} | Lam^2 cv/N {r['cvN']*r['Lam']**2:+.4e} | "
                   f"LC/N {r['LCN']:.3e} | lb {r['lb']:+.5f}")
    json.dump(res, open(os.path.join(HERE, 'cache', 'state.json'), 'w'))
    pr(f"\n# done in {time.time()-t0:.0f}s")
    open(os.path.join(HERE, 'cache', 'state.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
