"""
po_identity (PARITY_ODD_CLOSE.md sect. 1): the OO <-> ++ identity on the K2 rows of the 3x3 [even(2), odd(1)] action
(cache/oo_p3_{EE,C3H}_J<J>_D<D>.pkl, n0 = +1, Laplace orders M0, M2, M4 separate, after the shadow inverse).
For every record (state rows and their references, Delta_ref = 2e4):
  (a) single-S^0-point null vectors:  M_{O,c} = sum_a w_a M_{a,c},  M_{a,O} = sum_c M_{a,c} v_c  (w, v from M0; residual per order);
  (b) w, v vs the l8 helicity basis T0(Delta) (rank-1 decomposition of the order-0 EE 2x2 block over x >= 0.15, columns ++, +-):
      w . t_{+-} = 0 (w is the ++ dual vector), v likewise;
  (c) kappa(nu) = M_OO / (T0^-1 M_EE T0^-T)_{++}, all orders summed: nu-independent (=> Gamma_OO == Gamma_++ identically);
  (d) per-row constant  kappa(Delta)/kappa(Delta_ref)  vs HYPOTHESES_CLOSE2's a_O/a_++.
Read-only. Usage: python3 po_identity.py  -> cache/po_identity.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l8_analysis as L8

C = os.path.join(HERE, 'cache')
DREF = 20000.0
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


def nullvecs(M):
    A = M[:2, :].T; b = M[2, :]; w = np.linalg.lstsq(A, b, rcond=None)[0]
    A2 = M[:, :2]; b2 = M[:, 2]; v = np.linalg.lstsq(A2, b2, rcond=None)[0]
    return w, v


def rel(M, w, v):
    r1 = np.linalg.norm(M[2, :] - w @ M[:2, :])/np.abs(M).max()
    r2 = np.linalg.norm(M[:, 2] - M[:, :2] @ v)/np.abs(M).max()
    return max(r1, r2)


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def analyse(kind, J, Ds):
    db = pickle.load(open(os.path.join(C, 'oo_p3_%s_J%d_D%g.pkl' % (kind, J, Ds)), 'rb'))
    R = {(k[0], k[1]): v for k, v in db.items() if k[2] == 1.0}
    res = {}
    for D in (Ds, DREF):
        recs = {k: v for k, v in R.items() if k[0] == D}
        Ree = {k: dict(M0=v['M0'][:2, :2], M2=v['M2'][:2, :2], M4=v['M4'][:2, :2], logC=v['logC'], logK=v['logK'],
                       slog=v['slog'], sS=v['sS'], m=v['m']) for k, v in recs.items()}
        B = (L8.basis(J, D, Ree) or L8.basis(J, D, Ree, xs_min=0.0)) if kind == 'EE' else None
        ws, worst, kap = [], 0.0, []
        for k in sorted(recs):
            r = recs[k]
            w, v = nullvecs(r['M0'])
            ws.append((w, v))
            for o in ('M0', 'M2', 'M4'):
                worst = max(worst, rel(r[o], w, v))
        w0, v0 = ws[0]
        spread = max(max(np.abs(w - w0).max(), np.abs(v - v0).max()) for w, v in ws)/np.abs(w0).max()
        res[D] = dict(w=w0, v=v0, spread=spread, worst=worst, B=B, recs=recs)
    return res


def ident(kind, J, Ds, BEE=None):
    res = analyse(kind, J, Ds)
    pr(f"{kind} J={J} row Delta={Ds:g}:")
    kap = {}
    for D in (Ds, DREF):
        r = res[D]
        B = r['B'] if BEE is None else BEE[D]['B']
        T, Ti = B['T'], B['Ti']
        wt = r['w'] @ T; vt = Ti.T.T @ r['v'] if False else (r['v'] @ T)     # components along (++, +-) columns
        ks = []
        for k in sorted(r['recs']):
            rc = r['recs'][k]
            tot = rc['M0'] + rc['M2'] + rc['M4']
            pp = (Ti @ tot[:2, :2] @ Ti.T)[0, 0]
            ks.append(tot[2, 2]/pp)
        ks = np.array(ks)
        kap[D] = ks
        pr(f"   Delta={D:g}: w = {np.round(r['w'], 6)} v = {np.round(r['v'], 6)} (spread over nu {r['spread']:.1e}); null-vector "
           f"residual (all orders, all nu) {r['worst']:.1e}; w.T0 = {np.round(wt, 6)} [+- comp / ++ comp = {abs(wt[1]/wt[0]):.1e}], "
           f"v.T0 = {np.round(vt, 6)} [{abs(vt[1]/vt[0]):.1e}]; T0 sv3 {B['sv3']:.1e}")
        pr(f"      kappa = M_OO/(T0 M_EE T0^T)_++ (K2 total): mean {np.mean(ks):.8g}, rel spread over {len(ks)} nu "
           f"{np.abs(ks/np.mean(ks) - 1).max():.1e}; (w.t++)(v.t++) = {wt[0]*vt[0]:.8g}")
    ratio = np.mean(kap[Ds])/np.mean(kap[DREF])
    pr(f"   per-row constant kappa(Delta)/kappa(Delta_ref) = {ratio.real:.6g} (Im {ratio.imag:.1e})")
    return res, ratio


if __name__ == '__main__':
    rows = [(18, 336), (18, 400), (22, 336), (12, 400), (14, 400), (16, 400), (20, 400)]
    pr("# po_identity: single-S^0 null vectors, ++ duality and kappa constancy on the oo_p3 K2 rows (n0 = +1)")
    for J, D in rows:
        rEE, rat = ident('EE', J, D)
        ident('C3H', J, D, BEE=rEE)
    open(os.path.join(C, 'po_identity.txt'), 'w').write("\n".join(out) + "\n")
