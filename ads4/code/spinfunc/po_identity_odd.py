"""
po_identity_odd (PARITY_ODD_CLOSE.md sect. 1.5): the odd-J analogue of the OO <-> helicity identity, on the existing p3 K2 records with
the structure set [anti-symmetric parity-even (2), symmetric parity-odd (1)] (cache/p3_pts_C2_J{5,7}_anti.pkl, P3_ANTI=1 of
PARITY_ODD_ROUTEA sect. 4; n0 = +1).  As po_identity: null vectors w, v from M0; the order-0 rank-1 helicity basis of the anti-even 2x2
block (l8 rule; channels labelled by their shape: the one whose W-ratio to the odd entry is x-flat is the odd structure's partner);
kappa(nu) = M_OO / (dual-helicity diagonal) for all orders summed.  Read-only.  Usage: python3 po_identity_odd.py -> cache/po_identity_odd.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l7_dc2_analysis as An
import l8_analysis as L8

C = os.path.join(HERE, 'cache')
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


for J in (5, 7):
    db = pickle.load(open(os.path.join(C, 'p3_pts_C2_J%d_anti.pkl' % J), 'rb'))
    Ds = sorted(set(k[0] for k in db))
    pr(f"J={J} (anti-even(2) + odd(1)):")
    for D in Ds:
        recs = {k: v for k, v in db.items() if k[0] == D and v['M0'].shape[0] == 3}
        if len(recs) < 4:
            continue
        ks = sorted(recs)
        sv, vs = An.rank1([recs[k]['M0'][:2, :2] for k in ks])
        T = np.array([L8.realvec(v)[0] for v in vs]).T; Ti = np.linalg.inv(T)
        M0 = recs[ks[0]]['M0']
        w = np.linalg.lstsq(M0[:2, :].T, M0[2, :], rcond=None)[0]; v = np.linalg.lstsq(M0[:, :2], M0[:, 2], rcond=None)[0]
        wt = w @ T; vt = v @ T
        j = int(np.argmax(np.abs(wt)))
        kap = []
        for k in ks:
            r = recs[k]; tot = r['M0'] + r['M2'] + r['M4']
            kap.append(tot[2, 2]/(Ti @ tot[:2, :2] @ Ti.T)[j, j])
        kap = np.array(kap)
        pr(f"   Delta={D:g}: rank-1 sv3/sv1 {sv[2]/sv[0]:.1e}; w.T = {np.round(wt, 6)} (other channel / partner {abs(wt[1-j]/wt[j]):.1e}), "
           f"v.T = {np.round(vt, 6)} ({abs(vt[1-j]/vt[j]):.1e}); kappa over {len(ks)} nu: mean {np.mean(kap).real:.8g}, rel spread "
           f"{np.abs(kap/np.mean(kap) - 1).max():.1e}")
open(os.path.join(C, 'po_identity_odd.txt'), 'w').write("\n".join(out) + "\n")
