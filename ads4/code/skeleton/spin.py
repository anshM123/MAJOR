"""
Two-pair <JJJJ> exact tensor Regge moments (fourpt/tmoments4, unchanged) reduced to the covariant even / odd pair kernels on the exact
helicity waves W[hhat_h, hhat_h], h = same, opp (fourpt/tttt/t4_sign_l1 conventions), and their leading (factorised) prediction
(fourpt/predict4).  MFT normalisation: Q_h of fourpt/tttt/t2_mft_heavy_l1 (exact, Weyl-even).

   K_even = sum_{qL,qR in {(1,-1),(-1,1)}} (1/2)(1/2) Pihat^{(qL,qR)},   K_odd = the same with weights w(-1,1) = i/2, w(1,-1) = -i/2.
   R_{kern,h}(k, nu) = (-1)^J K_kern[exact] / K_kern[predicted]          (-> 1 + O(1/m^2) at fixed nu; SPINNING_4PT.md 4.5, 5)
   MFT-normalised exact moment / leading:   Q_h R_{kern,h}.
"""
import os, sys
from common import CACHE, CODE, info
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(CODE, 'fourpt', 'tttt'))
from predict4 import prediction
from tmoments4 import TensorMoments4, get_moments, default_N
from t5_analysis_l1 import hel_vectors
from t2_mft_heavy_l1 import factor as Qfactor

WE = {(1, -1): 0.5, (-1, 1): 0.5}
WO = {(-1, 1): 0.5j, (1, -1): -0.5j}
KERNELS = {'even': WE, 'odd': WO}

def load_case(Delta, J, grid=None, N=None):
    """TensorMoments4 on the standard grid (cache of fourpt if present, else ours) or on a named grid of common.grid_kw."""
    if grid is None:
        fn0 = os.path.join(CODE, 'fourpt', 'cache', 'modes4', f"modes4_l1_D{Delta:.6f}_J{J}.npz")
        if os.path.exists(fn0):
            return TensorMoments4.load(fn0)
    from common import grid_kw
    gname = grid or 'std'
    fn = os.path.join(CACHE, 'modes4', f"modes4_l1_D{Delta:.6f}_J{J}_{gname}_N{N or default_N(1, Delta)}.npz")
    if os.path.exists(fn):
        return TensorMoments4.load(fn)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    kw = GRIDS[gname]
    mom = TensorMoments4(1, Delta, J, N or default_N(1, Delta), **kw)
    mom.save(fn)
    return mom

def _grids():
    from common import grid_kw
    return {'std': grid_kw(), 'u2': grid_kw(fu=2), 'u15': grid_kw(fu=1.5)}
GRIDS = _grids()

def kernel_action(mom, w, k, nu, c, predicted=True, xfac=None):
    """(-1)^J K[exact] and K[pred] for the helicity wave W[c, c] (c on (v_1, v_2))."""
    tot_p = 0j; tot_e = 0j
    J = mom.J
    for qL, wl in w.items():
        for qR, wr in w.items():
            q = (qL[0], qL[1], qR[0], qR[1])
            if predicted:
                tot_p += wl*wr*complex(prediction(1, mom.Delta, J, k, nu, q, c, c))
            val = 0j
            for a in range(2):
                for b in range(2):
                    v, sh = mom.moment(k, nu, (a, b), q, xfac=xfac)
                    val += c[a]*c[b]*complex(v)*float(mp.exp(sh))
            tot_e += wl*wr*(-1)**J*val
    return tot_e.real, tot_p.real

def ratios(mom, k, nu, tail=True):
    """dict (kern, h) -> (R, tail_estimate, K_pred)"""
    o, s = hel_vectors(mom.Delta, mom.J)
    out = {}
    for h, c in (('same', s), ('opp', o)):
        for kern, w in KERNELS.items():
            ee, ep = kernel_action(mom, w, k, nu, c)
            te = np.nan
            if tail:
                xf = 1.5
                e2, _ = kernel_action(mom, w, k, nu, c, predicted=False, xfac=xf)
                te = abs(ee - e2)/abs(ee)/(xf**(k + 0.5) - 1)
            out[(kern, h)] = (ee/ep, te, ep)
    return out

def Q(Delta, J, h):
    return float(Qfactor(mp.mpf(Delta), J, h))
