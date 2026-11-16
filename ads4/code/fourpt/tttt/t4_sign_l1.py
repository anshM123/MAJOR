"""
Referee M2: sign of the opposite-sign term with physical labels, for covariantly defined kernels.
At each Regge pair (in, out) = (1,3), (2,4) with the actual polarisations of the operators and x_in - x_out:
   even:  D1 = z_in.z_out + (z_in.n)(z_out.n),    odd:  D5 = eps(z_in, z_out, n),   n = unit vector along x_in - x_out.
In the frame of fourpt/moments4: x_13 = -2y (n = -yhat), x_24 = 2x/x^2 = -2x/|x|^2 (n = -xhat), z_{2,4} = I(xhat) z'_{2,4}; det I = -1 and
I xhat = -xhat, so  D5_13 D5_24 = eps(z1,z3,yhat) eps(z2',z4',xhat)  and D1 is unchanged: on the circles z(th) = yhat + cos th e_par + sin th e_perp:
   D1 = cos(th_out - th_in),  D5 = sin(th_out - th_in)   (both pairs, det-rule frames).
Kernel action (real pairing int k G):  K[G] = sum_q w_L(q1,q3) w_R(q2,q4) Pihat^q,  w_even(+-1,-+1) = 1/2,  w_odd(-1,+1) = i/2, w_odd(+1,-1) = -i/2.
Flat-space analogue: D5 at a pair with preserved helicity measures the helicity of the incoming particle (same orientation rule at both
pairs), so D5 (x) D5 weighs h_1 h_2: +1 for CHL 'same' (J_z = h_1 - h_2 = 0, flat ++ of arXiv:2201.06602), -1 for 'opposite' (flat +-).
Computed: K_even and K_odd on W[hhat_same, hhat_same] and W[hhat_opp, hhat_opp], predicted (predict4) and exact (cached modes).
"""
import os, sys, glob
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from predict4 import prediction
from tmoments4 import TensorMoments4, DATADIR
from t5_analysis_l1 import hel_vectors
from cases import info
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
WE = {(1, -1): 0.5, (-1, 1): 0.5}
WO = {(-1, 1): 0.5j, (1, -1): -0.5j}

def action(w, k, nu, cL, cR, mom=None, Delta=None, J=None):
    tot_p = 0j; tot_e = 0j
    for qL, wl in w.items():
        for qR, wr in w.items():
            q = (qL[0], qL[1], qR[0], qR[1])
            pd = complex(prediction(1, Delta, J, k, nu, q, cL, cR))
            tot_p += wl*wr*pd
            if mom is not None:
                val = 0j
                for a in range(2):
                    for b in range(2):
                        if cL[a]*cR[b] == 0: continue
                        v, sh = mom.moment(k, nu, (a, b), q)
                        val += cL[a]*cR[b]*complex(v)*float(mp.exp(sh))
                tot_e += wl*wr*(-1)**J*val
    return tot_p, tot_e

if __name__ == "__main__":
    pr("# covariant even/odd pair kernels on the exact helicity waves (fourpt/tttt/t4_sign_l1.py)")
    files = sorted(glob.glob(os.path.join(DATADIR, 'modes4_l1_*.npz')))
    for f in files:
        mom = TensorMoments4.load(f)
        m, beta, _ = info(mom.Delta, mom.J)
        if float(m) < 60: continue
        o, s = hel_vectors(mom.Delta, mom.J)
        for k in (6,):
            for nu in (0.0, 1.0, 3.0):
                row = []
                for name, c in (('same', s), ('opp', o)):
                    ep, ee = action(WE, k, nu, c, c, mom, mom.Delta, mom.J)
                    op, oe = action(WO, k, nu, c, c, mom, mom.Delta, mom.J)
                    row.append(f"{name}: K_odd/K_even pred {op.real/ep.real:+.5f} exact {oe.real/ee.real:+.5f} (K_even exact/pred {ee.real/ep.real:.4f})")
                pr(f"  beta={float(beta):.2f} m={float(m):.1f} k={k} nu={nu}: " + ";  ".join(row))
    open(os.path.join(HERE, 'out_t4_sign_l1.txt'), 'w').write("\n".join(out) + "\n")
