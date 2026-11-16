"""
Absolute two-pair MFT constant for <TTTT> (l = 2), the analogue of t2_mft_heavy_l1 (same normalisation of the functional,
a_{Delta_O,Delta_O}/gamma_a^2 with Delta_O = 3, density q = 2 n_J (2 Delta - 3)/(pi m^2 p^MFT)).
Inputs:  p^MFT[X_h] = |C_2|^2 A_h P_FK(1,1)/2304  (t7b; |C_2|^2 = 2^11 2^J P(Delta)^2),
         leading heavy integrand of the diagonal wave  T^{hh} = lambda_h E/m^8  (t12/t14; E = 8 from the vertex data),
         lambda_h = 2^{2J+1} P(Delta~)/(A_h P(Delta)).
Then  p[X_h] 8^-J T^{hh} = 2^12 P(Delta) P(Delta~) P_FK(1,1) E/(2304 m^8):  A_h cancels (helicity independent), and
   Q_h := C_h / [(1/(4 m^{2k})) sum_{h in class} M M K]  = E/8 * (1 + O(1/m^2)).
Printed: Q_h/E*8 (the purely kinematic factor -> 1) at fixed beta, and Q_h with the measured finite-m E of t14 along the b0.4 ray.
Output: out_t15_mft_constant_l2.txt
"""
import os, sys
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from formulas import log_b
from t2_mft_heavy_l1 import PFK_mp
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def P(D, J): return mp.rf((2 - D + J)/2, 4)*mp.rf((D + J - 7)/2, 4)

def Q(D, J, E=8):
    D = mp.mpf(D); Ap, Am = D - mp.mpf(3)/2, J + mp.mpf(1)/2; m2 = Ap**2 - Am**2; m = mp.sqrt(m2)
    X = 2**12*P(D, J)*P(3 - D, -1 - J)*PFK_mp(1, 1, D, J)*E/(2304*m**8)            # p[X_h] 8^-J T^{hh}
    a33 = 4*(mp.gamma(3)*mp.gamma(mp.mpf(5)/2))**2
    nJ = 16*mp.pi*(2*J + 1)
    # identical algebra to t2_mft_heavy_l1.factor with l = 2: 4^{2(l+1)} = 4^6, m^{-2(2l+1)} = m^-10
    return a33*X*(mp.pi/4)*4**6*mp.exp(-log_b(D, J, 0))*m**-10*mp.pi*m2/(nJ*(2*D - 3))*4

if __name__ == "__main__":
    pr("# l = 2 absolute two-pair MFT constant (fourpt/tttt/t15_mft_constant_l2.py): Q_h = (E/8) x kinematic factor")
    for beta in (0.4, 1.0, 2.0):
        for mm in (20, 100, 1000, 10000):
            J = int(round((2*mm*mp.sinh(beta/2) - 1)/2)); m_ = (2*J + 1)/(2*mp.sinh(mp.mpf(beta)/2))
            D = (m_*mp.exp(-beta/2) + m_*mp.exp(beta/2) + 3)/2
            m = mp.sqrt((D - J - 2)*(D + J - 1))
            q = Q(D, J)
            pr(f"  beta={beta} m={mp.nstr(m, 6)} J={J}: Q_h(E = 8) = {mp.nstr(q, 12)},  m^2 (Q - 1) = {mp.nstr(m*m*(q - 1), 6)}")
    pr("# with the measured leading-vertex normalisation E(m) of t14 (ray b0.4): Q_h = E(m)/8 x kinematic factor")
    import pickle
    fn = os.path.join(HERE, 'cache', 't12_res_b0.4.pkl')
    if os.path.exists(fn):
        from t14_vertex_fit_l2 import lam, REP
        data = pickle.load(open(fn, 'rb'))
        for s, (Dv, Jv, m, R) in sorted(data.items()):
            D = mp.mpf(Dv.p)/Dv.q; J = mp.mpf(Jv.p)/Jv.q
            for h in ('same', 'opp'):
                E = mp.re(R[(h, h)][0][REP[(h, h)]]*m**8/lam(h, D, J))
                pr(f"  s={s} m={mp.nstr(m, 6)} {h}: E = {mp.nstr(E, 8)},  E/8 = {mp.nstr(E/8, 6)},  m^2 (E/8 - 1) = {mp.nstr(m*m*(E/8 - 1), 5)}")
    open(os.path.join(HERE, 'out_t15_mft_constant_l2.txt'), 'w').write("\n".join(out) + "\n")
