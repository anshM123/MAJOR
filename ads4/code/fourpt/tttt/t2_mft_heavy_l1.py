"""
MFT-normalised two-pair heavy formula for <JJJJ>, absolute constant [D from the exact Fock closed forms of t1_mft_l1].
Per MFT operator of helicity h (contribution p_h (-1/8)^J W[hhat_h, hhat_h]), with the functional normalised as CHMRS / SPINNING_MFT
(a_{Delta_O,Delta_O}/gamma_a^2 and the heavy density q = 2 n_J (2Delta-3)/(pi m^2 p)):
   C_h := a_{2,2}/gamma_a^2 * [p_h 8^{-J} (-1)^J Pihat^q[W[hhat_h,hhat_h]]] * pi m^2/(2 n_J (2Delta-3)),  n_J = 16 pi (2J+1).
Prediction of fourpt/predict4 for (-1)^J Pihat: (pi/4) 4^4/b^(0) gamma_a^2 m^{-2a} sum_h V V M M K.  Result: C_h = (1/(4 m^{2k})) sum_{h in class}
(unit helicity tensors) M M K x Q_h, Q_h := C_h/[(1/(4 m^{2k})) sum M M K] -> 1, Weyl-even.  Output: the factor  C_h / [(1/(8 m^{2k})) sum M M K]  exactly (no moments needed).
"""
import os, sys
import mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'e2e'))
from formulas import log_b
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def PFK_mp(D1, D2, Delta, J):
    n = (Delta - D1 - D2 - J)/2; h = mp.mpf(3)/2
    rf = lambda a, k: mp.gamma(a + k)/mp.gamma(a)
    return (rf(D1-h+1, n)*rf(D2-h+1, n)*rf(D1, J+n)*rf(D2, J+n) /
            (mp.gamma(n+1)*mp.factorial(J)*rf(J+h, n)*rf(D1+D2+n-2*h+1, n)*rf(D1+D2+2*n+J-1, J)*rf(D1+D2+n+J-h, n)))

def factor(Delta, J, h):
    D = mp.mpf(Delta)
    Ap, Am = D - mp.mpf(3)/2, J + mp.mpf(1)/2; m2 = Ap**2 - Am**2; m = mp.sqrt(m2)
    Y = D**2 - 3*D + J**2 + J
    Pi = (J-1)*J*(J+1)*(J+2) if h == 'opp' else D*(D-1)*(D-2)*(D-3)
    p8 = PFK_mp(2, 2, D, J)*0 + PFK_mp(1, 1, D, J)/(8*Y**2*Pi)          # p_h 8^{-J}
    V = Am**2*m2*(Ap**2 + Am**2)/2 if h == 'opp' else Ap**2*m2*(Ap**2 + Am**2)/2
    a22 = 4*(mp.gamma(2)*mp.gamma(mp.mpf(3)/2))**2
    nJ = 16*mp.pi*(2*J + 1)
    # C_h * 8 m^{2k} / (sum M M K) ;  m^{-2a} m^{2k} = m^{-6} (a = k + 3)
    return a22*p8*(mp.pi/4)*4**4*mp.exp(-log_b(D, J, 0))*V**2*m**(-6)*mp.pi*m2/(nJ*(2*D - 3))*4

if __name__ == "__main__":
    pr("# Q_h = C_h / [(1/(4 m^{2k})) sum_{h in class} M M K]  (fourpt/tttt/t2_mft_heavy_l1.py); must -> 1 with 1/m^2 corrections")
    for beta in (0.4, 1.0, 2.0):
        for mm in (20, 100, 1000, 10000):
            J = int(round((2*mm*mp.sinh(beta/2) - 1)/2)); m_ = (2*J + 1)/(2*mp.sinh(mp.mpf(beta)/2))
            Delta = (m_*mp.exp(-beta/2) + m_*mp.exp(beta/2) + 3)/2
            m = mp.sqrt((Delta - J - 2)*(Delta + J - 1))
            qo = factor(Delta, J, 'opp'); qs = factor(Delta, J, 'same')
            pr(f"  beta={beta} m={mp.nstr(m, 6)} J={J}: Q_opp = {mp.nstr(qo, 12)}  m(Q_opp-1) = {mp.nstr(m*(qo-1), 4)}, m^2(Q_opp-1) = {mp.nstr(m*m*(qo-1), 6)};  Q_same = {mp.nstr(qs, 12)}  m^2(Q_same-1) = {mp.nstr(m*m*(qs-1), 6)}")
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t2_mft_heavy_l1.txt'), 'w').write("\n".join(out) + "\n")
