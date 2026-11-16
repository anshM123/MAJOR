"""
Block-level test of the two-pair MFT normalisation (coordinator goal 4): leading-twist (tau = 2, n = 0) lightcone limit of the MFT
correlator <J_1 J_3><J_2 J_4> on the real section, decomposed with the exact doubly-spinning partial waves of fourpt/spinblock2
(imported, not modified) and the exact Fock-space coefficients (t1_mft_l1 closed form, n = 0: hhat_opp = hhat_same):
     F_MFT  =  sum_{J >= 2} p(N = J - 2, J) (-1/8)^J W[hhat, hhat]^{(Delta = J + 2, J)}  x (1 + O(eps)),
p = 8^J P_FK(1,1; J+2, J)/(8 Y^2 (J-1)J(J+1)(J+2)).  x_12 = T(1, (1-eps) n) nearly null; z_1 = z_2 = n_12 (the null direction of x_12),
which suppresses higher twists (as in SPINNING_MFT.md sect. 5.1); z_3, z_4 various.  Control: (+1/8)^J.
"""
import os, sys
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'mft'))
from spinblock2 import SpinBlock4
from t2_validate_l1 import W_eval          # fourpt validation helper (imported)
from mft_exact import PFK
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
d = lambda a, b: -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def hhat(Delta, J):
    Y = Delta**2 - 3*Delta + J**2 + J; m2 = (Delta - 1.5)**2 - (J + 0.5)**2
    v1 = {'AC': (m2 - 2)*(m2 + 2), 'AE': -(Delta + J)*(Delta - J - 1)*Y, 'BC': -(Delta + J)*(Delta - J - 1)*Y, 'BE': Y**2, 'H': 0.0}
    v2 = {'AC': 2.0, 'AE': 0.0, 'BC': 0.0, 'BE': 0.0, 'H': Y}
    f = -(Delta - 1)*(Delta - 2)*(m2 - 2)       # same; at n = 0 equal to the opp one
    return {w: v1[w] + f*v2[w] for w in v1}

def p_mft(J):
    D = J + 2; Y = D*D - 3*D + J*J + J
    return float(sp.Integer(8)**J*PFK(1, 1, D, J)/(8*Y**2*(J - 1)*J*(J + 1)*(J + 2)))

def F_mft(X, Z):
    x13 = X[1] - X[3]; x24 = X[2] - X[4]
    f1 = (d(Z[1], Z[3])*d(x13, x13) - 2*d(Z[1], x13)*d(Z[3], x13))/abs(d(x13, x13))**3
    f2 = (d(Z[2], Z[4])*d(x24, x24) - 2*d(Z[2], x24)*d(Z[4], x24))/abs(d(x24, x24))**3
    return f1*f2

if __name__ == "__main__":
    pr("# block-level test of the <JJJJ> MFT normalisation, leading twist (fourpt/tttt/t3_lightcone_l1.py)")
    bases = {'A': (np.zeros(3), np.array([5.0, 1.0, 0.2]), np.array([-4.0, 0.3, -0.5]), np.array([1.0, 0.0]), 1.0),
             'B': (np.array([0.2, -0.1, 0.3]), np.array([4.0, 0.5, 1.0]), np.array([-3.5, -0.6, 0.2]), np.array([0.6, 0.8]), 0.8)}
    Jmax = 30
    pols = [((1, 0, 0), (0, 1, 0)), ((1, np.cos(2.6), np.sin(2.6)), (1, np.cos(0.7), -np.sin(0.7))), ((1, -1, 0), (1, 0, 1))]
    for bname, (x2, x3, x4, nn, T) in bases.items():
        n12 = np.array([1.0, *nn])
        for eps in (1e-3, 1e-4, 1e-5):
            x1 = x2 + T*np.array([1.0, *((1 - eps)*nn)])
            X = {1: x1[None], 2: x2[None], 3: x3[None], 4: x4[None]}
            rows = []
            for (z3, z4) in pols:
                Z = {1: n12[None, None], 2: n12[None, None], 3: np.array(z3, float)[None, None], 4: np.array(z4, float)[None, None]}
                if (z3, z4) == pols[0]:
                    Z[3] = n12[None, None]*0 + np.array([1.0, *nn])[None, None]; Z[4] = Z[3]
                F = F_mft({i: X[i][0] for i in X}, {i: Z[i][0, 0] for i in Z})
                S = Sc = 0.0
                for J in range(2, Jmax + 1):
                    sb = SpinBlock4(1, float(J + 2), J, N=400)
                    h = hhat(float(J + 2), J)
                    W = W_eval(sb, {'w': (h, h)}, X, Z)['w'][0, 0]
                    p = p_mft(J)
                    S += p*(-1/8)**J*W; Sc += p*(1/8)**J*W
                    last = abs(p*W/8**J/S)
                rows.append((F/S, F/Sc, last))
            pr(f"  point {bname} eps={eps:.0e}: ratio F_MFT/sum = " + "  ".join(f"{r:.9f}" for r, _, _ in rows)
               + "   | ratio-1: " + " ".join(f"{r-1:.1e}" for r, _, _ in rows) + "   | control (+1/8)^J: " + " ".join(f"{c:.3f}" for _, c, _ in rows)
               + f"   | last J term {rows[-1][2]:.0e}")
    open(os.path.join(HERE, 'out_t3_lightcone_l1.txt'), 'w').write("\n".join(out) + "\n")
