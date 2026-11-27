r"""
m1_remainder/m1_n0.py  [D-sketch for the identification; DV for the number]
The T-block normalisation n_0 of Lemma 1 (no_eft/NOTE.md): p_T omega[T] = n_0 G [K tau^2 Lam^8/57600 - 16 L_A].
Its tau-independent part is the weight integrated against the graviton pole of row L1, -C'_{EE,nu}|_E = 16 pi G/(nu^2+9/4),
in the units of Eq. (fp-target) (whose B is Einstein-calibrated, sec_fullpos.tex):  n_0 = int_0^oo g(nu) 16 pi/(nu^2+9/4) dnu / (16 L_A).
Usage: python3 m1_n0.py  -> prints n_0 at Delta_g = 300, 1000.
"""
import sys, os, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m1_weight import Weight
L16 = {300: 10734.61, 1000: 14339.59}
for Dg in (300, 1000):
    W = Weight(Dg)
    xg, wg = np.polynomial.legendre.leggauss(16)
    nu = np.array([p + 0.5 + 0.5*x for p in range(int(W.Lam) + 60) for x in xg])
    w = np.array([0.5*x for p in range(int(W.Lam) + 60) for x in wg])
    g = -9216*(4*nu**2 + 1)/((4*nu**2 + 9)*(4*nu**2 + 49))*W.chi_np(nu)*W.theta_np(nu)/W.Lam
    IT = np.sum(w*g*16*np.pi/(nu**2 + 9/4))
    print(f'Delta_g={Dg}: int g 16pi/(nu^2+9/4) = {IT:.6e};  n_0 = {IT/L16[Dg]:.5e} = {IT/L16[Dg]*W.Lam**2:.5f}/Lam^2')
