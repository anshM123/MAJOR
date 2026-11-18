import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from scipy.integrate import quad
from nf_common import NF
from common import load
from rcommon import peval
F = NF(); FM, PM = load('M02')
for Dg in (100, 1000):
    Lam = np.sqrt((Dg-2.)*(Dg-1.))
    LN = quad(lambda t: float(F.psi2(t))/(t*t+2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
    LM = quad(lambda t: float(peval(FM.psi2, t))/(t*t+2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
    print(Dg, "L_NF/L_M02 =", LN/LM, "psi3_M02(0) =", float(peval(FM.psi3, 0.0)), "M02 k0-normalised?")
