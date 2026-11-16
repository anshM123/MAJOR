"""
Exact scalar-skeleton Regge moments (identical external scalars of dimension Df; e2e.moments.ScalarMoments, unchanged) and their
normalisations:
   R_skel(k, nu; Df) = (-1)^J Pihat_{k,nu}[exact] / lead,  lead = (pi/4) 4^{2Df}/b^(0) gamma_{2Df+k-1}(nu)^2 P_{-1/2+i nu}(eta_AdS)/m^{2(2Df+k-1)}
   Nrat(Df) = MFT normalisation factor (CHMRS q_{Delta,J} a_Df / ...; uniform/skel), nu- and k-independent, = 1 + O(1/m^2):
   MFT-normalised skeleton moment / (its leading Regge value) = Nrat * R_skel   (this is m^{2k} Pihat-type C / 2P of uniform/).
Also the CHMRS order-8 series of the same quantity (seriesPi), for comparison.
"""
import os, pickle
from common import CACHE, info, grid_kw
import numpy as np, mpmath as mp
from moments import ScalarMoments
from formulas import gamma_a, log_b, Pnu
import skel

def lead(Delta, J, Df, k, nu):
    m, beta, eta = info(Delta, J)
    a = 2*Df + k - 1
    return mp.re(mp.exp(mp.log(mp.pi/4) + 2*Df*mp.log(4) - log_b(mp.mpf(Delta), J, 0) - a*mp.log(m*m))*gamma_a(a, nu)**2*Pnu(nu, eta))

def Nrat(Delta, J, Df, k=4):
    m, beta, eta = info(Delta, J); m = float(m); eta = float(eta)
    lq, la = skel.log_q(Delta, J, Df), mp.log(skel.a_Df(Df))
    nu = 0.0
    pref = mp.exp(lq - la - k*mp.log(m*m))*skel.gamma_a(2*Df + k - 1, nu)**2
    Pd = skel.Pders(nu, eta)
    t = skel.series_terms('seriesPi', k, Df, nu, eta, Pd=Pd)
    return float(pref*t[0]/lead(Delta, J, Df, k, nu))

def series_ratio(Delta, J, Df, k, nu, j=4):
    m, beta, eta = info(Delta, J); m = float(m); eta = float(eta)
    Pd = skel.Pders(nu, eta)
    t = skel.series_terms('seriesPi', k, Df, nu, eta, Pd=Pd)
    return float(sum(t[i]/(m*m)**i for i in range(j + 1))/t[0])     # = Nrat*R_skel predicted by the series / Nrat

class Skeleton:
    def __init__(self, Delta, J, grid='std'):
        self.Delta, self.J = Delta, J
        N = int(min(12000, max(600, 45*Delta)))
        kw = {'std': grid_kw(), 'u2': grid_kw(fu=2)}[grid]
        self.moms = {}
        self.N, self.kw = N, kw
        self.fn = os.path.join(CACHE, f"skel_D{Delta:.6f}_J{J}_{grid}.pkl")
        self.res = pickle.load(open(self.fn, 'rb')) if os.path.exists(self.fn) else {}

    def R(self, Df, k, nu):
        key = (float(Df), k, round(float(nu), 6))
        if key not in self.res:
            if Df not in self.moms:
                self.moms[Df] = ScalarMoments(self.Delta, self.J, Df, Df, self.N, **self.kw)
            s, lg = self.moms[Df].moment(k, nu)
            self.res[key] = float(s*mp.exp(lg)/lead(self.Delta, self.J, Df, k, nu))
        return self.res[key]

    def save(self):
        pickle.dump(self.res, open(self.fn, 'wb'))
