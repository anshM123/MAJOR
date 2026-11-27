"""
Independent check of the t4-coefficient of the CDP d=3 shear sum rule using free fields.

S := (1/(pi eps)) Int_{-inf}^{inf} dw [rho_T(w) - rho_0(w)]/w  (rho_T = Im G_R^{xy,xy}(w, p=0), sign: rho/w -> eta > 0)

CDP prediction (1612.00609 eq. (1.5)/(4.3) at d = 3, P = eps/2):   S = 3/8 + t4/32.
Free Majorana / Dirac fermion: t4 = -4  -> S = 1/4.     Free real scalar: t4 = +4 -> S = 1/2.
(t4 normalisation: HM / Buchel et al. 0911.4257 = CDP = this paper; free scalar +4, Majorana -4.)

At p = 0 the spectral function of the free theory splits into
  (a) the pair ("interband") continuum, rho_T = rho_0 * (1 -/+ 2 n(w/2)) for fermions/bosons, and
  (b) the zero-frequency transport ("Drude") peak, (1/pi) Int rho_peak/w = W = Int_k |<k|T_xy|k>|^2 (-dn/dE).
In any weakly interacting deformation (b) becomes the hydrodynamic shear peak, rho/w = eta/(1+w^2 tau^2), whose
weight eta/tau = W is the kinetic-theory relaxation-time identity, so S = S_pair + W/eps.
T_xy at p = 0: improvement terms are total derivatives d_x d_y(...) and vanish, so T_xy = d_x phi d_y phi (scalar)
and T_xy = psi^dag (sigma_x k_y + sigma_y k_x)/2 psi (2-component Dirac, = 2 Majorana; ratios are per-species
independent).
All integrals are done analytically below and re-done by numerical quadrature as a check. Units T = 1.
"""
import numpy as np
from scipy.integrate import quad, dblquad
from mpmath import zeta

z3 = float(zeta(3))
nF = lambda x: 1.0/(np.exp(x) + 1.0)
nB = lambda x: 1.0/np.expm1(x)

# ---------- Dirac fermion (2-component), H = sigma.k ----------
# interband |<+|M|->|^2 = (k^2/4) cos^2 2theta ; intraband <s|M|s> = s (k/2) sin 2theta
# rho_pair(w>0) = pi Int d^2k/(2pi)^2 (k^2/4)cos^2(2th) [f(-k)-f(k)] delta(w-2k)  = (w^3/256) tanh(w/4)
rho0_F = lambda w: w**3/256.0
eps_F = 2*quad(lambda k: k*k*nF(k), 0, np.inf)[0]/(2*np.pi)          # particles + holes
pair_F = (2/np.pi)*quad(lambda w: rho0_F(w)*(np.tanh(w/4) - 1)/w, 0, np.inf)[0]
# W = Int d^2k/(2pi)^2 sum_s (k^2/4) sin^2 2th * beta f(1-f)
W_F = 2*(1/(4*np.pi**2))*np.pi*0.25*quad(lambda k: k**3*nF(k)*(1 - nF(k)), 0, np.inf)[0]
# direct 2D numerical check of the angular factors
W_F2 = 2*dblquad(lambda th, k: k*(k*k/4)*np.sin(2*th)**2*nF(k)*(1 - nF(k)), 0, 60, 0, 2*np.pi)[0]/(4*np.pi**2)

# ---------- real scalar ----------
# T_xy(p=0) = sum_k (k_x k_y/2E)(a a_{-k} + a^dag a^dag_{-k} + a a^dag + a^dag a)
# pair: rho(w>0) = (pi/2) Int d^2k/(2pi)^2 (k^2/4) sin^2 2th (1 + 2 n) delta(w - 2k) = (w^3/512)(1 + 2 n(w/2))
rho0_B = lambda w: w**3/512.0
eps_B = quad(lambda k: k*k*nB(k), 0, np.inf)[0]/(2*np.pi)
pair_B = (2/np.pi)*quad(lambda w: rho0_B(w)*(2*nB(w/2))/w, 0, np.inf)[0]
W_B = (1/(4*np.pi**2))*np.pi*0.25*quad(lambda k: k**3*nB(k)*(1 + nB(k)), 1e-12, np.inf)[0]
W_B2 = dblquad(lambda th, k: k*(k*k/4)*np.sin(2*th)**2*nB(k)*(1 + nB(k)), 1e-9, 60, 0, 2*np.pi)[0]/(4*np.pi**2)

if __name__ == '__main__':
    print('C_T-independent ratios (both species): rho_0 per eps at T=1')
    print('Dirac:  eps = %.10f (analytic 3 zeta3/(2pi) = %.10f)' % (eps_F, 3*z3/(2*np.pi)))
    print('        pair/eps = %.10f (analytic -1/8)' % (pair_F/eps_F))
    print('        W/eps    = %.10f, 2D-quad %.10f (analytic 3/8)' % (W_F/eps_F, W_F2/eps_F))
    SF = (pair_F + W_F)/eps_F
    print('        S = %.10f   CDP with t4=-4: %.10f' % (SF, 3/8 - 4/32))
    print('Scalar: eps = %.10f (analytic zeta3/pi = %.10f)' % (eps_B, z3/np.pi))
    print('        pair/eps = %.10f (analytic +1/8)' % (pair_B/eps_B))
    print('        W/eps    = %.10f, 2D-quad %.10f (analytic 3/8)' % (W_B/eps_B, W_B2/eps_B))
    SB = (pair_B + W_B)/eps_B
    print('        S = %.10f   CDP with t4=+4: %.10f' % (SB, 3/8 + 4/32))
    # eta/s cross-check of the kinetic identity is not needed; report the affine fit through the two free points
    slope = (SB - SF)/8.0
    icpt = (SB + SF)/2
    print('Line through free points: S = %.6f + %.6f t4   (CDP: 0.375 + 0.03125 t4)' % (icpt, slope))
    print('Einstein point t4 = 0 predicted by this line: %.6f ; holographic value checked in einstein_ads4_sumrule.py' % icpt)
