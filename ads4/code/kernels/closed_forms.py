"""
Closed forms (see KERNELS_D3.md, section 5).  Delta = dimension, J = spin of O;  all for the source |T,+2>.

Universal m-dependence (verified symbolically in Delta for all conserved structures, J <= 12):
   |ahat_{-k}|^2 / |ahat_0|^2 = W0(k) = (J+k)!/((J-k)! k!^2) * (4)_k^2 / ((Delta-1-k)_k (Delta-1)_k),   0 <= k <= J
   |ahat_{4+k}|^2 / |ahat_4|^2 = W4(k) = (J-4)!(J+4+k)!/(k!^2 (J-4-k)! (J+4)!) * (4)_k^2 / ((Delta-5-k)_k (Delta+3)_k)
   ahat_1 = ahat_2 = ahat_3 = 0.
Kernels for the canonical structures (normalisations: KERNELS_D3.md section 5.3 and data/structures_J*.json):
   |ahat_0|^2 = |lam_pp Kpp + i lam_odd Kodd|^2 / C_T,   |ahat_4|^2 = lam_pm^2 Kpm^2 / C_T.
"""
import mpmath as mp

def W0(k, Delta, J):
    if k < 0 or k > J:
        return mp.mpf(0)
    return (mp.factorial(J + k)/(mp.factorial(J - k)*mp.factorial(k)**2)*mp.rf(4, k)**2
            / (mp.rf(Delta - 1 - k, k)*mp.rf(Delta - 1, k)))

def W4(k, Delta, J):
    if k < 0 or k > J - 4:
        return mp.mpf(0)
    return (mp.factorial(J - 4)*mp.factorial(J + 4 + k)/(mp.factorial(k)**2*mp.factorial(J - 4 - k)*mp.factorial(J + 4))
            * mp.rf(4, k)**2/(mp.rf(Delta - 5 - k, k)*mp.rf(Delta + 3, k)))

def norm_O(Delta, J, m):
    return (mp.pi**2.5*mp.mpf(2)**(4 - 2*Delta - J)*mp.binomial(2*J, J + m)*mp.rf(Delta - J - 1, J - abs(m))
            * mp.rf(Delta - 1, abs(m))*mp.rgamma(Delta + J)*mp.rgamma(Delta - mp.mpf(1)/2))

NT_PLUS = mp.pi**2/3      # norm of T(e_+)|0> per unit C_T

# ---- amplitudes M_m(theta=0) per unit coefficient of the canonical structures (source T(e_+))
def M4_pm(Delta, J):      # even J >= 4, structure S^{+-}
    s = (Delta + J)/mp.mpf(2)
    return (-3j*(-1)**(J//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(2 - J)*mp.rf(Delta - J - 1, J - 4)
            * mp.rgamma(s)*mp.rgamma(Delta - 1)*mp.rgamma(1 - (Delta - J)/mp.mpf(2)))
def M0_pp(Delta, J):      # even J >= 4, structure S^{++}
    s = (Delta + J)/mp.mpf(2)
    return (3j*(-1)**(J//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(3 - J)*mp.rf(Delta - J - 1, J)
            * mp.rgamma(s)*mp.rgamma(Delta + 3)*mp.rgamma(1 - (Delta - J)/mp.mpf(2)))
def M0_odd(Delta, J):     # even J >= 4, parity-odd structure S^{odd}
    s = (Delta + J)/mp.mpf(2)
    return (3*(-1)**(J//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(4 - J)*mp.rf(Delta - J - 1, J)
            * mp.rgamma(s + mp.mpf(1)/2)*mp.rgamma(Delta + 3)*mp.rgamma((J + 3 - Delta)/mp.mpf(2)))
def M4_oddJ(Delta, J):    # odd J >= 5, the single structure
    s = (Delta + J)/mp.mpf(2)
    return (-3j*(-1)**((J - 5)//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(4 - J)*mp.rf(Delta - J - 1, J - 4)
            * mp.rgamma(s + mp.mpf(1)/2)*mp.rgamma(Delta - 1)*mp.rgamma((J + 3 - Delta)/mp.mpf(2)))

# ---- squared kernels  |K|^2 = 4 pi^2 |M|^2 / (N_O N_T^+)
def K2(M, Delta, J, m):
    return 4*mp.pi**2*abs(M)**2/(norm_O(Delta, J, m)*mp.pi**2/3)      # N_T^+ = pi^2/3 per unit C_T
def Kpm2(Delta, J):
    return K2(M4_pm(Delta, J), Delta, J, 4)
def Kpp2(Delta, J):
    return K2(M0_pp(Delta, J), Delta, J, 0)
def Kodd2(Delta, J):
    return K2(M0_odd(Delta, J), Delta, J, 0)
def KoddJ2(Delta, J):
    return K2(M4_oddJ(Delta, J), Delta, J, 4)

# ---- J = 0 (scalar): structures normalised as  coefficient(V1^2 V3^2) = 1  (even; KKSZ (5.110)),
#      beta_1 = Delta+1, beta_2 = Delta-3 (odd; CMT (D.3) with (beta_2-beta_1)/4 = -1)
def K0_even_scalar2(Delta):
    return (mp.mpf(27)/4*mp.pi**4.5*mp.gamma(Delta)*mp.gamma(Delta - mp.mpf(1)/2)
            * mp.rgamma(2 + Delta/mp.mpf(2))**4*mp.rgamma(3 - Delta/mp.mpf(2))**2*mp.rgamma((3 + Delta)/mp.mpf(2))**2)
def K0_odd_scalar2(Delta):
    return (mp.mpf(27)/4*mp.pi**4.5*mp.gamma(Delta)*mp.gamma(Delta - mp.mpf(1)/2)
            * mp.rgamma(2 + Delta/mp.mpf(2))**2*mp.rgamma((7 - Delta)/mp.mpf(2))**2*mp.rgamma((1 + Delta)/mp.mpf(2))**2
            * mp.rgamma((3 + Delta)/mp.mpf(2))**2)

# ---- top components (finite at the unitarity bound Delta = J + 1):  M_{J}[S^{+-}], M_{-J}[S^{++}], M_{-J}[S^{odd}]
def Mtop_pm(Delta, J):
    s = (Delta + J)/mp.mpf(2)
    return (-3j*(-1)**(J//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(2 - J)*mp.rf(4, J - 4)
            * mp.rgamma(s)*mp.rgamma(Delta - 1)*mp.rgamma(1 - (Delta - J)/mp.mpf(2)))
def Mtop_pp(Delta, J):
    s = (Delta + J)/mp.mpf(2)
    return (3j*(-1)**(J//2)*mp.exp(1j*mp.pi*s)*mp.pi**3*mp.mpf(2)**(3 - J)*mp.rf(4, J)
            * mp.rgamma(s)*mp.rgamma(Delta + 3)*mp.rgamma(1 - (Delta - J)/mp.mpf(2)))

# ---- J = 2 (non-conserved spin 2), structures normalised as in data/structures_J2.json
def M0_J2_even(Delta):
    s = Delta/mp.mpf(2) + 1
    return (3j*mp.pi**3*mp.exp(1j*mp.pi*s)*(Delta - 3)*(Delta - 2)*mp.rgamma(Delta + 2)*mp.rgamma(2 + Delta/mp.mpf(2))
            * mp.rgamma(3 - Delta/mp.mpf(2)))
def M0_J2_odd(Delta):
    s = Delta/mp.mpf(2) + 1
    return (mp.mpf(3)/4*mp.pi**3.5*mp.exp(1j*mp.pi*s)*mp.mpf(2)**(1 - Delta)*(Delta - 3)*(Delta - 2)
            * mp.rgamma(2 + Delta/mp.mpf(2))*mp.rgamma((Delta + 3)/mp.mpf(2))**2*mp.rgamma((7 - Delta)/mp.mpf(2)))

# ---- literature forms of the scalar kernel in d = 3
def q_KKSZ(Delta):        # 1904.05905 (5.132) at d = 3
    return (mp.mpf(2)**(-8)*mp.pi**4.5*4*36*mp.gamma(Delta)*mp.gamma(Delta - mp.mpf(1)/2)
            * mp.rgamma(3 - Delta/mp.mpf(2))**2*mp.rgamma(2 + Delta/mp.mpf(2))**4*mp.rgamma((3 + Delta)/mp.mpf(2))**2)
def h_even_CMT(Delta):    # 1710.03199 (4.11)
    return (12*mp.sqrt(6)*mp.pi**2*mp.sqrt(mp.gamma(2*Delta - 1))*mp.rgamma(2 + Delta/mp.mpf(2))*mp.rgamma(Delta + 3)
            * mp.rgamma(3 - Delta/mp.mpf(2)))
def h_odd_CMT(Delta):     # 1710.03199 (4.10)
    return (12*mp.sqrt(6)*mp.pi**2*mp.sqrt(mp.gamma(2*Delta - 1))*mp.rgamma((Delta + 1)/mp.mpf(2))*mp.rgamma(Delta + 3)
            * mp.rgamma((7 - Delta)/mp.mpf(2)))
