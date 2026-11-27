"""Numbers for the eta/s and thermal-free-energy corollaries.  Inputs: letter.tex / macros.tex bounds."""
from math import sqrt, log, pi

zeta3 = 1.2020569031595942


def Lam(D):
    return sqrt((D - 2) * (D - 1))


def tau_star(D):
    return 160 / (9 * (D - 2)) / sqrt(1 - 49 / (2 * D + 9) ** 2)


def tau_an1(D):   # analytic all-gap bound (needs one convergence input)
    return 280 * sqrt(log((D + 2) / 3)) / (D - 2) ** 2


def tau_an2(D):   # with the further k=3 input
    return 420 * sqrt(log(D)) / (D - 2) ** 2.5


C_ETA = 3 / 40          # |4 pi eta/s - 1| = (3/40)|t4| + ... ;  (3/40) tau = 18 |ghat3|
C_F = 1 / 720           # |f/f_E - 1| = |t4|/720 + ...
D = 1000
L = Lam(D)
print(f"Delta_g = {D}: Lambda = {L:.4f}, Lambda^4 = {L**4:.4e}, Lambda^6 = {L**6:.4e}, Lambda^8 = {L**8:.4e}")
print(f"  tau_* = {tau_star(D):.6e}  ->  X_* = (3/40) tau_* = {C_ETA*tau_star(D):.6e}")
print(f"  tau <= 2.8914e-8 (certified)  ->  X_cubic = {C_ETA*2.8914e-8:.4e}  (2.89e-8 -> {C_ETA*2.89e-8:.4e})")
for Q, lab in ((92.8, 'no scalar < 0.46 Lambda'), (164.5, 'with light scalars')):
    G4 = Q / L ** 6
    print(f"  quartic: g4/8piG L^6 <= {G4:.4e} ({lab});  216 G4 = {216*G4:.3e};  108 G4 = {108*G4:.3e}")
print(f"  type II reference: g4/8piG = 2 zeta3/(ML)^6 with ML ~ Lambda: 108 G4 = {108*2*zeta3/L**6:.3e}")
X_tot = C_ETA * 2.8914e-8 + 216 * 92.8 / L ** 6
print(f"  X_total (cubic + quartic, c1 <= G4/4) = {X_tot:.4e}")
print(f"  t4^2 term: (17/3840) tau^2 <= {17/3840*(2.8914e-8)**2:.1e}")
for c in (1, 20, 100):
    print(f"  remainder, natural coefficient c = {c}: I2 I3 ~ 1944 c/Lambda^8 = {1944*c/L**8:.1e}; "
          f"I3^2 ~ 6480 c/Lambda^10 = {6480*c/L**10:.1e}; I2^3 ~ 31104 c/Lambda^10 = {31104*c/L**10:.1e}")
print(f"  thermal free energy |f/f_E - 1| <= tau/720 + 4 G4 = {2.8914e-8/720 + 4*92.8/L**6:.3e}; "
      f"unconditional tau_*/720 = {tau_star(D)/720:.3e}")
print("\nGap scan:  Delta  X_*=(3/40)tau_*    (4/3)/(Delta-2)   X_an1        X_an2       X_ownpos (18*B/Lambda^4)")
B = {1000: 120, 2000: 128, 5000: 139}
for Dg in (10, 30, 100, 300, 1000, 2000, 5000, 1e4, 2e4, 1e5, 2e5, 1e6):
    own = f"{18*B[Dg]/Lam(Dg)**4:.3e}" if Dg in B else "   -"
    print(f"  {Dg:>8g}  {C_ETA*tau_star(Dg):.4e}   {4/(3*(Dg-2)):.4e}     {C_ETA*tau_an1(Dg):.3e}   "
          f"{C_ETA*tau_an2(Dg):.3e}   {own}")
print("\nReverse reading: deviation delta = |4 pi eta/s - 1| needs tau >= (40/3) delta, i.e. a higher-spin single trace"
      " below twist g_*((40/3)delta) ~ 2 + 17.8*3/(40 delta) = 2 + 1.335/delta")
for dl in (0.3, 0.1, 0.01, 1e-3):
    print(f"  delta = {dl}: Delta_g <~ {2 + 17.8 * 3 / (40 * dl):.1f}")
