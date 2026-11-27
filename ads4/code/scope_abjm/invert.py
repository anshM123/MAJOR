r"""Invert the quartic (R^4) bound of the Letter against the ABJM (type IIA) value of g4.

INPUTS (all cited, see README.md):
  paper  : B(Lambda) = upper end of  g4~ Lambda^6 / 8 pi G   [supplement Table tab:dg, eq:quarticL, eq:repair]
  CHZ    : nu := L^4/l_s^4 = 2 pi^2 (lambda - 1/24) + ...      [Chester-Hansen-Zhong 2412.08689, eq.(1.1); Bergman-Hirano shift]
  BCP    : flat-space R^4 coefficient  g4/(8 pi G) = zeta(3) alpha'^3/32  (A_IIA = A_SG(1 + l_s^6 zeta(3)/32 stu)),
           fixed in the 't Hooft limit by localization  [Binder-Chester-Pufu 1906.07195, eqs.(1.2),(1.4),(4.31)]
           => in AdS units (L=1):  g4/(8 pi G) = zeta(3)/(32 nu^{3/2}) = 2 zeta(3)/(M L)^6,  M L = 2 nu^{1/4}
  CHZ    : first massive level: Delta(Delta-3) = R^2 m^2 = 4 delta R^2/alpha', delta=1 => twist = 2 nu^{1/4} + O(1)
           [2412.08689 eqs.(3.3),(3.4),(1.6)]
COMPUTED here: Lambda_max(lambda), Delta_max, Y = Delta_max/Delta_true, lambda_*(Delta_g) for the fixed-gap inversion.
"""
import mpmath as mp
mp.mp.dps = 30
z3 = mp.zeta(3); twoz3 = 2*z3

# ---- the paper's bound -------------------------------------------------------------------------------
tab = {  # Delta_g : (Lambda, 16 L_A, 2 L_AdS)   [supplement tab:dg]
    300: (298.50, 10734.61, 72.041), 500: (498.50, 12259.23, 80.855), 1000: (998.50, 14339.59, 92.814),
    2000: (1998.5, 16427.19, 104.769), 5000: (4998.5, 19192.12, 120.571)}
def B_asym(Lam):   # eq:quarticL, all terms shown (accurate to 1e-3 absolute at 300, 1000)
    x = mp.log(2*Lam/3)
    return 17.2407*x - 19.2782 + 10.314/Lam + (860.3*x - 2508.6)/Lam**2
def B_repair(D):   # eq:repair: 2L_AdS + 2 eps L_A, eps = 0.040  ->  2L_AdS + 0.005*(16 L_A)
    Lam, L16, B = tab[D]; return B + 0.005*L16
Lam_of = lambda D: mp.sqrt((D-2)*(D-1))
D_of = lambda Lam: mp.mpf(3)/2 + mp.sqrt(Lam**2 + mp.mpf(1)/4)

print("1. The bound B(Delta_g) [paper values; 2b at 500/2000/5000 recomputed from tab:dg, NOT stated in the paper]")
for D, (Lam, L16, B) in tab.items():
    print(f"   Delta_g={D:5d}  Lambda={Lam:7.2f}  2a (scalar cond.) B={B:8.3f}  asym={float(B_asym(Lam)):8.3f}  2b (no scalar cond.)={float(B_repair(D)):8.2f}"
          f"   Y_2a=(B/2z3)^(1/6)={float((B/twoz3)**(mp.mpf(1)/6)):.4f}  Y_2b={float((B_repair(D)/twoz3)**(mp.mpf(1)/6)):.4f}")
print(f"   2 zeta(3) = {float(twoz3):.5f};  2z3/92.8135 = {float(twoz3/92.8135):.4f};  (92.8135/2z3)^(1/6) = {float((92.8135/twoz3)**(mp.mpf(1)/6)):.4f}")
print(f"   10d flat window [0,3.000] (CHMRS-swamp): Y_10d = (3/2z3)^(1/6) = {float((3/twoz3)**(mp.mpf(1)/6)):.4f}")
print(f"   AdS5 / N=4 (CHMRS 2106.10274 eq.5.49, alpha/8piG <= 10.9/Delta_gap^6): Y_AdS5 = (10.9/2z3)^(1/6) = {float((10.9/twoz3)**(mp.mpf(1)/6)):.4f}")

# ---- ABJM dictionary -----------------------------------------------------------------------------------
nu = lambda lam, shift=True: 2*mp.pi**2*(lam - (mp.mpf(1)/24 if shift else 0))
ML = lambda lam, shift=True: 2*nu(lam, shift)**mp.mpf(0.25)
print("\n2. Dictionary  ML = 2 nu^{1/4} = 2^{5/4} pi^{1/2} lambda^{1/4} =", mp.nstr(2**mp.mpf(1.25)*mp.sqrt(mp.pi), 6), "lambda^{1/4}")
for D in (300, 1000):
    Lam = Lam_of(D); n = (Lam/2)**4; lam = n/(2*mp.pi**2) + mp.mpf(1)/24
    # alternative: twist of spin-4 = ML - 4 (paper's own estimate) equal to Delta_g
    lam2 = ((D+4)/mp.mpf(2))**4/(2*mp.pi**2)
    print(f"   Delta_g={D}: Lambda=ML  -> lambda = {mp.nstr(lam, 4)};   twist(J=4)=ML-4=Delta_g -> lambda = {mp.nstr(lam2, 4)}")
print("   g4/(8 pi G L^6) = zeta(3)/(32 nu^{3/2}) =", mp.nstr(z3/(32*(2*mp.pi**2)**1.5), 6), "* lambda^{-3/2}  (leading)")
print("   8 pi G/L^2 = 256/c_T;  c_T = (64 sqrt2/(3 pi)) N^2/sqrt(lambda) (leading) -> g4/L^8 = 8 zeta(3)/(c_T nu^{3/2})")

# ---- (A) fixed-gap inversion: only the gaps at which the paper evaluates/claims the quartic bound ------
print("\n3A. Fixed-gap inversion (the only form the paper's claims support).")
print("    If twist gap >= Delta_g then g4 Lambda_g^6/8piG <= B(Delta_g); ABJM: g4 Lambda^6/8piG = 2z3 (Lambda/ML)^6.")
print("    Violated iff ML < Lambda_g (2z3/B)^{1/6}  <=>  lambda < lambda_*.  Then ABJM has a spin>=4 single trace in TxT with twist < Delta_g.")
for D, label, B in [(1000, "2a, claimed gap (HLJclosed=false)", 92.8135), (1000, "2b, no scalar condition", 164.51),
                    (300, "2a, HLJclosed=true (one hypothesis open)", 72.0406), (300, "2b at 300", 125.71)]:
    Lam = Lam_of(D); MLs = Lam*(twoz3/B)**(mp.mpf(1)/6)
    lam_s = (MLs/2)**4/(2*mp.pi**2) + mp.mpf(1)/24
    print(f"    Delta_g={D:4d} [{label:42s}] B={B:8.3f}: ML_* = {mp.nstr(MLs,5):>7}, lambda_* = {mp.nstr(lam_s,4):>9};"
          f"  true twist at lambda_* ~ ML_*-4 = {mp.nstr(MLs-4,5)};  ratio Delta_g/true = {mp.nstr(D/(MLs-4),4)}  (Lambda/ML = {mp.nstr(Lam/MLs,4)})")

# ---- (B) running-gap inversion: needs positivity of Omega^Q at the gap Delta_max(lambda) itself ----------
print("\n3B. Running-gap inversion (needs (Pos) for Omega^Q at Delta_max(lambda); established nowhere, examined only at 300, 1000).")
print("    Solve 2z3 (Lambda/ML)^6 = B(Lambda) with B = eq:quarticL (2a).")
def Lam_max(lam):
    m = ML(lam); f = lambda L: twoz3*(L/m)**6 - B_asym(L)
    return mp.findroot(f, 1.85*m)
print("      lambda        ML      Lambda_max   Delta_max   true twist~ML-4   Y=Lambda_max/ML   Delta_max/(ML-4)   Delta_max/lambda^{1/4}")
for lam in [1e4, 1e5, 1e6, 1e7, 2.75e8, 1e9, 3.15e9, 1e11, 1e14, 1e20]:
    lam = mp.mpf(lam); m = ML(lam); Lm = Lam_max(lam); Dm = D_of(Lm)
    print(f"    {mp.nstr(lam,3):>9} {mp.nstr(m,6):>9} {mp.nstr(Lm,7):>11} {mp.nstr(Dm,7):>11} {mp.nstr(m-4,6):>13}"
          f"        {mp.nstr(Lm/m,5):>8}         {mp.nstr(Dm/(m-4),5):>8}          {mp.nstr(Dm/lam**0.25,5)}")
print("    Asymptotics: Y(lambda) = [ (17.24 ln(2 Y ML/3) - 19.28)/(2 zeta3) ]^{1/6}  ~ (ln lambda)^{1/6}; "
      "Delta_max ~ 4.2156 Y lambda^{1/4}")
# sensitivity: Y as function of B
print("\n4. Sensitivity of Y to the bound: Y = (B/2z3)^{1/6}")
for B in [3.0, 10.9, 30, 50, 72.04, 92.81, 164.5, 500, 1000, 1e4]:
    print(f"    B={B:8.2f}  Y={float((B/twoz3)**(1/6.)):.3f}")
# a bound that was wrong by a factor F changes Y by F^{1/6}
print("    A factor-2 change in B moves Y by 2^{1/6} =", round(2**(1/6.),4), "; factor 10 by", round(10**(1/6.),4))

# ---- M-theory regime (for contrast only) ------------------------------------------------------------
print("\n5. M-theory regime (k fixed, N->inf), contrast only: 11d R^4 = l11^6 stu/(3*2^7) (CPY 1804.00949); no string scale.")
print("   g4 Lambda^6/8piG ~ (Lambda l11/L)^6/384 = B  =>  Lambda_max ~ (384 B)^{1/6} L/l11 =",
      mp.nstr((384*92.81)**(1/6.),4), "L/l11 at B=92.8: Planckian; A2 (tree-level, c_T->inf at fixed gap) fails there.")
