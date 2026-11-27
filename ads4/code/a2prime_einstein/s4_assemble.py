"""S4 [computed given S1-S3; constants set to 1]. The A2' detector sum on tree-level Einstein data.

Per light state O = [TT]_{n,J} (paper's normalisation (3.1), structure X, ladder b):
    |a_{O,b,0}|^2 = C_T * lamhat_X^2 * |K_X(Delta,J)|^2 ,   Delta = J + 2h,  h = h_n + gamma_{n,J}/2,
    lamhat_X^2   = A_n 4^{-Delta} 2^{J} s^{p_X}          (GFF, unit-normalised T; p = 1/2 even, 3/2 odd: S3),
    K_X           exact closed forms (code/kernels), which carry sin^2(pi h) / cos^2(pi h)  (S2),
    gamma_{n,J}   = kappa_n/(c_T J)                       (tree level, S1).
The weighted sum is  E_q(c_T) = sum_{n, J} (1+J)^q log(2+J) w_{O,b} |a_{O,b,0}|^2,  w = [T_Delta - b(b-1)]/3 + 3/2,
the paper's being q = 5/2.  Unknown O(1) constants (A_n, kappa_n) are set to 1: only exponents are claimed.
Also included: the twist-2 'N-structure' of S3(c), amplitude ~ B_J/c_T with no zero factor, for B_J = 1/J (estimated
tree-level decay) and for B_J = 1 (worst case: no decay)."""
import mpmath as mp
from common_a2 import Tee
import closed_forms as cf
mp.mp.dps = 25
log = Tee('s4_assemble.txt')

def a2(name, n_h, J, cT, kappa=1):
    gam = mp.mpf(kappa)/(cT*J)
    h = n_h + gam/2
    D = J + 2*h; s = (D + J)/2
    p = mp.mpf(1)/2 if name in ('pp', 'pm') else mp.mpf(3)/2
    lam2 = mp.mpf(4)**(-D)*mp.mpf(2)**J*s**p
    K2 = {'pp': cf.Kpp2, 'pm': cf.Kpm2, 'odd': cf.Kodd2}[name](D, J)
    return cT*lam2*K2, D

def a2_N(J, cT, B):
    """twist-2 N-structure: |a_4|^2 ~ C_T (B/C_T)^2 x (stripped pm envelope at h = 1)"""
    D = J + 2; s = (D + J)/2
    eps = mp.mpf('1e-12')
    strip = cf.Kpm2(J + 2 + 2*eps, J)/mp.sin(mp.pi*(1 + eps))**2
    return cT*(B/cT)**2*mp.mpf(4)**(-D)*mp.mpf(2)**J*s**mp.mpf(0.5)*strip, D

def weight(D, b):
    return ((D - 1)*(D - 2) - b*(b - 1))/3 + mp.mpf(3)/2

FAM = [('pp', 0, 1), ('pp', 0, 2), ('pp', 0, 3), ('pp', 0, 4), ('pm', 4, 1), ('pm', 4, 2), ('odd', 0, mp.mpf(3)/2),
       ('odd', 0, mp.mpf(5)/2)]

def E(cT, q, Jmax, kappa_power=0, extraN=None):
    tot = mp.mpf(0); tail = {}
    for name, b, nh in FAM:
        kap = (nh)**kappa_power
        for Jv in range(4, Jmax + 1, 2):
            v, D = a2(name, nh, Jv, cT, kap)
            t = (1 + Jv)**q*mp.log(2 + Jv)*weight(D, b)*v
            tot += t
            if Jv in (Jmax//2, Jmax): tail[(name, nh, Jv)] = t
    if extraN is not None:
        for Jv in range(4, Jmax + 1, 2):
            v, D = a2_N(Jv, cT, extraN(Jv))
            tot += (1 + Jv)**q*mp.log(2 + Jv)*weight(D, 4)*v
    return tot, tail

log("(i) per-state summand exponent at fixed n (q = 5/2): d log(term)/d log J between J and 2J")
for name, b, nh in FAM[:1] + FAM[4:5] + FAM[6:7]:
    for cT in (mp.mpf(10)**6,):
        row = []
        for J0 in (100, 1000, 10000):
            t = []
            for Jv in (J0, 2*J0):
                v, D = a2(name, nh, Jv, cT)
                t.append((1 + Jv)**2.5*mp.log(2 + Jv)*weight(D, b)*v)
            row.append(mp.log(t[1]/t[0])/mp.log(2))
        log(f"   {name:4s} h={float(nh)}: " + ", ".join(f"{float(x):+.3f}" for x in row)
            + "   (expected -7/2 + log: per state J^{5/2} log * J^2 * C_T J^{-6} * gamma^2)")

log("\n(ii) E_{5/2}(c_T) and its c_T scaling (families n = 0..3 even, both ladders; odd n = 0, 1)")
for Jmax in (2000, 8000):
    vals = []
    for cT in (mp.mpf(10)**4, mp.mpf(10)**6, mp.mpf(10)**8):
        e, _ = E(cT, 2.5, Jmax)
        vals.append(e)
        log(f"   Jmax={Jmax:5d}  c_T=1e{int(mp.log10(cT))}:  E = {mp.nstr(e, 8)}   c_T*E = {mp.nstr(cT*e, 8)}")
log("   => c_T*E converges in Jmax and is c_T-independent: E_light = O(1/c_T) -> 0 (A2' holds).")

log("\n(iii) with the twist-2 N-structure added (c_T = 1e6)")
for lab, Bf in (("B_J = 1/J (estimated)", lambda J: mp.mpf(1)/J), ("B_J = 1 (no decay, worst case)", lambda J: mp.mpf(1))):
    r = [E(mp.mpf(10)**6, 2.5, Jm, extraN=Bf)[0]*mp.mpf(10)**6 for Jm in (2000, 8000)]
    log(f"   {lab:32s}: c_T*E (Jmax 2000, 8000) = {mp.nstr(r[0], 8)}, {mp.nstr(r[1], 8)}")

log("\n(iv) how much spin weight is available: summand exponent for weight (1+J)^q, n = 0 pp, c_T = 1e6")
for q in (2.5, 4.0, 4.8, 5.0, 5.5):
    t = []
    for Jv in (10000, 20000):
        v, D = a2('pp', 1, Jv, mp.mpf(10)**6)
        t.append((1 + Jv)**q*weight(D, 0)*v)          # log(2+J) left out of the slope; it multiplies the result
    sl = mp.log(t[1]/t[0])/mp.log(2)
    log(f"   q = {q}: summand ~ J^{float(sl):+.3f} x log J  -> " + ("converges" if sl < -1.02 else "DIVERGES"))
log("   with gamma ~ 1/J: converges iff q < 5;  using only sup_J c_T gamma^2 -> 0 (no decay of gamma): iff q < 3.")

log("\n(v) eikonal n-growth kappa_n ~ n^3 (h_n^3) included: c_T*E at c_T = 1e6, Jmax = 4000")
e0, _ = E(mp.mpf(10)**6, 2.5, 4000, kappa_power=0)
e3, _ = E(mp.mpf(10)**6, 2.5, 4000, kappa_power=3)
log(f"   kappa = 1: {mp.nstr(e0*mp.mpf(10)**6, 6)};  kappa = h^3: {mp.nstr(e3*mp.mpf(10)**6, 6)}  (finite: n-sum is finite at fixed Delta_g)")
