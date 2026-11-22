"""
dgscan/bounds.py (DG_SCAN.md §1): explicit Delta_g dependence of the bound VALUES of Results D and Q, with the functionals held fixed.

Nothing here touches per-state positivity: these are the light-side numbers that the implication "positivity => bound" produces at a given
Lambda = sqrt((Dg-2)(Dg-1)).  Inputs (read-only):
  * NF02 weights psi2, psi3 (code/robust/functionals_noforward.json via omegau/nf_common.NF, k(0) = 1)   -> Result D
  * Q3-full-delta0.2 psi2 (code/resultq/functionals_q3.json)                                              -> Result Q upper end
Formulas (RESULT_D_STATEMENT §3, RESULT_Q_STATEMENT §0, REVIEW_UNIFIED item 34, H7_I14 §23.2, h3_H17_model.py):
  16 L_A(Lam)   = 16 int_0^1 psi2^NF(t) dt/(t^2 + 9/(4 Lam^2))
  cubic bound   = 240 sqrt(16 L_A/(k(0)(1+dk))) / Lam^4            (headline: k(0)(1+dk) -> 1)
  2 L_AdS(Lam)  = 2 int_0^1 psi2^Q(p) dp/(p^2 + 9/(4 Lam^2)) = 17.241 log(2 Lam/3) + B_eff(Lam)
  c1 (cubic)    : relative shift of the bound per unit c1 = (1/2) 16 Lam^2 t4* int psi2^NF / (16 L_A),  t4* = cubic bound
  c1 (quartic)  : |eps_c1| = 2 |c1 t4| Lam^2 int psi2^Q;  with t4 <= cubic bound: <= 2 int psi2^Q * 240 sqrt(16 L_A)/Lam^2 * |c1|
  dk(0) (EE)    : first-order symbol X - psi2 = Lam^-2 [psi2''/2 - (23/2)(psi2/p)'] (H7_I14 §23.2)  =>  Lam^2 dk(0) = -168 int psi2 p^4
  H17 band      : |dk_H17| = 36 sum_k |c_k mu_k| D_k / Lam^2 (h3 heuristic)  -> bound factor (1 -+ dk)^(-1/2)
Asymptotics: for L(a) = int_0^1 psi/(p^2 + a^2), psi = c1 p + c2 p^2 + c3 p^3 + R(p), a = 3/(2 Lam):
  L = c1 log(1/a) + [c2 + c3/2 + int R/p^2] - (pi/2) c2 a + a^2 [c1/2 + c2 - c3 log(1/a) - int R/p^4] + O(a^3)
Usage: python3 bounds.py  -> stdout + cache/bounds.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
import mpmath as mp
from nf_common import NF

mp.mp.dps = 40
OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def lam_of(Dg):
    return mp.sqrt((mp.mpf(Dg) - 2)*(mp.mpf(Dg) - 1))


# ------------------------------------------------------------------------------------------------ exact polynomials (monomial, low first)
F = NF()
PSI2_NF = [mp.mpf(c.numerator)/c.denominator for c in F.s2]            # exact rational NF02 psi2 / k0  (k(0) = 1)
PSI3_NF = [mp.mpf(c.numerator)/c.denominator for c in F.s3]
Qd = json.load(open(os.path.join(HERE, '..', 'resultq', 'functionals_q3.json')))['Q3-full-delta0.2']
PSI2_Q = [mp.mpf(c) for c in Qd['psi2']]


def poly(c, p):
    return mp.polyval(c[::-1], p)


def mom(c, j):
    return sum(ci/(i + j + 1) for i, ci in enumerate(c))


def L_exact(c, Lam):
    """int_0^1 psi(p)/(p^2 + a^2) dp in closed form (term by term, exact for a polynomial)."""
    a = mp.mpf(3)/(2*Lam)
    tot = mp.mpf(0)
    # I_n = int_0^1 p^n/(p^2 + a^2): I_0 = atan(1/a)/a, I_1 = log(1 + 1/a^2)/2, I_n = 1/(n-1) - a^2 I_{n-2}
    I = [mp.atan(1/a)/a, mp.log(1 + 1/a**2)/2]
    for n in range(2, len(c)):
        I.append(mp.mpf(1)/(n - 1) - a*a*I[n - 2])
    return sum(ci*I[n] for n, ci in enumerate(c))


def L_quad(c, Lam):
    a = mp.mpf(3)/(2*Lam)
    return mp.quad(lambda p: poly(c, p)/(p*p + a*a), [0, a, 10*a, 0.1, 1])


def asym(c):
    """coefficients of L = c1 log(1/a) + K0 + K1 a + a^2 (K2L log(1/a) + K2) + O(a^3)."""
    c1, c2, c3 = c[1], c[2], c[3]
    R = [mp.mpf(0)]*4 + list(c[4:])
    intR_p2 = sum(R[n]/(n - 1) for n in range(4, len(R)))          # int_0^1 R/p^2
    intR_p4 = sum(R[n]/(n - 3) for n in range(5, len(R))) if len(R) > 4 else 0
    # int R/p^4 diverges only if R has a p^4 term at p -> 0?  R/p^4 -> c4 (finite): int_0^1 c4 dp = c4  (n = 4 term: 1/(4-3) = 1)
    intR_p4 = sum(R[n]/(n - 3) for n in range(4, len(R)))
    K0 = c2 + c3/2 + intR_p2
    K1 = -mp.pi/2*c2
    K2L = -c3
    K2 = c1/2 + c2 - intR_p4
    return c1, K0, K1, K2L, K2


def L_asym(c, Lam, order=2):
    a = mp.mpf(3)/(2*Lam)
    c1, K0, K1, K2L, K2 = asym(c)
    v = c1*mp.log(1/a) + K0
    if order >= 1:
        v += K1*a
    if order >= 2:
        v += a*a*(K2L*mp.log(1/a) + K2)
    return v


GRID = (100, 150, 200, 250, 300, 500, 1000, 2000, 5000, 10000)
INT2_NF = mom(PSI2_NF, 0)
INT2_Q = mom(PSI2_Q, 0)
MU4 = mom(PSI2_NF, 4)
# H17 heuristic constant (h3_H17_model.py, route A'): 36 * (|2 mu(psi2, 4)| * 3 + |-2 mu(psi3, 2)| * 2)
CH17 = 36*(abs(2*mom(PSI2_NF, 4))*3 + abs(-2*mom(PSI3_NF, 2))*2)


def main():
    pr("# dgscan/bounds.py: bound values of Results D and Q vs Delta_g, functionals fixed (NF02 for D; Q3-full-delta0.2 for Q)")
    pr(f"# checks: k(0) = -2 mu(psi2,6) + 2 mu(psi3,4) = {mp.nstr(-2*mom(PSI2_NF, 6) + 2*mom(PSI3_NF, 4), 12)} (must be 1)")
    A_D = 16*PSI2_NF[1]; A_Q = 2*PSI2_Q[1]
    pr(f"#   A_D = 16 psi2'(0) = {mp.nstr(A_D, 10)} (doc 3022.34);  A_Q = 2 psi2^Q'(0) = {mp.nstr(A_Q, 10)} (doc 17.241)")
    pr(f"#   int psi2^NF = {mp.nstr(INT2_NF, 8)} (doc 10.386);  int psi2^Q = {mp.nstr(INT2_Q, 8)} (doc 1.016);  mu4 = int psi2^NF p^4 = {mp.nstr(MU4, 8)}")
    pr(f"#   H17 heuristic: |dk_H17| = {mp.nstr(CH17, 6)}/Lam^2;  first-order-symbol dk(0)_EE = -168 mu4/Lam^2 = {mp.nstr(-168*MU4, 6)}/Lam^2")
    # asymptotic coefficients
    for name, c, fac in (("16 L_A (D)", PSI2_NF, 16), ("2 L_AdS (Q)", PSI2_Q, 2)):
        c1, K0, K1, K2L, K2 = asym(c)
        pr(f"# asymptotics of {name} = A log(2Lam/3) + B_inf + b1/Lam + (b2L log(2Lam/3) + b2)/Lam^2 + O(Lam^-3):")
        pr(f"#   A = {mp.nstr(fac*c1, 10)}, B_inf = {mp.nstr(fac*K0, 10)}, b1 = {mp.nstr(fac*K1*1.5, 8)}, b2L = {mp.nstr(fac*K2L*2.25, 8)},"
           f" b2 = {mp.nstr(fac*K2*2.25, 8)}   [log(Lam/Lam0) form: Lam0 = 1.5 exp(-B_inf/A) = {mp.nstr(1.5*mp.exp(-K0/c1), 6)}]")
        pr(f"#   (B_inf + A log(2/3) = {mp.nstr(fac*K0 + fac*c1*mp.log(mp.mpf(2)/3), 9)}: the 'A log Lam + B' form)")
    pr(f"# Q's stored flat B = {Qd['B']:.6f}  (compare B_inf of 2 L_AdS)")
    pr("")
    pr("| Dg | Lam | 16L_A exact | (quad check) | asym O(a^2) | cubic bound (k=1) | dk(0) factor | H17 band | per unit c1 (D) |"
       " 2L_AdS exact | B_eff = 2L - A log(2Lam/3) | asym B_eff (O(1/Lam)) | (O(a^2)) | eps_c1 bound (Q, with D's t4) | ABJM 2zeta3/2L |")
    pr("|" + "---|"*15)
    rows = []
    for Dg in GRID:
        Lam = lam_of(Dg)
        L16 = 16*L_exact(PSI2_NF, Lam)
        L16q = 16*L_quad(PSI2_NF, Lam)
        L16a = 16*L_asym(PSI2_NF, Lam, 2)
        bD = 240*mp.sqrt(L16)/Lam**4
        dk = -168*MU4/Lam**2
        dkf = (1 + dk)**-0.5
        h17 = CH17/Lam**2
        h17b = ((1 - h17)**-0.5 - 1)
        c1D = mp.mpf(1)/2*16*Lam**2*bD*INT2_NF/L16
        L2 = 2*L_exact(PSI2_Q, Lam)
        Beff = L2 - A_Q*mp.log(2*Lam/3)
        Beff1 = 2*L_asym(PSI2_Q, Lam, 1) - A_Q*mp.log(2*Lam/3)
        Beff2 = 2*L_asym(PSI2_Q, Lam, 2) - A_Q*mp.log(2*Lam/3)
        epsQ = 2*INT2_Q*240*mp.sqrt(L16)/Lam**2
        rows.append(dict(Dg=Dg, Lam=float(Lam), L16=float(L16), bD=float(bD), dk=float(dk), h17=float(h17), c1D=float(c1D), L2=float(L2),
                         Beff=float(Beff), epsQ=float(epsQ)))
        pr(f"| {Dg} | {mp.nstr(Lam, 6)} | {mp.nstr(L16, 8)} | {mp.nstr(L16q - L16, 2)} | {mp.nstr(L16a, 8)} | {mp.nstr(bD, 5)} |"
           f" x{mp.nstr(dkf, 6)} (dk {mp.nstr(dk, 4)}) | +-{mp.nstr(100*h17b, 3)} % | -+{mp.nstr(100*c1D, 3)} % | {mp.nstr(L2, 7)} |"
           f" {mp.nstr(Beff, 7)} | {mp.nstr(Beff1, 7)} | {mp.nstr(Beff2, 7)} | {mp.nstr(epsQ, 4)} |c1| | {mp.nstr(100*2*mp.zeta(3)/L2, 3)} % |")
    pr("")
    # documented values
    pr("# reproduction of the stated numbers: 300 / 1000")
    for r in rows:
        if r['Dg'] in (100, 300, 1000):
            pr(f"  Dg={r['Dg']}: 16L_A={r['L16']:.2f} (doc 7508.51/10734.61/14339.59)  bound={r['bD']:.5e} (doc 2.2094e-4/3.1321e-6/2.8913e-8)"
               f"  2L_AdS={r['L2']:.3f} (doc 53.107/72.041/92.814)  B_eff={r['Beff']:.4f} (doc -19.0378/-19.2197/-19.2647)"
               f"  c1(D)={100*r['c1D']:.4f}% (doc 0.22/0.017)  eps_c1(Q)={r['epsQ']:.4f}|c1| (doc 0.57/0.059)  dk={r['dk']:.4e} (EE doc -7.245e-3/-7.666e-4/-6.830e-5)")
    # leading-order asymptotic form of the cubic bound in Delta_g
    pr("")
    pr("# closed asymptotic forms (checked against exact):")
    c1, K0, K1, K2L, K2 = asym(PSI2_NF)
    A = 16*c1; B = 16*K0; b1 = 16*K1*1.5
    L0D = 1.5*mp.exp(-K0/c1)
    c1q, K0q, K1q, _, _ = asym(PSI2_Q)
    L0Q = 1.5*mp.exp(-K0q/c1q)
    for Dg in GRID:
        Lam = lam_of(Dg)
        ex = 240*mp.sqrt(16*L_exact(PSI2_NF, Lam))/Lam**4
        f0 = 240*mp.sqrt(A*mp.log(Lam/L0D))/Lam**4
        f1 = 240*mp.sqrt(A*mp.log(Lam/L0D) + b1/Lam)/Lam**4
        g0 = 240*mp.sqrt(A*mp.log(Dg/L0D))/mp.mpf(Dg)**4*(1 + mp.mpf(6)/Dg)
        q_ex = 2*L_exact(PSI2_Q, Lam)
        q0 = 2*c1q*mp.log(Lam/L0Q)
        q1 = q0 + 2*K1q*1.5/Lam
        pr(f"  Dg={Dg:5d}: cubic exact {mp.nstr(ex, 6)} | 240 sqrt(A log(Lam/{mp.nstr(L0D, 4)}))/Lam^4: rel {mp.nstr(f0/ex - 1, 3)} | + b1/Lam: rel {mp.nstr(f1/ex - 1, 3)}"
           f" | in Dg (240 sqrt(A log(Dg/Lam0))(1+6/Dg)/Dg^4): rel {mp.nstr(g0/ex - 1, 3)} || quartic exact {mp.nstr(q_ex, 7)} | A log(Lam/{mp.nstr(L0Q, 4)}): "
           f"{mp.nstr(q0 - q_ex, 3)} | + b1/Lam: {mp.nstr(q1 - q_ex, 3)}")
    open(os.path.join(HERE, 'cache', 'bounds.txt'), 'w').write("\n".join(OUT) + "\n")
    json.dump(rows, open(os.path.join(HERE, 'cache', 'bounds.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
