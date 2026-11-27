"""S2 [computed]. Large-J decay of the stripped detector kernels at the light double-twist twists.

|a_{O,0,0}|^2 = |lam_pp Kpp + i lam_odd Kodd|^2/C_T, |a_{O,4,0}|^2 = lam_pm^2 Kpm^2/C_T (code/kernels/closed_forms.py).
Kpp, Kpm carry the zero factor sin^2(pi h), Kodd, KoddJ carry cos^2(pi h), h = (Delta-J)/2.  We divide it out
("stripped kernel") and measure, at FIXED h on each light family (parity even: h = n+1; odd: h = n+3/2), the exponent of
     k_X(h,J) := 4^{-Delta} 2^{J} |K_X|^2 / zero factor      as J -> infinity.
Expected from the factorised forms (angproofs/ap6_envelopes.py): -13/2 (pp, pm) and -15/2 (odd, oddJ).
The factorised forms are imported read-only and re-checked against the closed forms at off-zero points."""
import mpmath as mp
from common_a2 import Tee
import ap6_envelopes as ap6
mp.mp.dps = 50
log = Tee('s2_kernel_decay.txt')

def stripped(name, h, J):
    D = J + 2*h
    zf = mp.sin(mp.pi*h)**2 if name in ('pp', 'pm') else mp.cos(mp.pi*h)**2
    return ap6.K2_fact(name, D, J)/zf*mp.mpf(4)**(-D)*mp.mpf(2)**J

log("check: factorised = closed form (rel. dev. at h off the zeros)")
for name, Jv in (('pp', 40), ('pm', 40), ('odd', 40), ('oddJ', 41)):
    D = Jv + 2*mp.mpf('1.37')
    log(f"   {name:5s} J={Jv}: {mp.nstr(abs(ap6.K2_fact(name, D, Jv)/ap6.K2_closed(name, D, Jv) - 1), 3)}")
eps = mp.mpf('1e-20')
log("\nlocal exponent d log k_X / d log J at fixed h (J -> 2J):")
log(" X     h        J=100     J=1000    J=10^4    J=10^5")
for name, hs in (('pp', [1, 2, 3, 4, 11, 101]), ('pm', [1, 2, 4, 101]),
                 ('odd', [1.5, 2.5, 3.5, 4.5, 101.5]), ('oddJ', [1.5, 2.5, 4.5])):
    for h0 in hs:
        h = mp.mpf(h0) + eps
        row = []
        for J0 in (100, 1000, 10**4, 10**5):
            Ja = J0 + (1 if name == 'oddJ' else 0)
            Jb = 2*J0 + (1 if name == 'oddJ' else 0)
            row.append(mp.log(stripped(name, h, Jb)/stripped(name, h, Ja))/mp.log(mp.mpf(Jb)/Ja))
        log(f" {name:5s} {float(h0):6.1f}  " + "  ".join(f"{float(x):+.5f}" for x in row))
log("=> k_pp, k_pm ~ J^{-13/2};  k_odd, k_oddJ ~ J^{-15/2}  at every fixed light twist (incl. h = 1, 3/2).")
