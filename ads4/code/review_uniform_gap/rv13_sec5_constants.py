"""rv13_sec5_constants.py -- numerical/exact checks of Sec.5 elementary constants: psi2^(j) = (1-p)^(5-j) Q_j with ||Q_j||_1 < 2^40;
|(log F)^(j)| < 100 on the real axis (j=1,2,3); |F^(j)|/F <= 2^21."""
import sympy as sp, mpmath as mp
from rv_common import PSI2, Fz
p = sp.symbols('p')
psi = sum(sp.Rational(c.numerator, c.denominator)*p**i for i, c in enumerate(PSI2))
for j in range(4):
    dj = sp.diff(psi, p, j)
    Qj, rem = sp.div(sp.Poly(dj, p), sp.Poly((1-p)**(5-j), p))
    l1 = sum(abs(c) for c in Qj.all_coeffs())
    print(f"j={j}: remainder zero: {rem.is_zero}; ||Q_j||_1 = {float(l1):.4e} (2^40 = {2.0**40:.4e}) -> {l1 < 2**40}")
mp.mp.dps = 30
logF = lambda x: mp.log(mp.re(Fz(x)))
mx = [0, 0, 0, 0]; arg = [None]*4
for i in range(0, 4001):
    x = mp.mpf(i)/100 if i <= 2000 else mp.mpf(20) + (i-2000)*5
    for j in (1, 2, 3):
        v = abs(mp.diff(logF, x, j))
        if v > mx[j]:
            mx[j], arg[j] = v, x
print("max |(log F)^(j)| on x in [0,10020] for j=1,2,3:", [(j, mp.nstr(mx[j], 5), mp.nstr(arg[j], 4)) for j in (1, 2, 3)], "(claim < 100)")
# |F^(j)|/F: via Faa di Bruno from log-derivatives, bounded by (1+100)^3 < 2^21 trivially; direct check
mxr = 0
for i in range(0, 2001, 5):
    x = mp.mpf(i)/100
    F0 = mp.re(Fz(x))
    for j in (1, 2, 3):
        mxr = max(mxr, abs(mp.diff(lambda y: mp.re(Fz(y)), x, j))/F0)
print("max |F^(j)|/F on [0,20], j<=3:", mp.nstr(mxr, 5), "(claim <= 2^21 =", 2**21, ")")
