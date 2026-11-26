"""rv6_compare_fpH.py -- check (3a): the note's closed form vs the paper's Eq. (fp-H), the latter evaluated with the independent
xcheck_ladder implementation (xl_common.H_ang, read-only import) fed with moments Z_k, Y_k of the NEW family chi = B q
(rv5 table, uncut; also with theta(L+2000,10)), and Q_k exact.  Also validates the g[k,p] recurrence against mpmath."""
import sys, os, time, math
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'xcheck_ladder')))
import xl_common as xl
from flint import arb, fmpq, ctx
from rv_closed import gpolys, closed_form, AG
from rv_common import PSI2, mom

# --- recurrence vs Taylor coefficients of (1-z)^-p P_{-1/2+i nu}((1+z)/(1-z))
mp.mp.dps = 40
worst = 0
for nu in [mp.mpf('0.7'), mp.mpf('3.1')]:
    for p in (3, 4):
        fz = lambda zz: (1-zz)**(-p)*mp.legenp(-mp.mpf(1)/2+1j*nu, 0, (1+zz)/(1-zz), type=3)
        tay = mp.taylor(fz, mp.mpf('0'), 8, method='quad', radius=mp.mpf('0.5'))
        gp = gpolys(p, 8)
        for k in range(9):
            val = sum(mp.mpf(c.numerator)/c.denominator*nu**(2*i) for i, c in enumerate(gp[k]))
            worst = max(worst, abs(val-tay[k])/max(1, abs(tay[k])))
print("g[k,p] recurrence vs mpmath Taylor coefficients (k<=8, p=3,4, nu=0.7,3.1): max rel err", mp.nstr(worst, 3))

lam2 = int(sys.argv[1]); prec = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
states = [(int(a.split(',')[0]), int(a.split(',')[1]), a.split(',')[2]) for a in sys.argv[3:]]
d = np.load(f'out/chi_{lam2}.npz')
nu, wt, g, th = d['nu'], d['wt'], d['g'], d['theta']
L = float(d['L'])
Kmax = max(J for J, b, _ in states) + 2
t0 = time.time()
Sun = xl.power_moments(nu, wt*g, Kmax+6, prec)[0]
Scut = xl.power_moments(nu, wt*g*th, Kmax+6, prec)[0]
Z, Y = xl.moments_from_power(Sun, Kmax+1)
Zc, Yc = xl.moments_from_power(Scut, Kmax+1)
Q = xl.Q_exact(Kmax+1, lam2)
print(f"moments built in {time.time()-t0:.1f}s; Z_0/3 = {Z[0]/3}  (2 int psi2 = {float(2*mom(PSI2,0)):.12f})")
print("J b mass | H_fpH/w^3 (new chi, uncut) | closed form | rel diff | rel diff (theta-cut moments)")
for J, b, mspec in states:
    # mspec: 'tw' = twist edge at the gap with Lambda^2 = (Dg-1)(Dg-2) [only for lam2 = 997002], or a decimal m/Lambda
    if mspec == 'tw':
        Dg = int(round((3 + (1 + 4*lam2)**0.5)/2)); assert (Dg-1)*(Dg-2) == lam2
        Dl = J + Dg
        T = Fr((Dl-1)*(Dl-2))
    else:
        T = Fr(mspec)**2*lam2 + J*(J+1)
    m2 = T - J*(J+1)
    ctx.prec = prec
    Ta = arb(fmpq(T.numerator, T.denominator))
    H = xl.H_ang(J, b, Ta, Z, Y, Q, lam2)
    Hc = xl.H_ang(J, b, Ta, Zc, Yc, Q, lam2)
    w3 = (1-arb(fmpq(J*(J+1)))/Ta)**3
    cf = closed_form(J, b, m2, lam2)
    cfa = arb(fmpq(cf.numerator, cf.denominator))
    r = (H/w3-cfa)/cfa; rc = (Hc/w3-cfa)/cfa
    print(f"{J} {b} {mspec} | {(H/w3).str(12)} | {cfa.str(12)} | {r.str(3)} | {rc.str(3)}", flush=True)
