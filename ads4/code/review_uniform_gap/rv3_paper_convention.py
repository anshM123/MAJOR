"""rv3_paper_convention.py -- check (1c): the paper's shifted density E(mu) = 2 Re[chi(mu-i) R_-(mu-i)] of the PUBLISHED Delta_g=1000
continuum weight (k2_continuum_D_1000_rep.json, chi_re part) against psi2(mu/Lambda) (NF02 with 1/k0), i.e. the same equation and
normalisation as the note's R_-(x-i)chi(x-i)+R_+(x+i)chi(x+i) = f(x).  Also compares the sum-of-two-terms form."""
import json, os
import mpmath as mp
from rv_common import *
mp.mp.dps = 30
d = json.load(open(os.path.join(DATA, 'k2', 'k2_continuum_D_1000_rep.json')))
L = mp.sqrt(997002)
cs = [(mp.mpf(c), mp.mpf(s)) for c, s in d['centres_sigma']]
xr = [mp.mpf(v) for v in d['x_re']]; yr = [mp.mpf(v) for v in d['y_re']]

def chi_re(nu):
    Q = nu**2/(nu**2+mp.mpf(1)/4)**2
    S = 0
    for (c, s), x in zip(cs, xr):
        if abs(mp.re(nu)-c) > 40*s and abs(mp.re(nu)+c) > 40*s:
            continue
        g = mp.exp(-(nu-c)**2/(2*s*s))
        if c > 0:
            g += mp.exp(-(nu+c)**2/(2*s*s))
        S += x*g
    l = (nu**2+1)*mp.log(nu**2+1)
    return Q*S + (yr[0]+yr[1]*nu**2+yr[2]*nu**4)*l*mp.exp(-nu**2/2)

mx = 0
psimax = max(psi_mp(PSI2, mp.mpf(i)/200) for i in range(201))
print("mu | psi2(mu/L) | 2Re[R_-chi](mu-i) | R_-chi(mu-i)+R_+chi(mu+i) | diff/max psi")
for mu in [0.05, 0.5, 2.0, 7.3, 50.0, 300.0, 700.0, 950.0, 998.0]:
    mu = mp.mpf(mu)
    a = Rm(mp.mpc(mu, -1))*chi_re(mp.mpc(mu, -1))
    b = Rp(mp.mpc(mu, 1))*chi_re(mp.mpc(mu, 1))
    E = 2*mp.re(a)
    p = psi_mp(PSI2, mu/L)
    mx = max(mx, abs(E-p)/psimax)
    print(mp.nstr(mu, 6), mp.nstr(p, 10), mp.nstr(E, 10), mp.nstr(a+b, 10), mp.nstr(abs(E-p)/psimax, 3))
print("max |E - psi2(mu/L)| / max psi2 =", mp.nstr(mx, 3), "  (max psi2 =", mp.nstr(psimax, 6), ")")
