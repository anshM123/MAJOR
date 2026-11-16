"""
Goal 2 (sect. 6 of SPINNING_MFT.md): spin dependence of the O(1/m^2) coefficient of the helicity-preserving component
after MFT normalisation.  Input: the measured coefficients c_2 of SPINNING_E2E.md sect. 6.2 (k = 6; fits (B) there;
copied verbatim below, together with the scalar-skeleton coefficients of the same table).  Nothing is fitted here.

MFT normalisation multiplies the e2e ratio R by Q_E (mft_heavy.py), whose 1/m^2 coefficient is
   q_2 = kappa_l [1/sinh^2(b/2) + 1/cosh^2(b/2)] - (8/3) B_3((Dphi-1)/2) cosh b,   kappa_1 = 1/4, kappa_2 = 5/2
(derived), so c_2^MFT = c_2 + q_2.  The MFT-normalised scalar skeleton (D1, Dphi) has c_2^MFT = c_2^skel - (8/3)B_3((D1+Dphi-2)/2) cosh b.
Tests:  (a) c_2^MFT(l=2) - c_2^MFT(l=1) = -6 P'/P ?          (flat-space helicity term -2 l^2 P'/P, difference of l = 2 and 1)
        (b) c_2^MFT(l) + 2 l^2 P'/P - c_2^MFT,skel(2, Dphi) = ?  (residual, both l)
"""
import mpmath as mp
mp.mp.dps = 20
Dphi = mp.mpf(7)/4
# beta, nu : (c2 pres l=1, c2 pres l=2, skeleton (2,Dphi), skeleton (3,Dphi))      [SPINNING_E2E.md sect. 6.2, k = 6]
T = {(0.4, 0): (-96.94, -153.85, -90.95, -72.28), (0.4, 1): (-89.53, -143.50, -84.56, -66.05), (0.4, 3): (-25.47, -50.01, -30.32, -13.47),
     (1.0, 0): (-127.43, -136.87, -126.66, -100.23), (1.0, 1): (-102.01, -108.60, -102.38, -77.27), (1.0, 3): (-205.37, -235.23, -197.89, -162.09),
     (2.0, 0): (-295, -297, -293.6, -230.4), (2.0, 1): (-60, -59, -63.0, -14.9), (2.0, 3): (-730, -741, -723.3, -618.1)}
B3 = lambda s: s*(s - mp.mpf(1)/2)*(s - 1)
print("| beta | nu | -6P'/P | c2MFT(l=2)-c2MFT(l=1) | difference | residual l=1 | residual l=2 |")
print("|---|---|---|---|---|---|---|")
for (b, nu), (c1, c2, s2, s3) in sorted(T.items()):
    eta = mp.cosh(b)
    P = lambda e: mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, (1 - e)/2))
    PpP = mp.diff(P, eta)/P(eta)
    pol = 1/mp.sinh(b/2)**2 + 1/mp.cosh(b/2)**2
    g = -mp.mpf(8)/3*B3((Dphi - 1)/2)*eta
    m1 = c1 + mp.mpf(1)/4*pol + g; m2 = c2 + mp.mpf(5)/2*pol + g
    sk2 = s2 - mp.mpf(8)/3*B3(Dphi/2)*eta
    print("| %.1f | %d | %s | %s | %s | %s | %s |" % (b, nu, mp.nstr(-6*PpP, 5), mp.nstr(m2 - m1, 5), mp.nstr(m2 - m1 + 6*PpP, 2),
          mp.nstr(m1 + 2*PpP - sk2, 3), mp.nstr(m2 + 8*PpP - sk2, 3)))

# per-beta fits, as in SPINNING_E2E.md sect. 6.3 ("direct test"), but with the l = 1 data (MFT-normalised) as reference instead of the
# scalar skeleton (3, Dphi):  c2MFT(l=2) - [c2MFT(l=1) + 2 P'/P] = c_beta + d_beta P'/P ;  prediction of the flat-space term: d_beta = -8
print()
print("| beta | d_beta (fit c + d P'/P over nu = 0,1,3) | c_beta | max residual with d fixed to -8 (c refitted) |")
print("|---|---|---|---|")
for b in (0.4, 1.0, 2.0):
    xs, ys = [], []
    for nu in (0, 1, 3):
        c1, c2, s2, s3 = T[(b, nu)]
        eta = mp.cosh(b)
        P = lambda e: mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, (1 - e)/2))
        PpP = mp.diff(P, eta)/P(eta)
        pol = 1/mp.sinh(b/2)**2 + 1/mp.cosh(b/2)**2
        ys.append((c2 + mp.mpf(5)/2*pol) - (c1 + mp.mpf(1)/4*pol + 2*PpP)); xs.append(PpP)
    n = 3; sx = sum(xs); sy = sum(ys); sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    d = (n*sxy - sx*sy)/(n*sxx - sx*sx); c = (sy - d*sx)/n
    c8 = (sy + 8*sx)/n
    res = max(abs(y - (c8 - 8*x)) for x, y in zip(xs, ys))
    print("| %.1f | %s | %s | %s |" % (b, mp.nstr(d, 4), mp.nstr(c, 3), mp.nstr(res, 2)))
