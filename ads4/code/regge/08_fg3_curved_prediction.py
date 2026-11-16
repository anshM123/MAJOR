r"""Check 8 (conditional prediction, NOT a derivation): the functional F_g3 of arXiv:2201.06602 (3.3a) (reading C of
06_fg3_reading.py) with every flat partial wave dtilde^J(1 - 2p^2/m^2) replaced by the scalar harmonic kernel
P_{-1/2+i nu}(cosh beta), nu = M p (pair weights (0,0) for all three sum rules B_2^(1), B_3^(1), B_4^(1); see
REGGE_FROM_BLOCKS.md sections 5-7).  m is in units of the gap, M = Delta_gap in AdS units.  Valid as a model of the
regime m >> 1 only (the Regge formula needs nu << m M); corrections of relative order 1/m^2 are not known.
    F^{++/+-} = \int_0^1 dp [ psi_2(p) (2m^2-p^2) +/- psi_3(p) ] P(M p; beta)/m^6 + (5/m^10) [ 2 m^2 M^2 (-d/dnu^2 P)|_0 - P(0; beta) ].
Printed: m^4 M^2 F  (so that the leading term is 2 (65 Q_1 + 40 |Q_2|/M))."""
import mpmath as mp
mp.mp.dps = 15
P = lambda nu, be: mp.re(mp.hyp2f1(0.5+1j*nu, 0.5-1j*nu, 1, -mp.sinh(be/2)**2))
psi2 = lambda p: (1-p)**3*p*(65+155*p+47*p*p)
psi3 = lambda p: 2*(1-p)**3*(81-1927*p+4601*p*p)
def F(M, m, be, sgn):
    n = max(4, int(M*be/3) + 4)
    I = mp.quad(lambda p: (psi2(p)*(2*m*m-p*p) + sgn*psi3(p))*P(M*p, be), mp.linspace(0, 1, n))/m**6
    dP = -mp.diff(lambda x: P(x, be), 0, 2)/2
    return I + 5*(2*m*m*M*M*dP - P(0, be))/m**10
for M in (20, 100):
    print(f"M = {M}:   beta      m     m^4 M^2 F^(++)    m^4 M^2 F^(+-)")
    for be in map(mp.mpf, ('0.3', '1.0', '2.0')):
        for m in map(mp.mpf, (3, 10, 30, 100, 1000)):
            a, b = F(M, m, be, +1), F(M, m, be, -1)
            print(f"           {mp.nstr(be,2):>4} {mp.nstr(m,4):>6} {mp.nstr(a*m**4*M*M, 6):>16} {mp.nstr(b*m**4*M*M, 6):>16}")
