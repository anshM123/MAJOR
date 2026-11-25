"""ANGULAR_MATCH item 1: convergence bookkeeping for weighted light-ray products (KKSZ 1904.05905 §4.2-4.6), and the
impact parameter carried by the abstract defect examples of ANGULAR_ROUTE §10, §12, §13.

KKSZ (4.86): L^n = int dv v^{n+1} T_vv, so E = L^{-1} and the first moment E^(1) = L^0 (4.88).  In the alpha-parametrisation
(4.17) the weight v^{n+1} is |alpha|^{n+1} at the ends alpha -> +-oo, where the two contours meet.  Every convergence condition
that comes from an alpha -> oo end therefore shifts:
  Regge corner (both ends):        J1 + J2 - (n1+1) - (n2+1) > 1 + J0             (= KKSZ (4.87))
  OPE corner (both ends):          J1 + J2 - (n1+1) - (n2+1) > 2 - Delta0'
  light-cone edge, alpha_i -> oo:  J1 + J2 +- (Delta1 - Delta2) - 2(n_i+1) > 2 - tau  (KKSZ (4.27)/(4.35) shifted by our count)
with tau = tau0' (first sheet, rigorous; second sheet under KKSZ's asymptotic light-cone assumption (4.66)),
     tau = 0     (second sheet, Rindler positivity only: KKSZ §4.2.7 'conditions with tau0' = 0'),
     tau = tau_* (double-commutator integral, KKSZ §4.3; tau_* >= d-2 = 1).
Output: cache/am1_kksz.txt
"""
import os, math
out = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'am1_kksz.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

d = 3; J1 = J2 = 2; D1 = D2 = 3
tau0p_min = 0.5        # smallest nonzero twist allowed in T x T in d = 3 (scalar unitarity bound)
Delta0p_min = 0.5
tau_star_min = d - 2
P('Margins (left side minus right side; > 0 means the KKSZ sufficient condition holds), T detectors in d = 3')
P(f'{"pair":10s} {"Regge J0<":>10s} {"m(J0=1)":>8s} {"m(J0=2)":>8s} {"OPE":>6s} {"LC tau0p":>9s} {"LC tau=0":>9s} {"LC tau*":>8s}')
for name, n1, n2 in [('E E', -1, -1), ('E E(1)', -1, 0), ('E(1) E(1)', 0, 0)]:
    w = (n1+1) + (n2+1)
    regge_bound = J1 + J2 - w - 1            # need J0 < this
    m1, m2 = regge_bound - 1, regge_bound - 2
    ope = J1 + J2 - w - (2 - Delta0p_min)
    lc = lambda tau: min(J1 + J2 + s*(D1-D2) - 2*(ni+1) - (2 - tau) for s, ni in ((+1, n1), (-1, n2)))
    P(f'{name:10s} {regge_bound:>10d} {m1:>8d} {m2:>8d} {ope:>6.1f} {lc(tau0p_min):>9.1f} {lc(0):>9.1f} {lc(tau_star_min):>8.1f}')
P('check against KKSZ (4.87) J1+J2 > J0+3+n+m: E E -> J0<3, E E(1) -> J0<2, E(1) E(1) -> J0<1 (reproduced above).')
P('Reading: the mixed product passes Regge with margin 1 at J0 <= 1 and is exactly marginal at J0^planar = 2;')
P('its light-cone edge has margin tau0\' >= 1/2 (or tau_* >= 1) but is marginal (0) with Rindler positivity alone.')

P('\nImpact parameter e^beta = (Delta+J)/(Delta-J) of the states carrying the defect in the abstract examples:')
for n in [4, 6, 8, 10]:
    L = 2**n
    Ja, Da = L**3, L**3 + L**2 + 2          # §10
    Jb, Db = L**2, L**2 + L + 2             # §12
    l = 2*n; Jc = 2**l; Dc = Jc + math.ceil(Jc/l) + 2   # §13 (ell = 2n)
    eb = lambda J, D: (D+J)/(D-J)
    P(f'  L=2^{n:<2d}: §10 e^beta={eb(Ja,Da):10.2f} (~2L={2*L});  §12 e^beta={eb(Jb,Db):10.2f} (~2L);  §13 (ell={l}) e^beta={eb(Jc,Dc):8.2f} (~2 ell={2*l})')
P('All three defect constructions put the weight at e^beta -> oo: the defect is a large-impact-parameter (Regge-corner) effect.')
