"""
Free-field values of the canonical couplings of the higher-spin currents J_s (s = 4, 6), and check of the closed
forms (top components) at the unitarity bound:  |ahat^{J_s}_{s,2}|^2 = lam_pm^2 |K^{+-}_s|^2/C_T = 1,
|ahat^{J_s}_{-s,2}|^2 = lam_pp^2 |K^{++}_{-s}|^2/C_T = 1.   (lam for unit-normalised J_s and canonical T.)
"""
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
import wick_free as wf, kernels_lib as kl, closed_forms as cf, free_checks as fc
from canonical import canonical
mp.mp.dps = 30
S = kl.S
fourpi = 4*mp.pi
CT = mp.mpf(3)/(32*mp.pi**2)
nuS = {4: 645120, 6: 30656102400}; nuF = {4: 1290240, 6: 61312204800}
for s_ in (4, 6):
    can = canonical(s_)
    sig = sp.Rational(2*s_ + 1, 2)
    monos = can['pp'][1]
    vpp = [sp.sympify(v).subs(S, sig) for v in can['pp'][2]]; vpm = [sp.sympify(v).subs(S, sig) for v in can['pm'][2]]
    for theory in ('scalar', 'fermion'):
        if theory == 'scalar':
            _, f = fc.fit('scalar', s_, s_ + 1, wf.c_scalar_J(s_))
            # lam (unit-normalised O, canonical T) = l (4 pi)^-3 / sqrt(nu_O (4 pi)^-2)
            scale = fourpi**(-3)/mp.sqrt(nuS[s_]*fourpi**(-2)); ct = CT
        else:
            _, f = fc.fit('fermion', s_, s_ + 1, wf.c_fermion_J(s_), 'z')
            # Dirac: lam = kappa^2 c_F^3 l/(c_F sqrt(nu_O)), (kappa c_F)^2 = C_T^Dirac/12; C_T^Dirac = 2 C_T
            scale = (2*CT/12)/mp.sqrt(nuF[s_]); ct = 2*CT
        m2, lam = f['even']
        assert m2 == monos
        A = sp.Matrix([vpp, vpm]).T
        sol, params = A.gauss_jordan_solve(sp.Matrix(lam))
        assert len(params) == 0
        xpp, xpm = [mp.mpf(sp.Rational(v).p)/sp.Rational(v).q*scale for v in sol]
        Delta = mp.mpf(s_ + 1)
        a_top = xpm**2*cf.K2(cf.Mtop_pm(Delta, s_), Delta, s_, s_)/ct
        a_bot = xpp**2*cf.K2(cf.Mtop_pp(Delta, s_), Delta, s_, -s_)/ct
        print(f"{theory:7s} J_{s_}: lam_pp/sqrt(C_T) = {mp.nstr(xpp/mp.sqrt(ct), 12)}, lam_pm/sqrt(C_T) = {mp.nstr(xpm/mp.sqrt(ct), 12)};"
              f"  |ahat_({s_},2)|^2 = {mp.nstr(a_top, 15)}, |ahat_(-{s_},2)|^2 = {mp.nstr(a_bot, 15)}")
