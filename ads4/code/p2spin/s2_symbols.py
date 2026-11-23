"""
p2spin/s2_symbols.py (P2_SPINNING.md §2): heavy symbols of the conserved <T T O_{Delta,J}> structures in CHMRS's symmetric frame.

Frame (CHMRS 2106.10274 app. A, eq. otherframe): x1 = -y, x2 = -x/(-x^2), x3 = y, x4 = x/(-x^2), x5 = e; rescaled x, y -> eps x, eps y,
eps = 1/m, with  Delta + J - 1 = m q,  Delta - J - 2 = m/q,  q = e^{beta/2}   (so J + 1/2 = m sinh(beta/2), m^2 = (Delta-J-2)(Delta+J-1)).
Light-cone components (u, v, y) as tto_conserved; e = (1, 1, 0) (unit time), n = (2, 0, 0) (null, = (1, e_z) in Cartesian (t, z, y)).
Polarisations at the far points are referred by the inversion tensor: z2 = I(x) z2', z4 = I(x) z4', I(x) = 1 - 2 x x^T g / x^2
(eps-independent), so the leading symbols are finite constant tensors.

Left symbol:  R_L = sum_mu c_mu(Delta, J) N_mu(x1, z1, x2, z2, e, n) / X12^2      (N_mu = CPPR monomial / V5^J; p2_core)
= <T T O>_{spinning} / <phi phi O>_{scalar seed, Delta_phi = 3}  (the seed carries V5^J and the X-powers; CHMRS's T(x, y)).
Right symbol: R_R = sum_mu c_mu(3 - Delta, -1 - J) N_mu(x3, z3, x4, z4, e, n) / X34^2   (Weyl-reflected family, W1a).

Output per basis structure: Laurent coefficients R^(0), R^(1), R^(2) (exact, polynomial in x, y and the polarisations).
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import p2_core as C
import p2_series as PS
from p2_series import TS, EPS

Dl, Jl = sp.symbols('Delta J')
q = sp.Symbol('q', positive=True)
xs = sp.symbols('xu xv xy'); ys = sp.symbols('yu yv yy')
E = (sp.Integer(1), sp.Integer(1), sp.Integer(0))
NN = (sp.Integer(2), sp.Integer(0), sp.Integer(0))


def CAN():
    return {k: [[sp.sympify(eval(x, sp.__dict__)) for x in v] for v in B]
            for k, B in pickle.load(open(os.path.join(HERE, 'cache', 's1b_canonical.pkl'), 'rb')).items()}


def inv_tensor(x):
    """I(x)^mu_nu = delta - 2 x^mu x_nu / x^2 (acting on upper-index vectors)."""
    xl = [sum(C.GL[m][n]*x[n] for n in range(3)) for m in range(3)]
    x2 = C.dot(x, x)
    return [[(1 if m == n else 0) - 2*x[m]*xl[n]/x2 for n in range(3)] for m in range(3)]


def apply(Mx, z):
    return tuple(sum(Mx[m][n]*z[n] for n in range(3)) for m in range(3))


def frame_points(side, xv, yv):
    """(xa, xb) = (x1, x2) for 'L', (x3, x4) for 'R', with x, y -> eps x, eps y."""
    ex = [EPS*c for c in xv]; ey = [EPS*c for c in yv]
    exx = C.dot(ex, ex)
    if side == 'L':
        return tuple(-c for c in ey), tuple(c/exx for c in ex)          # x1 = -y, x2 = -x/(-x^2) = x/x^2
    return tuple(ey), tuple(-c/exx for c in ex)                         # x3 = y, x4 = x/(-x^2)


def block_series(side, xv, yv, za, zbp, n=NN, e=E):
    xa, xb = frame_points(side, xv, yv)
    zb = apply(inv_tensor(xv), zbp)
    B = C.blocks(xa, za, xb, zb, e, n)
    out = {k: TS.from_expr(B[k]) for k in ('V1', 'V2', 'V5', 'H12', 'H15', 'H25', 'X12', 'e12', 'e15', 'e25')}
    # parity-odd common factor (|X12| |X15| |X25|)^(-1/2) relative to the even seed (p2_core docstring); |.| as in CHMRS
    X12X25 = TS.from_expr(B['X12']*B['X25']); mX15 = TS.from_expr(-B['X15'])
    xx = C.dot(xv, xv)                                   # x^2 < 0 (timelike): leading of X12 X25 is 1/(eps^4 x^4)
    out['oddfac'] = (X12X25.sqrt(-1/xx)*mX15.sqrt(sp.Integer(1))).inv()
    return out


def mono_series(Bs, sector, t):
    if sector == 'even':
        a, b, c = t
        return Bs['V1']**(2 - a - b)*Bs['V2']**(2 - a - c)*Bs['V5']**(-b - c)*Bs['H12']**a*Bs['H15']**b*Bs['H25']**c*Bs['X12']**(-2)
    typ, a, b, c = t
    pref = {'12': (Bs['e12'], 1, 1, 0), '15': (Bs['e15'], 1, 2, 1), '25': (Bs['e25'], 2, 1, 1)}[typ]
    eb, p1, p2, s5 = pref
    return Bs['oddfac']*eb*Bs['V1']**(p1 - a - b)*Bs['V2']**(p2 - a - c)*Bs['V5']**(-s5 - b - c)*Bs['H12']**a*Bs['H15']**b*Bs['H25']**c*Bs['X12']**(-2)


def coeff_series(cexpr, side, qv):
    """c_mu(Delta, J) with Delta, J from (m = 1/eps, q); right side: Weyl-reflected (3 - Delta, -1 - J)."""
    m = 1/EPS
    D = (m*(qv + 1/qv) + 3)/2; J = (m*(qv - 1/qv) - 1)/2
    if side == 'R':
        D, J = 3 - D, -1 - J
    return TS.from_expr(sp.cancel(cexpr.subs({Dl: D, Jl: J}, simultaneous=True)))


def symbol(side, sector, coeffs, xv, yv, za, zbp, qv, monos=None):
    monos = monos or (C.EVEN if sector == 'even' else C.ODD)
    Bs = block_series(side, xv, yv, za, zbp)
    tot = None
    for cf, t in zip(coeffs, monos):
        if cf == 0:
            continue
        term = coeff_series(cf, side, qv)*mono_series(Bs, sector, t)
        tot = term if tot is None else tot + term
    return tot
