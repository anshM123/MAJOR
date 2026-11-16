"""Exact data and helpers for the F_{g3} heavy-action certificates (FLAT_GRAVITON_BOUND.md).

Units M = 1.  Heavy state: mass m >= 1, spin J, helicity channel h in {'++','+-'}.
u = 1/m^2 in (0,1].  Normalised heavy action (section 2 of the .md):

  FF_h(J,u) := m^4 F_{g3}[m,J,h]
             = int_0^1 Psi_h(p;u) D_h^J(1-2 u p^2) dp + Phi_h(J) u^3,

  Psi_h(p;u) = psi2(p) (2 - u p^2) + sigma_h u psi3(p),  sigma_{++}=+1, sigma_{+-}=-1,
  psi2 = (1-p)^3 p (65+155p+47p^2),   psi3 = 2 (1-p)^3 (81-1927p+4601p^2),
  D_{++}^J = Legendre P_J,  D_{+-}^J = dtilde^J_{44} = Jacobi P^{(0,8)}_{J-4},
  Phi:  literal functional  ++: 10 J(J+1) - 5,   +-: 10 J(J+1) - 205
        g5-free functional  ++: 10 J(J+1),       +-: 10 J(J+1) - 210.
All certificates use the PESSIMISTIC choice Phi_{++}=10J(J+1)-5, Phi_{+-}=10J(J+1)-210, which is
smaller than (or equal to) both variants in each channel, so positivity covers both functionals.
"""
from fractions import Fraction as Fr
from math import factorial, comb
import flint

PHI_CONST = {'++': 5, '+-': 210}          # pessimistic: Phi = 10*JJ - PHI_CONST
SIGMA = {'++': 1, '+-': -1}

# ---------------------------------------------------------------- exact polynomials
def pmul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r

def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def pscale(a, c):
    return [c * x for x in a]

ONE_MINUS_P_CUBED = [Fr(1), Fr(-3), Fr(3), Fr(-1)]
PSI2 = pmul(ONE_MINUS_P_CUBED, [Fr(0), Fr(65), Fr(155), Fr(47)])        # degree 6
PSI3 = pmul(ONE_MINUS_P_CUBED, [Fr(162), Fr(-3854), Fr(9202)])           # degree 5

def peval(a, x):
    r = 0
    for c in reversed(a):
        r = r * x + c
    return r

def pder(a, k=1):
    for _ in range(k):
        a = [i * a[i] for i in range(1, len(a))] or [Fr(0)]
    return a

def mom(a, j):
    """int_0^1 a(p) p^j dp, exact."""
    return sum(Fr(c) / (i + j + 1) for i, c in enumerate(a))

def mu2(j): return mom(PSI2, j)
def mu3(j): return mom(PSI3, j)

def e_h(h, i):
    """D^J(1-2y) = sum_k (-1)^k prod_{i<=k} (JJ - e_h(i)) y^k/(k!)^2 ,  JJ = J(J+1)."""
    return i * (i - 1) if h == '++' else (i + 3) * (i + 4)

def T0(n): return 2 * mu2(2 * n)
def T1(n, h): return -mu2(2 * n + 2) + SIGMA[h] * mu3(2 * n)

# ---------------------------------------------------------------- per-J exact polynomial in u
def poly_u(J, h):
    """Exact rational coefficients (list, index = power of u) of FF_h(J,u) (pessimistic Phi)."""
    JJ = J * (J + 1)
    kmax = J if h == '++' else J - 4
    co = {}
    prod = Fr(1)
    for k in range(0, kmax + 1):
        if k > 0:
            prod *= Fr(-(JJ - e_h(h, k)), k * k)
        co[k] = co.get(k, 0) + prod * T0(k)
        co[k + 1] = co.get(k + 1, 0) + prod * T1(k, h)
    co[3] = co.get(3, 0) + (10 * JJ - PHI_CONST[h])
    deg = max(co)
    return [Fr(co.get(k, 0)) for k in range(deg + 1)]

def to_fmpz_poly(c):
    from math import lcm
    L = 1
    for x in c:
        L = lcm(L, x.denominator)
    return flint.fmpz_poly([int(x * L) for x in c]), L

# ---------------------------------------------------------------- Clebsch-Gordan weights
# (1-y)^8 dtilde^J_44(1-2y) = sum_{j=-4}^{4} c_j(J) P_{J+j}(1-2y),  c_j = <J 4; 4 -4 | J+j 0>^2
def cg_coeffs(J):
    """c_{-4..4}(J) for J (Fraction / arb / float / complex); rational functions of J."""
    return [
        (J+1)*(J+2)*(J+3)*(J+4)/(16*(2*J-5)*(2*J-3)*(2*J-1)*(2*J+1)),
        (J+2)*(J+3)*(J+4)/(4*(2*J-3)*(2*J-1)*(2*J+1)),
        7*(J-3)*(J+2)*(J+3)*(J+4)/(4*(2*J-5)*(2*J-1)*(2*J+1)*(2*J+3)),
        7*(J-3)*(J+3)*(J+4)/(4*(2*J-3)*(2*J+1)*(2*J+3)),
        35*(J-3)*(J-2)*(J+3)*(J+4)/(8*(2*J-3)*(2*J-1)*(2*J+3)*(2*J+5)),
        7*(J-3)*(J-2)*(J+4)/(4*(2*J-1)*(2*J+1)*(2*J+5)),
        7*(J-3)*(J-2)*(J-1)*(J+4)/(4*(2*J-1)*(2*J+1)*(2*J+3)*(2*J+7)),
        (J-3)*(J-2)*(J-1)/(4*(2*J+1)*(2*J+3)*(2*J+5)),
        J*(J-3)*(J-2)*(J-1)/(16*(2*J+1)*(2*J+3)*(2*J+5)*(2*J+7)),
    ]
# factor structure of c_j as (numerator shifts a: (J+a)), (denominator shifts b: (2J+b)), constant
CG_FACTORS = [
    (Fr(1, 16), [1, 2, 3, 4], [-5, -3, -1, 1]),
    (Fr(1, 4), [2, 3, 4], [-3, -1, 1]),
    (Fr(7, 4), [-3, 2, 3, 4], [-5, -1, 1, 3]),
    (Fr(7, 4), [-3, 3, 4], [-3, 1, 3]),
    (Fr(35, 8), [-3, -2, 3, 4], [-3, -1, 3, 5]),
    (Fr(7, 4), [-3, -2, 4], [-1, 1, 5]),
    (Fr(7, 4), [-3, -2, -1, 4], [-1, 1, 3, 7]),
    (Fr(1, 4), [-3, -2, -1], [1, 3, 5]),
    (Fr(1, 16), [0, -3, -2, -1], [1, 3, 5, 7]),
]

def cg_abs_bound(Zmin):
    """Upper bound for sum_j |c_j(J)| over complex J with |J+1/2| >= Zmin (> 4).
    (J+a) = Z + (a-1/2), (2J+b) = 2Z + (b-1); each ratio bound is decreasing in |Z|."""
    tot = Fr(0)
    for c, num, den in CG_FACTORS:
        t = c
        for a in num:
            t *= (Zmin + abs(Fr(a) - Fr(1, 2)))
        for b in den:
            t /= (2 * Zmin - abs(Fr(b) - 1))
        tot += t
    return tot

# ---------------------------------------------------------------- arb helpers
def arb_of(x):
    if isinstance(x, Fr):
        return flint.arb(flint.fmpq(x.numerator, x.denominator))
    return flint.arb(x)

def arb_interval(lo, hi):
    """Ball containing the interval [lo, hi] (lo, hi Fractions or ints)."""
    lo = arb_of(Fr(lo)); hi = arb_of(Fr(hi))
    return lo.union(hi)

def upper(x):
    """Rigorous float-free upper bound (as arb with zero radius rounded up): returns arb.upper()."""
    return x.upper()

def set_prec(bits):
    flint.ctx.prec = bits
