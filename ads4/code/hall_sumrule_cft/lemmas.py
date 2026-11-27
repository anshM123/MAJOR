"""
Algebraic lemmas used in README.md (all exact / sympy unless stated).

 L1  Energy matrix of CDP 1707.03007 in the helicity basis: diagonal = E/4pi (t4, alpha_T independent);
     helicity-flip element = -(E/4pi)(t4 +- i alpha_T)/4.  So the collider fixes the relative normalisation of the
     parity-even homogeneous structure (coefficient t4/4) and the parity-odd one (coefficient alpha_T/4, multiplying
     the helicity-rotated structure R h, R = +-i on helicity +-), and the disc t4^2 + alpha_T^2 <= 16.
 L2  Rank/parity lemma: an SO(3)-invariant tensor of EVEN rank is reflection-even (contains no epsilon); an invariant
     of ODD rank is reflection-odd.  Checked numerically for ranks 4..7 by computing the invariant subspace.
     Consequence: a T-dependent local (contact) term c * p^{k} <O> in <T T>_beta (rank 4 + spin(O) + k) can be
     parity-odd only for odd k, and dimension matching needs k = 3 - Delta_O; with unitarity this leaves only
     scalars with Delta_O = 2 (k = 1) [and the identity, k = 3: the T-independent gravitational CS term].
     In particular NO parity-odd contact term proportional to <T>_beta exists (k = 0).
 L3  Fourier transforms (Euclidean time, zero spatial momentum): an operator of dimension Delta contributes
     |tau|^{Delta-4} (even channel) or sgn(tau)|tau|^{Delta-4} (odd channel) to the plane-integrated correlator.
       Int e^{i w tau} sgn(tau) |tau|^{-a} dtau = 2 i Gamma(1-a) cos(pi a/2) |w|^{a-1} sgn(w)
       Int e^{i w tau}          |tau|^{-a} dtau = 2   Gamma(1-a) sin(pi a/2) |w|^{a-1}
     a = 1 (T, Delta = 3): odd -> i pi sgn(w) (finite constant);  even -> pole (log|w|), so the T term in g has NO
     1/|tau| tail and its UV constant is a delta(tau) (scheme-sensitive), while the Hall constant is a genuine
     separated-point 1/tau tail.  a = 2 (Delta = 2): odd -> pole (w log|w| unless the residue vanishes; at N = oo
     CS-fermion it is purely analytic w, see cs_fermion_exact.py), even -> finite |w|.
 L4  Assembly: S_+- := S_E -+ i sigma S_H = 3/8 + (t4 -+ i sigma alpha_T)/32 and the CS-fermion / bulk numbers.
"""
import itertools

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- L1
E, t4, aT = sp.symbols('E t4 alpha_T', real=True)
M = E/(4*sp.pi)*sp.Matrix([[1 - t4/4, aT/4], [aT/4, 1 + t4/4]])
vp = sp.Matrix([1, sp.I])/sp.sqrt(2)
vm = sp.Matrix([1, -sp.I])/sp.sqrt(2)
Mh = sp.Matrix(2, 2, lambda i, j: sp.simplify(([vp, vm][i].H*M*[vp, vm][j])[0]))
print('L1  energy matrix in helicity basis (v+, v-):')
sp.pprint(sp.simplify(Mh*4*sp.pi/E))
assert sp.simplify(Mh[0, 1] + E/(4*sp.pi)*(t4 + sp.I*aT)/4) == 0
print('    eigenvalues:', [sp.simplify(ev) for ev in M.eigenvals()], ' -> positivity <=> t4^2 + alpha_T^2 <= 16')

# ---------------------------------------------------------------- L2
Lgen = [np.zeros((3, 3)) for _ in range(3)]
for i, j, k in itertools.permutations(range(3)):
    Lgen[i][j, k] = np.linalg.det(np.eye(3)[[i, j, k]])       # (L_i)_{jk} = eps_{ijk}
refl = np.diag([1.0, 1.0, -1.0])


def rep(Mat, n):
    """generator Mat acting on the n-fold tensor product (sum over slots)."""
    out = np.zeros((3**n, 3**n))
    for s in range(n):
        op = np.array([[1.0]])
        for t in range(n):
            op = np.kron(op, Mat if t == s else np.eye(3))
        out += op
    return out


print('\nL2  SO(3)-invariant tensors of rank n and their reflection parity:')
for n in (4, 5, 6, 7):
    A = np.vstack([rep(L, n) for L in Lgen])
    _, s, Vt = np.linalg.svd(A)
    null = Vt[np.sum(s > 1e-9):].T                              # invariant subspace
    R = np.array([[1.0]])
    for _ in range(n):
        R = np.kron(R, refl)
    par = np.linalg.eigvalsh(null.T @ R @ null) if null.shape[1] else np.array([])
    print('    rank %d: %3d invariants, reflection eigenvalues in {%s}' % (n, null.shape[1],
          ', '.join(sorted(set('%+d' % round(x) for x in par)))))
print('    => even rank: only epsilon-free (parity-even) invariants; odd rank: only epsilon-containing ones.')
print('    contact terms c p^k <O_{Delta,l}> in <T_A T_B>: rank = 4 + l + k, k = 3 - Delta;  parity-odd needs k odd:')
for l, Delta in [(2, 3), (0, 3), (0, 2), (0, 1), (4, 5)]:
    k = 3 - Delta
    ok = k >= 0 and (4 + l + k) % 2 == 1
    print('      O with spin %d, Delta = %d: k = %+d -> parity-odd local term %s' % (l, Delta, k, 'ALLOWED' if ok else 'impossible'))

# ---------------------------------------------------------------- L3
a, w, tau = sp.symbols('a w tau', positive=True)
odd = 2*sp.I*sp.gamma(1 - a)*sp.cos(sp.pi*a/2)
even = 2*sp.gamma(1 - a)*sp.sin(sp.pi*a/2)
chk_odd = sp.integrate(sp.sin(w*tau)*tau**(-sp.Rational(1, 3)), (tau, 0, sp.oo))
rhs = (odd/(2*sp.I)).subs(a, sp.Rational(1, 3))*w**(-sp.Rational(2, 3))
print('\nL3  check of the Mellin formula at a = 1/3:  Int_0^oo sin(w t) t^-1/3 dt =', sp.simplify(chk_odd),
      ' vs', sp.simplify(rhs), ' (difference at w = 1: %.1e)' % abs(sp.N((chk_odd - rhs).subs(w, 1), 30)))
print('    odd  channel, a -> 1 (T):        ', sp.limit(odd, a, 1), '* sgn(w)       [finite: Hall UV constant]')
print('    even channel, a -> 1 (T):         residue', sp.limit((a - 1)*even, a, 1), '-> log|w| unless the 1/|tau| tail vanishes')
print('    odd  channel, a -> 2 (Delta = 2): residue', sp.limit((a - 2)*odd, a, 2), '-> w log|w| unless position-space residue vanishes')
print('    even channel, a -> 2 (Delta = 2):', sp.limit(even, a, 2), '* |w|   [finite: CS-fermion a1 |q3|]')
print('    PV Int e^{i w tau}/tau dtau =', sp.integrate(2*sp.I*sp.sin(w*tau)/tau, (tau, 0, sp.oo)), '(w > 0)')

# ---------------------------------------------------------------- L4
print('\nL4  assembly')
lam = sp.symbols('lambda', positive=True)
SE = sp.Rational(3, 8) + t4/32
SH = aT/32
print('    S_E - 3/8 + i S_H = (t4 + i alpha_T)/32;  |.| = tau/32 <= 1/8 on the collider disc')
cs = {t4: -4*sp.cos(sp.pi*lam), aT: 4*sp.sin(sp.pi*lam)}
print('    CS-fermion (t4, alpha_T) = (-4 cos pi lam, 4 sin pi lam):  S_E - 3/8 =', sp.simplify(SE.subs(cs) - sp.Rational(3, 8)),
      ', S_H =', SH.subs(cs))
print('    cs_fermion_exact.py: Hall constant/eps = sin(pi lam)/8 exactly; shear constant/eps = -1/8 + cos(pi lam)/8 exactly')
print('    bulk (bulk_hall_exact.py): S_H = -45 lam_Q = 7.5 alpha3~;  |alpha_T| = 240|alpha3~|  ->  |S_H| = |alpha_T|/32')
