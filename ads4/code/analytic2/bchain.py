"""
bchain.py: Lemma B (beta-chain by chained similarities) -- certified bounds of |S^(beta0 - j)|, j = 1..J, from bounds of
|S^(beta0)| and |S^(beta0 + 1)|, uniformly over a parameter box.  ANALYTIC_SLIVER2 §2.

Exact recurrence (TAIL_CERT E3, d-normalised; S^(beta) = d^beta S(beta), g := 1/(d^2 D), D = N(N+s)):
    S^(beta - 1) = A(beta) S^(beta) + B(beta) S^(beta + 1),
    A = [1 + beta(s - 2beta - 1) g d] / [1 - beta(beta - s) g d^2],   B = beta(beta + 1) g / [1 - beta(beta - s) g d^2].
With x_j := S^(beta0 - j) (x_{-1} = S^(beta0+1)), V_j = (x_j, x_{j-1})^T:  V_{j+1} = M_j V_j,  M_j = [[A_j, B_j], [1, 0]],
A_j = A(beta0 - j).  For ANY invertible point matrices T_j (chosen numerically: eigenvectors, or a scaled Schur basis near a
double root) and c_j := T_j^{-1} V_j:
    ||c_{j+1}||_inf <= ||T_{j+1}^{-1} M_j T_j||_inf ||c_j||_inf,    |x_j| <= ||row_1(T_j)||_1 ||c_j||_inf,
where M_j is the Arb ball matrix over the whole parameter box (tau, d, g, s) and T^{-1} is the Arb inverse.  This is rigorous for
every T_j; only its sharpness depends on the choice.  (Lemma B.)

State relations (exact): t := a/N solves t^2 - (2+e)t + d(1+e) = 0 with e = s/N, and g = (t/d)^2/(a^2(1+e)).
"""
import math, sys
import numpy as np
import scipy.linalg as sla
from flint import arb, acb, ctx, acb_mat
ctx.prec = 128


def up(x):
    return float(arb(x).upper())


def lo_(x):
    return float(arb(x).lower())


def arb_iv(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return arb((lo + hi)/2, ((hi - lo)/2).upper()*(1 + 1e-15) + 1e-300)


def state_box(a_lo, a_hi, d_lo, d_hi, s_list=(1, 9)):
    """Arb intervals for d and g = 1/(d^2 D) over a in [a_lo, a_hi] (a_hi may be math.inf), d in [d_lo, d_hi], s in s_list."""
    d = arb_iv(d_lo, d_hi)
    smax = max(s_list)
    # t = a/N <= 0.51 d (for d <= 0.1), e = s/N = s t / a
    e = arb_iv(0, smax*0.51*d_hi/a_lo)
    def t_of(dd):
        return ((2 + e) - ((2 + e)**2 - 4*dd*(1 + e)).sqrt())/2
    r_lo = t_of(arb(d_lo))/arb(d_lo); r_hi = t_of(arb(d_hi))/arb(d_hi)     # t/d increasing in d
    tq_lo = min(lo_(r_lo), lo_(r_hi)); tq_hi = max(up(r_lo), up(r_hi))
    g_hi = up(arb(tq_hi)**2/arb(a_lo)**2)
    g_lo = 0.0 if a_hi == math.inf else lo_(arb(tq_lo)**2/(arb(a_hi)**2*(1 + arb(smax*0.51*d_hi/a_lo))))
    g = arb_iv(g_lo, g_hi)
    s = arb_iv(min(s_list), max(s_list))
    return d, g, s


def coef(beta, d, g, s):
    den = 1 - beta*(beta - s)*g*d*d
    A = (1 + beta*(s - 2*beta - 1)*g*d)/den
    B = beta*(beta + 1)*g/den
    return A, B


def _mid(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))


def _frame(Mmid, eta=0.15, cond_max=60.0):
    """point matrix T (numpy complex) for the frozen matrix Mmid: normalised eigenvectors if well conditioned, else scaled Schur."""
    w, V = np.linalg.eig(Mmid)
    V = V/np.linalg.norm(V, axis=0)
    if np.linalg.cond(V) <= cond_max:
        # order: larger modulus first
        idx = np.argsort(-np.abs(w))
        return V[:, idx]
    Tm, Z = sla.schur(Mmid, output='complex')
    if abs(Tm[0, 0]) < abs(Tm[1, 1]):
        # reorder so that the larger eigenvalue is first (Schur with a swapped target): use the other eigenvector as first column
        v = V[:, np.argmax(np.abs(w))]
        v = v/np.linalg.norm(v)
        q = np.array([-np.conj(v[1]), np.conj(v[0])])
        Z = np.column_stack([v, q])
    return Z @ np.diag([1.0, eta])


def _to_acb(T):
    return acb_mat([[acb(complex(T[i, j])) for j in range(2)] for i in range(2)])


def _inf_norm(Mm):
    return max(up(abs(Mm[i, 0]) + abs(Mm[i, 1])) for i in range(2))


def chain_bounds(sig0, t_lo, t_hi, M0, M1, J, d_lo, d_hi, a_lo, a_hi, s_list=(1, 9), eta=0.15):
    """bounds b[j] >= |S^(sig0 - j + i tau)| for j = -1, 0, ..., J (returned as list indexed j+1), uniformly over
    tau in [t_lo, t_hi] (both signs of tau are covered by conjugation symmetry only if the caller passes |tau|), the state box."""
    d, g, s = state_box(a_lo, a_hi, d_lo, d_hi, s_list)
    tau = arb_iv(t_lo, t_hi)
    Ms = []; Mm = []
    for j in range(J):
        beta = acb(arb(sig0 - j), tau)
        A, B = coef(beta, d, g, s)
        Ms.append(acb_mat([[A, B], [acb(1), acb(0)]]))
        Mm.append(np.array([[_mid(A), _mid(B)], [1, 0]], dtype=complex))
    Ts = [_frame(Mm[0], eta)] + [_frame(Mm[j], eta) for j in range(J)]
    # T_j for j = 0..J: use frame of M_j (T_J: frame of M_{J-1})
    Ts = [_frame(Mm[min(j, J - 1)], eta) for j in range(J + 1)]
    Ta = [_to_acb(T) for T in Ts]
    Tinv = [T.inv() for T in Ta]
    # initial
    c = max(up(abs(Tinv[0][r, 0])*M0 + abs(Tinv[0][r, 1])*M1) for r in range(2))
    out = [M1, M0]
    for j in range(J):
        E = Tinv[j + 1]*Ms[j]*Ta[j]
        c = c*_inf_norm(E)
        row = up(abs(Ta[j + 1][0, 0]) + abs(Ta[j + 1][0, 1]))
        out.append(min(row*c, 1e300))
    return out


def _schur_frame(Mmid, eta):
    """unitary Schur basis (larger eigenvalue first) with the second column scaled by eta."""
    w, V = np.linalg.eig(Mmid)
    v = V[:, np.argmax(np.abs(w))]; v = v/np.linalg.norm(v)
    q = np.array([-np.conj(v[1]), np.conj(v[0])])
    return np.column_stack([v, q]) @ np.diag([1.0, eta])


def chain_core(Ms, Ts, M0, M1, J):
    Ta = [_to_acb(T) for T in Ts]
    Tinv = [T.inv() for T in Ta]
    c = max(up(abs(Tinv[0][r, 0])*M0 + abs(Tinv[0][r, 1])*M1) for r in range(2))
    out = [M1, M0]
    for j in range(J):
        E = Tinv[j + 1]*Ms[j]*Ta[j]
        c = c*_inf_norm(E)
        row = up(abs(Ta[j + 1][0, 0]) + abs(Ta[j + 1][0, 1]))
        out.append(min(row*c, 1e300))
    return out


def chain_bounds_multi(sig0, t_lo, t_hi, M0, M1, J, d_lo, d_hi, a_lo, a_hi, s_list=(1, 9)):
    """elementwise minimum of Lemma B bounds over several frame families (each one rigorous): the default family of chain_bounds,
    and scaled Schur bases T_j = U_j diag(1, eta_j) with eta_0 = 1, eta_{j+1} = max(eta_min, q eta_j) (no initial conversion cost)."""
    d, g, s = state_box(a_lo, a_hi, d_lo, d_hi, s_list)
    tau = arb_iv(t_lo, t_hi)
    Ms = []; Mm = []
    for j in range(J):
        beta = acb(arb(sig0 - j), tau)
        A, B = coef(beta, d, g, s)
        Ms.append(acb_mat([[A, B], [acb(1), acb(0)]]))
        Mm.append(np.array([[_mid(A), _mid(B)], [1, 0]], dtype=complex))
    best = chain_core(Ms, [_frame(Mm[min(j, J - 1)]) for j in range(J + 1)], M0, M1, J)
    for q in (0.5, 0.7):
        for emin in (0.03, 0.15):
            eta = 1.0; Ts = []
            for j in range(J + 1):
                Ts.append(_schur_frame(Mm[min(j, J - 1)], eta)); eta = max(emin, q*eta)
            b = chain_core(Ms, Ts, M0, M1, J)
            best = [min(x, y) for x, y in zip(best, b)]
    return best


if __name__ == '__main__':
    # sanity vs exact (probe3 numbers): a = 998, L = 4.5, tau = 1031.5: |S^_j| = 1546, 808, 420, 218, 112, 57.8, 29.6, ..., 3.87 (j=9)
    d = math.exp(-4.5)
    b = chain_bounds(2.5, 1031.0, 1032.0, 1546.5*1.0, 2960.0, 40, d*0.999, d*1.001, 998.0, 999.0)
    print(' '.join(f'{x:.3g}' for x in b[:20]))
