"""
NOTE: the background/first-law half (general N(r), f(r) curvature + high-order Euler-Lagrange + dsolve) did not finish
within 40 min and was stopped; out_gradient_direct.txt contains the shear-side results (e1, U identity) only.
The entropy side is done instead by Iyer-Wald (boxC_wald.py) and cross-checked by the identity route (README item 2/3).

Items 2-3: direct (no integration-by-parts shortcut, no on-shell identity) shear responses for derivative operators.
  dV_E = e1 - k - 2 tau1, where
     e1   = eta_1/eta_0 at fixed horizon r_h = 1 (membrane flux of the full quadratic action, up to phi''' jets),
     k    = int_1^oo r^2 O dr   (first-order on-shell action: f = f0 (1 + k lam) at fixed T),
     tau1 = temperature shift at fixed r_h (O(lam) background from the reduced action),
  using the first law (w = k + 2 tau1 for the fixed-r_h entropy), so no Wald formula for derivative operators is assumed.
  For the Wald check, the claimed w = 4 A_I(1) + 6 A_J(1) is compared with k + 2 tau1.            [computed]
Operators: (grad I)^2, (grad J)^2, U = grad C . grad C, I U; plus the on-shell identity U = 1/2 box I + 6 I + 3 J
on the background.
"""
import sys
import sympy as sp
from static_engine import *

f0 = r ** 2 - 1 / r
lam = sp.Symbol('lam')

# ---------- quadratic (shear) side on the Einstein brane
gP = metric(1, f0)
gi, Gam, Rl, Ric, Rs, W = curvature(gP)
I, J = invariants(gP, gi, W)
Is, Js, Rss = eps_series(I), eps_series(J), eps_series(Rs)
rsp = 1 / f0
eta0 = flux_at_horizon(reduce_quadratic(r ** 2 * Rss[2]), rsp)
print("eta0 =", eta0)
U = gradC2(gP, gi, Gam, W)
Us = eps_series(U)
print("U background =", sp.factor(Us[0]), "   U first order =", Us[1])
boxI = sp.simplify(box(metric(1, f0, pert=False), sp.diag(-1 / f0, f0, r ** -2, r ** -2), Is[0]))
print("box I (background) =", sp.expand(boxI), ";  1/2 box I + 6 I + 3 J - U =", sp.simplify(boxI / 2 + 6 * Is[0] + 3 * Js[0] - Us[0]))


def quad_e1(L2):
    al = reduce_quadratic(L2)
    return sp.simplify(flux_at_horizon(al, rsp) / eta0), al


ops2 = {}
# (grad I)^2 = f (D_r I)^2 in the static diagonal metric (g^rr = f exactly; the shear does not touch g^rr)
DI = [Dr(Is[k]) for k in range(3)]
DJ = [Dr(Js[k]) for k in range(3)]
ops2['(gradI)^2'] = r ** 2 * f0 * (2 * DI[0] * DI[2] + DI[1] ** 2)
ops2['(gradJ)^2'] = r ** 2 * f0 * (2 * DJ[0] * DJ[2] + DJ[1] ** 2)
ops2['U'] = r ** 2 * Us[2]
ops2['I U'] = r ** 2 * (Is[0] * Us[2] + Is[2] * Us[0] + Is[1] * Us[1])
E1 = {}
for k_, L2 in ops2.items():
    E1[k_], al = quad_e1(sp.expand(L2))
    print(f"[{k_}] alphas = {{{', '.join(f'{a}: {sp.factor(b)}' for a, b in al.items())}}}  ->  e1 = {E1[k_]}")

# ---------- background side: O(lam) solution, tau1, k
N, f = sp.Function('N')(r), sp.Function('f')(r)
gB = sp.diag(-N ** 2 * f, 1 / f, r ** 2, r ** 2)
giB, GamB, RlB, RicB, RsB, WB = curvature(gB)
IB, JB = invariants(gB, giB, WB)
UB = gradC2(gB, giB, GamB, WB) if any(w in ('U', 'I U') for w in (sys.argv[1:] or ['U'])) else None
XII = f * sp.diff(IB, r) ** 2
XJJ = f * sp.diff(JB, r) ** 2
bg_ops = {'(gradI)^2': XII, '(gradJ)^2': XJJ, 'U': UB, 'I U': (IB * UB if UB is not None else None)}


def EL(Lag, F):
    out = 0
    for kk in range(0, 7):
        Dk = F if kk == 0 else sp.diff(F, r, kk)
        out += (-1) ** kk * sp.diff(sp.diff(Lag, Dk), r, kk)
    return out


f1, n1 = sp.Function('f1')(r), sp.Function('n1')(r)
# claimed Wald (Euler-derivative) coefficients: A_I, A_J for F(y, X)
boxB = lambda s: sp.diff(r ** 2 * f0 * sp.diff(s, r), r) / r ** 2
A = {'(gradI)^2': (-2 * boxB(Is[0]), 0), '(gradJ)^2': (0, -2 * boxB(Js[0]))}
which = [w.replace('_', ' ') for w in sys.argv[1:]] or list(bg_ops)
for name in which:
    O = bg_ops[name]
    Lag = N * r ** 2 * (RsB + 6 + lam * O)
    sub = {N: 1 + lam * n1, f: f0 + lam * f1}
    eqs = [sp.simplify(sp.diff(EL(Lag, F).subs(sub).doit(), lam).subs(lam, 0)) for F in (N, f)]
    n1s = sp.dsolve(eqs[1], n1).rhs.subs(sp.Symbol('C1'), sp.Symbol('cn'))
    f1s = sp.dsolve(sp.simplify(eqs[0].subs(n1, n1s).doit()), f1).rhs
    C = sorted((f1s.free_symbols | n1s.free_symbols) - {r, lam}, key=str)
    sol = sp.solve([sp.Eq(f1s.subs(r, 1), 0), sp.Eq(sp.limit((2 * n1s * f0 + f1s) / r ** 2, r, sp.oo), 0)], C, dict=True)[0]
    f1s, n1s = sp.simplify(f1s.subs(sol)), sp.simplify(n1s.subs(sol))
    tau1 = sp.simplify(sp.limit((n1s * sp.diff(f0, r) + sp.diff(f1s, r)) / sp.diff(f0, r), r, 1))
    k_ = sp.simplify(sp.integrate(sp.expand(r ** 2 * O.subs({N: 1, f: f0}).doit()), (r, 1, sp.oo)))
    dVE = sp.simplify(E1[name] - k_ - 2 * tau1)
    line = f"[{name}] f1 = {f1s}, n1 = {n1s};  tau1 = {tau1}, k = {k_};  e1 = {E1[name]}  =>  dV_E = {dVE}"
    if name in A:
        AI, AJ = A[name]
        w_claim = sp.simplify(4 * sp.sympify(AI).subs(r, 1) + 6 * sp.sympify(AJ).subs(r, 1))
        dV_claim = sp.simplify(-108 * sp.sympify(AJ).subs(r, 1) - 6 * sp.diff(sp.sympify(AI), r).subs(r, 1)
                               + 18 * sp.diff(sp.sympify(AJ), r).subs(r, 1))
        line += f";  claimed Wald w = {w_claim} vs first-law k + 2 tau1 = {k_ + 2 * tau1};  claimed dV_E = {dV_claim}"
    print(line)
