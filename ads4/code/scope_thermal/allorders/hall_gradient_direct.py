"""
Item 2 (Hall): direct EF computation for grad I . grad P (no integration by parts by hand).        [computed]
L2 = sqrt(-g) * [ g^{ab} d_a I0 d_b P^(2) ]  with I0 = 12/r^6 (v-independent), g^{vr} = 1, g^{rr} = f, sqrt(-g) = r^2:
     = r^2 I0' ( D_v P^(2) + f D_r P^(2) ),  D = total derivatives on jets (third-order jets introduced).
Euler-Lagrange with third-order jets, O(omega) classification as in hall_ef.py; claimed 4 pi eta_H/s = -15552.
"""
import io, contextlib
import sympy as sp
with contextlib.redirect_stdout(io.StringIO()):
    import ef_engine as E
r, v, phi, chi = E.r, E.v, E.phi, E.chi
J = dict(E.J)
names = []
for nm in ('p', 'q'):
    for k in ('vvv', 'vvr', 'vrr', 'rrr'):
        J[nm + k] = sp.Symbol(nm + k)
order = ['', 'v', 'r', 'vv', 'vr', 'rr', 'vvv', 'vvr', 'vrr', 'rrr']


def key(s):  # canonical sort of derivative label
    return ''.join(sorted(s, key=lambda ch: 0 if ch == 'v' else 1))


def D(expr, var):
    out = sp.diff(expr, var)
    for nm in ('p', 'q'):
        for lab in order[:6]:
            sym = J[nm + lab]
            new = key(lab + ('v' if var == v else 'r'))
            out += J[nm + new] * sp.diff(expr, sym)
    return sp.expand(out)


P2 = sp.expand(E.P_[2])
I0p = sp.diff(12 / r ** 6, r)
f = E.f
L2 = sp.expand(r ** 2 * I0p * (D(P2, v) + f * D(P2, r)))
L2c = E.cross_part(L2) if False else L2


def EL3(L2, field):
    F = phi if field == 'phi' else chi
    nm = 'p' if field == 'phi' else 'q'
    rev = {}
    for fn, n_ in ((phi, 'p'), (chi, 'q')):
        for lab in order:
            vs = [v if ch == 'v' else r for ch in lab]
            rev[J[n_ + lab]] = sp.diff(fn, *vs) if vs else fn
    out = 0
    for lab in order:
        if not L2.has(J[nm + lab]):
            continue
        dl = sp.diff(L2, J[nm + lab]).subs(rev)
        for ch in lab:
            dl = -sp.diff(dl, v if ch == 'v' else r)
        out += dl
    return sp.expand(out)


Echi = EL3(L2, 'chi')
PV = sp.Derivative(phi, v)
sig, bad, w2 = 0, [], 0
for t in sp.Add.make_args(Echi):
    if not t.has(phi):
        continue
    ders = [a for a in t.atoms(sp.Derivative) if a.expr == phi]
    if ders == [PV] and not (t / PV).has(phi):
        sig += t / PV
    elif not ders or all(all(vv == r for vv in a.variables) for a in ders):
        bad.append(t)
    else:
        w2 += 1
sig = sp.simplify(sig)
print("sigma(r) =", sig)
print("terms that would be O(1) or O(w) without phi_v (must be absent):", bad, ";  O(w^2) terms dropped:", w2)
print("4 pi eta_H/s = int_1^oo sigma =", sp.integrate(sig, (r, 1, sp.oo)), "  (claim -15552)")
