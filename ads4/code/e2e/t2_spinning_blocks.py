"""
Validation of the spinning partial waves W_diff for <J phi J phi> (l=1) and <T phi T phi> (l=2)  (task 1):
 (a) conservation at x_1 (Ward identity of the conserved operator) by finite differences in quad precision;
 (b) OPE limit x_12 -> 0 against an independent computation: leading OPE tensor of the left three-point structure
     contracted (spin-J index contraction on the null circle) with the right three-point structure;
     this fixes the normalisation  W[L, R] = (-1/8)^J W_diff  for unit-normalised O;
 (c) double versus quad pipelines.
Usage: python3 t2_spinning_blocks.py
"""
import sys, numpy as np, mpmath as mp, sympy as sp
from spinblock import SpinBlock, invariants, block_coords_from_uv, mdot, _mpfun
from structures import rel_structure, logMref, Dphi as symDphi, get as get_struct
from embed import Delta as symDelta, J as symJ
mp.mp.dps = 34
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def mpa(v): return np.array([mp.mpf(x) for x in np.ravel(v)], dtype=object).reshape(np.shape(v))
ETA = np.diag([-1, 1, 1])

def wfull(sb, x1, x2, x3, x4, z1, z3):
    """W_diff at one configuration (mp), polarisation arrays z1 [n1,3], z3 [n3,3] -> [n1,n3]"""
    obj = sb.prec == 'q'
    cv = (lambda v: mpa(v)) if obj else (lambda v: np.array(v, dtype=float))
    X = [cv(v)[None] for v in (x1, x2, x3, x4)]
    Z1 = cv(z1)[None, :, None, :] + 0*cv(z3)[None, None, :, :]; Z3 = cv(z3)[None, None, :, :] + 0*cv(z1)[None, :, None, :]
    s, z = invariants(X[0], X[1], X[2], X[3], Z1, Z3)
    r, chi, A0, E0 = block_coords_from_uv(s, obj=obj)
    ls, M = sb.eval(s, z, r, chi, A0, E0)
    ex = _mpfun(mp.exp) if obj else np.exp
    lg = _mpfun(mp.log) if obj else np.log
    pref = ex(ls - (sb.l+1)*lg(s['s13']) - sb.Dphi*lg(s['s24']))
    return (pref[:, None, None]*M)[0]

def tensor1(sb, x1, x2, x3, x4, z3, l):
    """covariant tensor components at point 1 (index-free in z3) reconstructed from null polarisations"""
    if l == 1:
        th = [mp.mpf(0), 2*mp.pi/3, 4*mp.pi/3]
    else:
        th = [2*mp.pi*k/5 for k in range(5)]
    z1 = [[mp.mpf(1), mp.cos(t), mp.sin(t)] for t in th]
    vals = wfull(sb, x1, x2, x3, x4, z1, [z3])[:, 0]
    if l == 1:
        A = mp.matrix(z1); sol = mp.lu_solve(A, mp.matrix(list(vals)))
        return [sol[i] for i in range(3)]
    # W(z) = W_mn z^m z^n, traceless: unknowns (W01, W02, W11, W12, W22), W00 = W11 + W22
    A = mp.matrix(5, 5)
    for k, zz in enumerate(z1):
        A[k, 0] = 2*zz[0]*zz[1]; A[k, 1] = 2*zz[0]*zz[2]; A[k, 2] = zz[1]**2 + zz[0]**2; A[k, 3] = 2*zz[1]*zz[2]; A[k, 4] = zz[2]**2 + zz[0]**2
    sol = mp.lu_solve(A, mp.matrix(list(vals)))
    W = mp.matrix(3, 3)
    W[0, 1] = W[1, 0] = sol[0]; W[0, 2] = W[2, 0] = sol[1]; W[1, 1] = sol[2]; W[1, 2] = W[2, 1] = sol[3]; W[2, 2] = sol[4]; W[0, 0] = sol[2] + sol[4]
    return W

def divergence(sb, x1, x2, x3, x4, z3, l, h=mp.mpf('1e-4')):
    def T(x):
        return tensor1(sb, x, x2, x3, x4, z3, l)
    def d(mu):
        e = [mp.mpf(0)]*3; e[mu] = h
        xs = [[x1[i] + k*e[i] for i in range(3)] for k in (-2, -1, 1, 2)]
        Ts = [T(x) for x in xs]
        if l == 1:
            return [(Ts[0][i] - 8*Ts[1][i] + 8*Ts[2][i] - Ts[3][i])/(12*h) for i in range(3)]
        return (Ts[0] - 8*Ts[1] + 8*Ts[2] - Ts[3])/(12*h)
    ds = [d(mu) for mu in range(3)]
    T0 = T(x1)
    if l == 1:
        div = -ds[0][0] + ds[1][1] + ds[2][2]
        grad = max(abs(ds[mu][i]) for mu in range(3) for i in range(3))
        return abs(div)/grad
    div = [-ds[0][0, nu] + ds[1][1, nu] + ds[2][2, nu] for nu in range(3)]
    grad = max(abs(ds[mu][i, j]) for mu in range(3) for i in range(3) for j in range(3))
    return max(abs(x) for x in div)/grad

y0 = [mp.mpf('0.5'), mp.mpf('0.05'), mp.mpf('0.1')]; xx = [mp.mpf('0.6'), mp.mpf('-0.2'), mp.mpf('0.1')]
x2n = sum((-xx[0]**2, xx[1]**2, xx[2]**2))
conf = dict(x1=[-v for v in y0], x2=[v/x2n for v in xx], x3=y0, x4=[-v/x2n for v in xx])
conf2 = dict(x1=[mp.mpf('-0.31'), mp.mpf('0.07'), mp.mpf('-0.02')], x2=[mp.mpf('-2.6'), mp.mpf('0.4'), mp.mpf('0.9')],
             x3=[mp.mpf('0.45'), mp.mpf('-0.11'), mp.mpf('0.06')], x4=[mp.mpf('3.1'), mp.mpf('1.1'), mp.mpf('-0.3')])
z3 = [mp.mpf(1), mp.cos(mp.mpf('1.9')), mp.sin(mp.mpf('1.9'))]

pr("# (a) conservation at x_1:  |d^mu W_{mu...}| / max |d_mu W_{nu...}|  (quad precision, 4th-order differences, h = 1e-4)")
for l, cases in ((1, [(9.75, 3, 1.75), (6.5, 4, 2.5)]), (2, [(9.75, 3, 1.75), (8.5, 4, 2.5)])):
    for (D, Jv, Dp) in cases:
        sb = SpinBlock(l, D, Jv, Dp, 260, prec='q')
        for nm, c in (('frame-like', conf), ('generic', conf2)):
            res = divergence(sb, c['x1'], c['x2'], c['x3'], c['x4'], z3, l)
            pr(f"  l={l} Delta={D} J={Jv} Dphi={Dp} [{nm}]: {mp.nstr(res, 3)}")
        if l == 1:
            # control: a non-conserved combination (coefficients swapped) must fail
            cw = dict(sb.c); sb.c = {'a': cw['b'], 'b': cw['a']}
            res = divergence(sb, conf['x1'], conf['x2'], conf['x3'], conf['x4'], z3, l)
            pr(f"     control (coefficients c_a <-> c_b swapped): {mp.nstr(res, 3)}")

pr("# (b) OPE limit x_12 = delta*v -> 0:  W_diff (-1/8)^J / <C_L, R>  - 1   (must vanish linearly in delta)")
def struct_fun(l, D, Jv, Dp):
    rs = rel_structure(l)
    sy, V5 = rs[0][2], rs[0][3]
    cs, se = get_struct(l)
    sub = {symDelta: sp.nsimplify(D), symJ: Jv, symDphi: sp.nsimplify(Dp)}
    expr = sum(cs[w]*Q for (w, Q, _, _) in rs).subs(sub)
    lm = logMref(l, sy, V5).subs(sub)
    names = ['s12', 's15', 's25', 'z12', 'z15', 'z51', 'z52', 'z1z5']
    f = sp.lambdify([sy[n] for n in names], expr, modules='mpmath')
    fl = sp.lambdify([sy[n] for n in names], sp.exp(lm), modules='mpmath')
    def L(x1, x2, x5, z1, z5):
        d = lambda u, v: -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
        df = lambda u, v: [u[i] - v[i] for i in range(3)]
        s12 = -d(df(x1, x2), df(x1, x2)); s15 = -d(df(x1, x5), df(x1, x5)); s25 = -d(df(x2, x5), df(x2, x5))
        args = [s12, s15, s25, d(z1, df(x2, x1)), d(z1, df(x5, x1)), d(z5, df(x1, x5)), d(z5, df(x2, x5)), d(z1, z5)]
        return fl(*args)*f(*args)
    return L
def refl(x, n):
    d = lambda u, v: -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    x2 = d(x, x); c = d(x, n)
    return [n[i] - 2*x[i]*c/x2 for i in range(3)]
def ope_check(l, D, Jv, Dp, delta):
    Lf = struct_fun(l, D, Jv, Dp)
    x2 = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    v = [mp.mpf('1.0'), mp.mpf('0.3'), mp.mpf('-0.2')]
    x1 = [x2[i] + delta*v[i] for i in range(3)]
    x3 = [mp.mpf('1.3'), mp.mpf('0.2'), mp.mpf('0.4')]; x4 = [mp.mpf('2.9'), mp.mpf('-0.3'), mp.mpf('0.5')]
    x5 = [mp.mpf('1e12'), mp.mpf('3e11'), mp.mpf('-2e11')]
    z1 = [mp.mpf(1), mp.cos(mp.mpf('0.7')), mp.sin(mp.mpf('0.7'))]
    z3 = [mp.mpf(1), mp.cos(mp.mpf('2.1')), mp.sin(mp.mpf('2.1'))]
    M = 2*Jv + 1
    phis = [2*mp.pi*k/M for k in range(M)]
    d = lambda u, w: -u[0]*w[0] + u[1]*w[1] + u[2]*w[2]
    x25 = [x2[i] - x5[i] for i in range(3)]
    s25 = -d(x25, x25)
    CL = []; RR = []
    for ph in phis:
        n = [mp.mpf(1), mp.cos(ph), mp.sin(ph)]
        CL.append(s25**mp.mpf(D)*Lf(x1, x2, x5, z1, refl(x25, n)))
        RR.append(Lf(x3, x4, x2, z3, n))         # O at x_2 with polarisation n
    am = {m: sum(CL[k]*mp.exp(-1j*m*phis[k]) for k in range(M))/M for m in range(-Jv, Jv+1)}
    bm = {m: sum(RR[k]*mp.exp(-1j*m*phis[k]) for k in range(M))/M for m in range(-Jv, Jv+1)}
    pair = sum(am[-m]*bm[m]*2**Jv*(-1)**(Jv-m)/mp.binomial(2*Jv, Jv-m) for m in range(-Jv, Jv+1))
    sb = SpinBlock(l, D, Jv, Dp, 60, prec='q') if l > 0 else None
    W = wfull(sb, x1, x2, x3, x4, [z1], [z3])[0, 0]
    return W*(-mp.mpf(1)/8)**Jv/pair - 1
for l, cases in ((1, [(9.75, 3, 1.75), (6.5, 4, 2.5)]), (2, [(9.75, 3, 1.75), (8.5, 4, 2.5)])):
    for (D, Jv, Dp) in cases:
        row = [mp.nstr(ope_check(l, D, Jv, Dp, mp.mpf(dl)), 4) for dl in ('1e-3', '1e-4', '1e-5')]
        pr(f"  l={l} Delta={D} J={Jv} Dphi={Dp}:  delta = 1e-3, 1e-4, 1e-5: {row}")

pr("# (c) double vs quad evaluation of the stripped tensor at a point of the dispersive region (frame, r_u = 0.002, eta_u = 1.4)")
for l, (D, Jv, Dp) in ((1, (30.25, 6, 1.75)), (2, (30.25, 6, 1.75))):
    q = np.sqrt(0.002); xi = np.arccosh(1.4)
    y = np.array([q, 0, 0]); xv = q*np.array([np.cosh(xi), np.sinh(xi), 0])
    x2v = xv/mdot(xv, xv)
    th = np.array([0.3, 1.1, 2.9]); zz = np.stack([np.ones(3), np.cos(th), np.sin(th)], axis=1)
    vals = {}
    for prec in ('d', 'q'):
        sb = SpinBlock(l, D, Jv, Dp, 1500, prec=prec)
        if prec == 'q':
            W = wfull(sb, [mp.mpf(float(v)) for v in -y], [mp.mpf(float(v)) for v in x2v], [mp.mpf(float(v)) for v in y], [mp.mpf(float(v)) for v in -x2v],
                      [[mp.mpf(float(c)) for c in row] for row in zz], [[mp.mpf(float(c)) for c in row] for row in zz])
        else:
            W = wfull(sb, -y, x2v, y, -x2v, zz, zz)
        vals[prec] = W
    rel = max(abs(float(vals['d'][i, j]/vals['q'][i, j] - 1)) for i in range(3) for j in range(3))
    pr(f"  l={l} Delta={D} J={Jv}: max |double/quad - 1| over 9 polarisation pairs = {rel:.2e};  sample value {mp.nstr(vals['q'][0,1], 12)}")
open('out_t2_spinning_blocks.txt', 'w').write("\n".join(out) + "\n")
