"""
Conservation of <O_1(x) O_2(0) X(inf, zeta)> in the invariants  AX = a.x, BX = b.x, WX = zeta.x, XX = x^2, AB, AW, BW  (a, b, zeta null),
exact at rational (Delta, J).  F = WX^J XX^{-al} R (constant (-1/2)^J dropped).
Divergence at x_1 (polarisation a; for spin 2 of the harmonic part h = R - (a.a) Lap_a R / 6), by the chain rule:
   d_{a_mu} f = f_AX x_mu + f_AB b_mu + f_AW w_mu (+ 2 a_mu f_AA),   d_{x_mu} g = g_AX a_mu + g_BX b_mu + g_WX w_mu + 2 g_XX x_mu,
   sum_mu d_mu(g x_mu) = 3g + x.dg,  x.dg = g_AX AX + g_BX BX + g_WX WX + 2 g_XX XX,  b.dg = g_AX AB + g_WX BW + 2 g_XX BX,
   w.dg = g_AX AW + g_BX BW + 2 g_XX WX,  a.dg = g_AX AA + g_BX AB + g_WX AW + 2 g_XX AX,  d log(WX^J XX^-al) = J w/WX - 2 al x/XX;
a.a = AA kept until the end, then set to 0.  At x_2: the same with a <-> b and an overall sign.
The result must vanish for all null a, b, zeta and x in d = 3: evaluated at random exact 3d configurations (Gaussian rationals).
"""
import random
import sympy as sp
from ft_helicity import R_inv, AX, BX, WX, XX, AB, AW, BW
from infframe import words
AA, BB = sp.symbols('AA BB')

def _dirs(g, which):
    gAX, gBX, gWX, gXX = (sp.diff(g, s) for s in (AX, BX, WX, XX))
    xd = gAX*AX + gBX*BX + gWX*WX + 2*gXX*XX
    bd = gAX*AB + gBX*BB + gWX*BW + 2*gXX*BX
    wd = gAX*AW + gBX*BW + 2*gXX*WX
    ad = gAX*AA + gBX*AB + gWX*AW + 2*gXX*AX
    return xd, ad, bd, wd

def divergence(R, l, al, Jv, pt):
    if pt == 2:      # swap roles a <-> b
        swap = {AX: BX, BX: AX, AW: BW, BW: AW}
        R = R.xreplace(swap)
    if l == 2:
        lap = (sp.diff(R, AX, 2)*XX + 2*sp.diff(R, AX, AB)*BX + 2*sp.diff(R, AX, AW)*WX + 2*sp.diff(R, AB, AW)*BW
               + sp.diff(R, AB, 2)*BB)
        h = R - AA*lap/6
    else:
        h = R
    hAX, hAB, hAW, hAA = (sp.diff(h, s) for s in (AX, AB, AW, AA))
    # g_mu = hAX x_mu + hAB b_mu + hAW w_mu + 2 hAA a_mu ;  div = sum_mu (d_mu + dlogG_mu) g_mu
    div = 3*hAX + _dirs(hAX, "x")[0] + _dirs(hAB, "b")[2] + _dirs(hAW, "w")[3] + 2*_dirs(hAA, "a")[1]
    # d log G . g
    xg = hAX*XX + hAB*BX + hAW*WX + 2*hAA*AX
    wg = hAX*WX + hAB*BW + 2*hAA*AW
    div += Jv*wg/WX - 2*al*xg/XX
    div = div.xreplace({AA: 0, BB: 0})
    if pt == 2:
        div = -div.xreplace({AX: BX, BX: AX, AW: BW, BW: AW})
    return div

TRIP = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (9, 40, 41), (12, 35, 37), (11, 60, 61)]
def null_vec(rng):
    a, b, c = rng.choice(TRIP); v = [sp.Integer(a), sp.Integer(rng.choice([1, -1])*b), sp.I*c]; rng.shuffle(v)
    k = sp.Rational(1, rng.randint(1, 3))          # ONE scale per vector (keeps it null)
    return [x*k for x in v]

def configs(n, seed=1):
    rng = random.Random(seed)
    d = lambda u, v: sum(u[i]*v[i] for i in range(3))
    out = []
    for _ in range(n):
        x = [sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)) or 1 for _ in range(3)]
        a, b, w = null_vec(rng), null_vec(rng), null_vec(rng)
        out.append({AX: d(a, x), BX: d(b, x), WX: d(w, x), XX: d(x, x), AB: d(a, b), AW: d(a, w), BW: d(b, w)})
    return out

def conserved_space(l, Dv, Jv, ncfg=None):
    ws = words(l)
    al = sp.Rational(4*l + 2, 2) - sp.nsimplify(Dv)/2 + sp.Rational(Jv, 2)
    cfg = configs(ncfg or (8 if l == 1 else 16))
    rows = []
    cols = []
    for w in ws:
        R = R_inv(l, {w: 1}, Dv, Jv)
        e = [divergence(R, l, al, Jv, pt) for pt in (1, 2)]
        cols.append(e)
    for c in cfg:
        for pt in range(2):
            rows.append([sp.expand(sp.nsimplify(cols[j][pt].xreplace(c))) for j in range(len(ws))])
    M = sp.Matrix(rows)
    # word values themselves (to find the d = 3 relations)
    vals = sp.Matrix([[sp.expand(R_inv(l, {w: 1}, Dv, Jv).xreplace(c)) for w in ws] for c in cfg])
    return ws, M.nullspace(), vals.nullspace()
