"""
Direct numerical test of the Lorentzian shadow representation (no large-Delta limit, no functional):
   scalar seeds:    K g^{a,b}_{Delta,J}  =?  (-1)^J/(4 b^{(D34)}) int_{x3<x5<x4} d^3x5 int D n  S3(1,2;5,n) S3^S(5,n;3,4)
   spinning:        W_diff               =?  (-1)^J/(4 b^{(1-Dphi)}) int d^3x5 int D n  L(z_1) Fs(z_3)
in the symmetric frame (x_1,x_2,x_3,x_4) = (-y, x/x^2, y, -x/x^2).  4-dimensional Gauss quadrature.
"""
import sys, numpy as np, mpmath as mp, sympy as sp
from blocks import CBlock
from spinblock import SpinBlock, invariants, block_coords_from_uv
from structures import rel_structure, logMref, Dphi as symDphi, get as get_struct
from embed import Delta as symDelta, J as symJ
from formulas import log_b
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def md(u, v): return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]

def diamond_nodes(x3, x4, nuv=70, nphi=40):
    """nodes and weights for int d^3x5 over the causal diamond x3 < x5 < x4"""
    d = x4 - x3; L = np.sqrt(-md(d, d)); e0 = d/L
    # orthonormal spacelike vectors orthogonal to e0
    a = np.array([0.0, 0.0, 1.0]); a = a + md(a, e0)*e0; e2 = a/np.sqrt(md(a, a))
    b = np.array([0.0, 1.0, 0.0]); b = b + md(b, e0)*e0; b = b - md(b, e2)*e2; e1 = b/np.sqrt(md(b, b))
    g, w = np.polynomial.legendre.leggauss(nuv)
    # u_ = t - r, v_ = t + r in [-L/2, L/2], u_ <= v_ ;  substitution concentrating nodes at the ends
    U = 0.5*L*g; V = 0.5*L*g
    UU, VV = np.meshgrid(U, V, indexing='ij'); WW = np.outer(w, w)*(0.5*L)**2
    mask = VV > UU
    t = 0.5*(UU + VV)[mask]; r = 0.5*(VV - UU)[mask]; wuv = WW[mask]*0.5*r       # d^3x = r dr dt dphi, dr dt = du dv/2
    phi = 2*np.pi*(np.arange(nphi) + 0.5)/nphi
    c = 0.5*(x3 + x4)
    X5 = c[None, None, :] + t[:, None, None]*e0 + r[:, None, None]*(np.cos(phi)[None, :, None]*e1 + np.sin(phi)[None, :, None]*e2)
    W5 = wuv[:, None]*(2*np.pi/nphi)*np.ones(nphi)[None, :]
    return X5.reshape(-1, 3), W5.reshape(-1)

def frame(q, xi):
    y = np.array([q, 0.0, 0.0]); x = q*np.array([np.cosh(xi), np.sinh(xi), 0.0])
    x2 = x/md(x, x)
    return -y, x2, y, -x2

def seed_struct():
    s12, s15, s25, z51, z52 = sp.symbols('s12 s15 s25 z51 z52', positive=True)
    return None

def S3num(Dl, Jl, d1, d2, x1, x2, x5, n):
    s12 = -md(x1-x2, x1-x2); s15 = -md(x1-x5, x1-x5); s25 = -md(x2-x5, x2-x5)
    v5 = 0.5*(md(n, x1-x5)*s25 - md(n, x2-x5)*s15)
    al = (d1 + d2 - Dl + Jl)/2; be = (Dl + Jl + d1 - d2)/2; ga = (Dl + Jl - d1 + d2)/2
    return np.exp(Jl*np.log(np.abs(v5)) - al*np.log(s12) - be*np.log(s15) - ga*np.log(s25)), np.sign(v5)

def scalar_test(Delta, J, d, q, xi, N=400, nn=64):
    d1, d2, d3, d4 = d
    x1, x2, x3, x4 = frame(q, xi)
    X5, W5 = diamond_nodes(x3, x4)
    tot = 0.0
    sg = set()
    for k in range(nn):
        ph = 2*np.pi*(k + 0.5)/nn
        n = np.array([1.0, np.cos(ph), np.sin(ph)])
        SL, sL = S3num(Delta, J, d1, d2, x1[None], x2[None], X5, n)
        SR, sR = S3num(3 - Delta, -1 - J, d3, d4, x3[None], x4[None], X5, n)
        tot += np.sum(W5*SL*SR)*(2*np.pi/nn)
        sg |= set(np.unique(sL)) | set(np.unique(sR))
    a, b = -(d1 - d2)/2, (d3 - d4)/2
    sh = float((-1)**J*mp.exp(-log_b(mp.mpf(Delta), J, d3 - d4))/4)*tot
    # K g
    s = {'s12': -md(x1-x2, x1-x2), 's13': -md(x1-x3, x1-x3), 's14': -md(x1-x4, x1-x4), 's23': -md(x2-x3, x2-x3), 's24': -md(x2-x4, x2-x4), 's34': -md(x3-x4, x3-x4)}
    u = s['s12']*s['s34']/(s['s13']*s['s24']); v = s['s14']*s['s23']/(s['s13']*s['s24'])
    A0 = 1/np.sqrt(u); E0 = np.sqrt(v/u)
    r = np.exp(-np.arccosh(A0 + E0)); chi = np.arccosh(E0 - A0)
    blk = CBlock(Delta, J, -a, b, N)
    gt = (4*r)**Delta*blk.series([r], [chi], sig=+1)[0, 0, 0]
    g = (-1)**J*v**(-b)*gt
    K = (s['s24']/s['s14'])**((d1-d2)/2)*(s['s14']/s['s13'])**((d3-d4)/2)/(s['s12']**((d1+d2)/2)*s['s34']**((d3+d4)/2))
    return K*g, sh, sg

_sf = {}
def struct_np(l):
    if l not in _sf:
        rs = rel_structure(l); sy, V5 = rs[0][2], rs[0][3]
        cs, se = get_struct(l)
        expr = sum(cs[w]*Q for (w, Q, _, _) in rs)
        names = ['s12', 's15', 's25', 'z12', 'z15', 'z51', 'z52', 'z1z5']
        args = [sy[n] for n in names] + [symDelta, symJ, symDphi]
        f = sp.lambdify(args, expr, modules='numpy')
        fl = sp.lambdify(args, logMref(l, sy, V5).subs(sp.log(V5), sp.log(sp.Abs(V5))), modules='numpy')
        _sf[l] = (f, fl)
    f, fl = _sf[l]
    def F(D, Jv, Dp, x1, x2, x5, z1, n):
        s12 = -md(x1-x2, x1-x2); s15 = -md(x1-x5, x1-x5); s25 = -md(x2-x5, x2-x5)
        a = [s12, s15, s25, md(z1, x2-x1), md(z1, x5-x1), md(n, x1-x5), md(n, x2-x5), md(z1, n), D, Jv, Dp]
        return np.exp(fl(*a))*f(*a)
    return F

def spinning_test(l, Delta, J, Dphi, q, xi, th1, th3, N=400, nn=64):
    x1, x2, x3, x4 = frame(q, xi)
    z1 = np.array([1.0, np.cos(th1), np.sin(th1)]); z3 = np.array([1.0, np.cos(th3), np.sin(th3)])
    X5, W5 = diamond_nodes(x3, x4)
    F = struct_np(l)
    tot = 0.0
    for k in range(nn):
        ph = 2*np.pi*(k + 0.5)/nn
        n = np.array([1.0, np.cos(ph), np.sin(ph)])
        Lv = F(Delta, J, Dphi, x1[None], x2[None], X5, z1[None], n)
        Rv = F(3 - Delta, -1 - J, Dphi, x3[None], x4[None], X5, z3[None], n)
        tot += np.sum(W5*Lv*Rv)*(2*np.pi/nn)
    sh = float((-1)**J*mp.exp(-log_b(mp.mpf(Delta), J, 1 - Dphi))/4)*tot
    sb = SpinBlock(l, Delta, J, Dphi, N)
    s, z = invariants(x1[None], x2[None], x3[None], x4[None], z1[None], z3[None])
    r, chi, A0, E0 = block_coords_from_uv(s)
    ls, M = sb.eval(s, z, r, chi, A0, E0)
    W = np.exp(ls - (l+1)*np.log(s['s13']) - Dphi*np.log(s['s24']))*M
    return W[0], sh

if __name__ == "__main__":
    pr("# spinning partial waves: W_diff (differential operators on exact blocks)  vs  shadow integral of L(z_1) Fs(z_3)")
    for l in (1, 2):
        for (Delta, J, Dphi) in ((7.3, 2, 1.75), (8.6, 3, 1.75), (9.1, 4, 2.5)):
            for (th1, th3) in ((0.3, 1.9), (2.2, 0.4), (4.0, 5.1)):
                W, sh = spinning_test(l, Delta, J, Dphi, 0.25, 0.5, th1, th3)
                pr(f"  l={l} Delta={Delta} J={J} Dphi={Dphi} pol=({th1},{th3}):  W_diff = {W:.8e}   shadow = {sh:.8e}   ratio = {W/sh:.6f}")
    pr("# scalar seeds: K g  vs  shadow integral  (q = |x| = |y|, xi);  third column: ratio")
    for (Delta, J) in ((7.3, 2), (8.6, 3)):
        for d in ((1.75, 1.75, 1.75, 1.75), (3.0, 1.75, 3.0, 1.75), (3.0, 1.75, 2.0, 2.75), (4.0, 1.75, 2.0, 3.75), (2.0, 3.75, 4.0, 1.75)):
            for (q, xi) in ((0.25, 0.5),):
                Kg, sh, sg = scalar_test(Delta, J, d, q, xi)
                pr(f"  Delta={Delta} J={J} dims={d} q={q} xi={xi}:  K g = {Kg:.8e}   shadow = {sh:.8e}   ratio = {Kg/sh:.6f}   signs of v5: {sorted(sg)}")
    open('out_t5_shadow_direct.txt', 'w').write("\n".join(out) + "\n")
