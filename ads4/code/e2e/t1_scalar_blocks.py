"""Validation of the scalar blocks with general external dimensions (task 1, scalar seeds)."""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, '../reggemoment')
from blocks import CBlock, block_mp
mp.mp.dps = 34
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

pr("# (1) a = b = 0: C (double) vs the recursion of reggemoment/blocks3d.py (independent derivation: depth 3 in r^2)")
from blocks3d import Block3d
for (Delta, J, rho, rhob) in ((7.3, 2, 0.3, 0.1), (20.3, 6, -0.9, -0.8), (40.3, 12, -0.95, -0.93)):
    N = 400 if abs(rho) < 0.5 else 2000
    g1 = CBlock(Delta, J, 0, 0, N).g(rho, rhob)[0]
    g0 = Block3d(Delta, J, N).g(rho, rhob, tol=1e-17)[0]
    pr(f"  Delta={Delta} J={J} rho={rho},{rhob}:  new {g1:.15e}  old {g0:.15e}  rel.diff {abs(g1/g0-1):.1e}")

pr("# (2) general (a, b): C double vs C quad vs independent mpmath recursion")
for (Delta, J, a, b, rho, rhob, N) in ((5.3, 2, 0.5, -0.5, 0.31, 0.12, 60), (9.75, 3, -0.625, 0.375, -0.45, -0.2, 70), (6.5, 0, 1.0, -0.25, -0.5, -0.33, 70),
                                       (12.25, 5, 0.625, -1.625, 0.4, 0.37, 80)):
    gd = CBlock(Delta, J, a, b, N).g(rho, rhob)[0]
    gq = CBlock(Delta, J, a, b, N, prec='q').g(rho, rhob)[0]
    gm = block_mp(Delta, J, a, b, N, rho, rhob)
    pr(f"  Delta={Delta} J={J} a={a} b={b} rho={rho},{rhob}: quad {mp.nstr(gq, 25)}  |double/quad-1| {float(abs(gd/gq-1)):.1e}  |mp/quad-1| {float(abs(gm/gq-1)):.1e}")

pr("# (3) Casimir equation in (z, zbar), Dolan-Osborn form, by numerical differentiation of the quad block (residual / |C g/2|)")
def gz(blk, z, zb):
    rho = z/(1+mp.sqrt(1-z))**2; rhob = zb/(1+mp.sqrt(1-zb))**2
    r = mp.sqrt(rho*rhob); chi = abs(mp.log(rho/rhob))/2
    sig = -1 if rho < 0 else 1
    S = blk.series([r], [chi], sig=sig)[0, 0, 0]
    return ((-1)**blk.J if sig < 0 else 1)*(4*r)**mp.mpf(blk.Delta)*S
for (Delta, J, a, b, z, zb, N) in ((5.3, 2, 0.5, -0.5, 0.4, 0.25, 160), (9.75, 3, -0.625, 0.375, -1.3, -0.6, 400), (6.5, 0, 1.0, -0.25, -2.5, -1.7, 500),
                                   (12.25, 5, 0.625, -1.625, -0.8, -0.3, 300)):
    blk = CBlock(Delta, J, a, b, N, prec='q')
    f = lambda x, y: gz(blk, x, y)
    h = mp.mpf('1e-5')
    def d(fun, i, order):
        return mp.diff(fun, (z, zb), tuple(order if q == i else 0 for q in range(2)), h=h)
    g0 = f(mp.mpf(z), mp.mpf(zb))
    gz1, gz2 = d(f, 0, 1), d(f, 0, 2); gb1, gb2 = d(f, 1, 1), d(f, 1, 2)
    zz, zzb = mp.mpf(z), mp.mpf(zb)
    Dz = zz**2*(1-zz)*gz2 - (1+a+b)*zz**2*gz1 - a*b*zz*g0
    Db = zzb**2*(1-zzb)*gb2 - (1+a+b)*zzb**2*gb1 - a*b*zzb*g0
    cas = Dz + Db + zz*zzb/(zz-zzb)*((1-zz)*gz1 - (1-zzb)*gb1)
    C = mp.mpf(Delta)*(mp.mpf(Delta)-3) + J*(J+1)
    pr(f"  Delta={Delta} J={J} a={a} b={b} z={z},{zb}:  g = {mp.nstr(g0, 15)}   residual {mp.nstr(abs(cas - C*g0/2)/abs(C*g0/2), 3)}")

pr("# (4) leading small-z term: g / [z^{(Delta-J)/2} k_{Delta+J}(zbar)] -> 1,  k_beta(x) = x^{beta/2} 2F1(beta/2+a, beta/2+b; beta; x)   (convention check for a, b)")
for (Delta, J, a, b) in ((5.3, 2, 0.5, -0.5), (9.75, 3, -0.625, 0.375), (6.5, 0, 1.0, -0.25)):
    blk = CBlock(Delta, J, a, b, 300, prec='q')
    row = []
    for z in ('1e-4', '1e-6', '1e-8'):
        z = mp.mpf(z); zb = mp.mpf('0.35')
        bt = mp.mpf(Delta) + J
        k = zb**(bt/2)*mp.hyp2f1(bt/2 + a, bt/2 + b, bt, zb)
        row.append(mp.nstr(gz(blk, z, zb)/(z**((mp.mpf(Delta)-J)/2)*k) - 1, 3))
    pr(f"  Delta={Delta} J={J} a={a} b={b}:  ratio-1 at z=1e-4,1e-6,1e-8: {row}")

pr("# (5) derivative sums: (r d_r)^k d_eta^l of the series vs numerical differentiation (quad), at rho,rhobar<0")
blk = CBlock(9.75, 3, -0.625, 0.375, 500, prec='q')
r0, chi0 = mp.mpf('0.8'), mp.mpf('0.1')
S = blk.series([r0], [chi0], sig=-1, kmax=4, lmax=4)[0]
def Ffun(L, eta):   # r^Delta sum B (-1)^n r^n P_j(x) as a function of L = log r and x = |eta|
    r = mp.exp(L); chi = mp.acosh(eta)
    return r**mp.mpf(blk.Delta)*blk.series([r], [chi], sig=-1)[0, 0, 0]
F0 = Ffun(mp.log(r0), mp.cosh(chi0))
row = []
for dl in ('4e-3', '2e-3', '1e-3'):
    dl = mp.mpf(dl); dL, de = dl, -dl/3
    tay = sum(r0**mp.mpf(blk.Delta)*S[k, l]*dL**k*de**l/(mp.factorial(k)*mp.factorial(l)) for k in range(5) for l in range(5 - k))
    row.append(mp.nstr((Ffun(mp.log(r0) + dL, mp.cosh(chi0) + de) - tay)/F0, 3))
pr(f"  (F(L+d, eta-d/3) - 4th-order Taylor polynomial built from the derivative sums)/F at d = 4e-3, 2e-3, 1e-3: {row}  (must scale as d^5)")
pr("# (6) convergence in N at a point of the dispersive region (Delta=40.3, J=12, a=0.5, b=-0.5; 1-|rho| = 0.6/Delta, 0.3/Delta)")
D = 40.3
ref = None
for N in (600, 1200, 1800, 2400):
    v = CBlock(D, 12, 0.5, -0.5, N).g(-(1-0.6/D), -(1-0.3/D), cut=60.0)[0]
    ref = v if N == 2400 else ref
    out.append((N, v))
for N, v in [o for o in out if isinstance(o, tuple)]:
    pr(f"  N={N}: g = {v:.14e}   rel. to N=2400: {abs(v/ref-1):.1e}")
out = [o for o in out if not isinstance(o, tuple)]
open('out_t1_scalar_blocks.txt', 'w').write("\n".join(out) + "\n")
