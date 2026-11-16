"""Validation of blocks3d.py.
 (1) float64 vectorised recursion vs independent mpmath implementation;
 (2) the block satisfies the quadratic Casimir PDE in (z, zbar) (finite differences, mpmath);
 (3) normalisation: g -> z^{(Delta-J)/2} zbar^{(Delta+J)/2} for z << zbar << 1;
 (4) convergence in N near |rho| -> 1 for both rho, rhobar < 0 (the region needed for the Regge moments).
"""
import numpy as np, mpmath as mp
from blocks3d import Block3d, block_mp

mp.mp.dps = 40
def rho_of_z(z):
    return z/(1+mp.sqrt(1-z))**2

print("(1) float64 vs mpmath")
for (Delta, J, rho, rhob) in ((7.3, 2, -0.31, -0.52), (21.7, 6, -0.62, -0.81), (33.1, 9, -0.2, -0.7), (12.4, 0, 0.35, 0.22)):
    N = 120
    blk = Block3d(Delta, J, N)
    a = blk.g(rho, rhob)[0]; b = block_mp(Delta, J, N, rho, rhob)
    print(f"  Delta={Delta} J={J} rho=({rho},{rhob}): float {a:.14e}  mp {mp.nstr(b, 15)}  rel.diff {float(abs(a-b)/abs(b)):.2e}")

print("(2) Casimir PDE residual (mpmath, N = 160)")
def gz(Delta, J, z, zb, N=160):
    return block_mp(Delta, J, N, rho_of_z(z), rho_of_z(zb))
for (Delta, J, z, zb) in ((5.3, 2, mp.mpf('-0.8'), mp.mpf('-2.1')), (8.7, 3, mp.mpf('0.21'), mp.mpf('0.43')), (6.1, 4, mp.mpf('-0.35'), mp.mpf('-0.9'))):
    f = lambda a, b: gz(Delta, J, a, b)
    h = mp.mpf('1e-8')
    fz = (f(z+h, zb) - f(z-h, zb))/(2*h); fzz = (f(z+h, zb) - 2*f(z, zb) + f(z-h, zb))/h**2
    fb = (f(z, zb+h) - f(z, zb-h))/(2*h); fbb = (f(z, zb+h) - 2*f(z, zb) + f(z, zb-h))/h**2
    Dz = z**2*(1-z)*fzz - z**2*fz; Db = zb**2*(1-zb)*fbb - zb**2*fb
    lhs = 2*Dz + 2*Db + 2*z*zb/(z-zb)*((1-z)*fz - (1-zb)*fb)
    Cas = mp.mpf(Delta)*(mp.mpf(Delta)-3) + J*(J+1)
    print(f"  Delta={Delta} J={J} z=({z},{zb}):  (D g)/(C g) - 1 = {mp.nstr(lhs/(Cas*f(z, zb)) - 1, 5)}")

print("(3) normalisation")
for (Delta, J) in ((4.2, 0), (6.6, 3), (9.1, 6)):
    z, zb = mp.mpf('1e-9'), mp.mpf('1e-4')
    val = gz(Delta, J, z, zb, N=8)
    print(f"  Delta={Delta} J={J}:  g/(z^((D-J)/2) zb^((D+J)/2)) = {mp.nstr(val/(z**((mp.mpf(Delta)-J)/2)*zb**((mp.mpf(Delta)+J)/2)), 10)}")

print("(4) convergence in N for rho, rhobar -> -1")
Delta, J = 40.3, 12
for (rb, rs) in ((0.9, 0.8), (0.97, 0.9), (0.99, 0.97), (0.995, 0.6)):
    vals = []
    for N in (400, 800, 1600, 3200):
        blk = Block3d(Delta, J, N)
        vals.append(blk.h_abs(rb, rs)[0])
    print(f"  |rho|=({rb},{rs}): h_abs for N=400,800,1600,3200:", " ".join(f"{v:.10e}" for v in vals))
