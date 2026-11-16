"""
Stage 1 numerics.
(a) Scalar d = 3 blocks (exact, blocks3d.py) against the leading large-Delta form at fixed J
        g_{Delta,J}(r, eta) ~ (4r)^Delta * N_J P_J(eta) / [ (1-r^2)^{1/2} ((1+r^2)^2 - 4 r^2 eta^2)^{1/2} ],  N_J = J!/(1/2)_J
    [Kos-Poland-Simmons-Duffin, arXiv:1406.4858, h_infinity for equal external dimensions; quoted from memory, tested here],
    at fixed r < 1 (also for J ~ Delta), and its failure in the region relevant for the Regge moments,
    1 - |rho_s| ~ 1/Delta (u-channel Regge limit of the s-channel block), where the plane-wave form of CHMRS (3.13) takes over.
(b) Conserved <J phi O_l>, <T phi O_l> (threept.py) in the collinear frame x2 < x1 < x5 on the time axis, z5 = (1, cos al, sin al):
    SO(2) (rotations about the axis) charge content of the vertex at fixed l, Delta -> infinity.
    Prediction (Chang-Landau-Simmons-Duffin, arXiv:2311.04271, section 4.1: conservation -> v_{1,23}.<...> = 0): only charges
    q_1 = +-l_1 survive, i.e. the q-basis label of a conserved operator becomes its helicity.
"""
import numpy as np, mpmath as mp
from blocks3d import Block3d
from threept import *
from s2_three_point_heavy import conserved_coeffs_J

print("(a) exact / asymptotic  for rho, rhobar > 0")
for (rho, rhob) in ((0.35, 0.2), (0.6, 0.45)):
    r = np.sqrt(rho*rhob); eta = (rho+rhob)/(2*r)
    for J in (0, 2, 6):
        row = []
        for Delta in (10.3, 40.3, 160.3, 640.3):
            blk = Block3d(Delta, J, 400)
            ex = blk.h_abs(max(rho, rhob), min(rho, rhob))[0]
            NJ = float(mp.factorial(J)/mp.rf(mp.mpf(1)/2, J))
            asym = NJ*float(mp.legendre(J, eta))/np.sqrt((1-r*r)*((1+r*r)**2 - 4*r*r*eta*eta))
            row.append(ex/asym)
        print(f"   rho=({rho},{rhob}) J={J}:  Delta=10.3,40.3,160.3,640.3 -> " + "  ".join(f"{x:.6f}" for x in row))
    row = []
    for Delta in (10.3, 40.3, 160.3):
        J = int(round(0.4*Delta))
        blk = Block3d(Delta, J, 400)
        ex = blk.h_abs(max(rho, rhob), min(rho, rhob))[0]
        NJ = float(mp.factorial(J)/mp.rf(mp.mpf(1)/2, J))
        asym = NJ*float(mp.legendre(J, eta))/np.sqrt((1-r*r)*((1+r*r)**2 - 4*r*r*eta*eta))
        row.append(ex/asym)
    print(f"   rho=({rho},{rhob}) J=0.4*Delta (Regge corner): " + "  ".join(f"{x:.4f}" for x in row) + "   (still -> 1: at fixed r < 1 the form holds also for J ~ Delta)")

print("    region 1-|rho| = c/Delta, 1-|rhobar| = cbar/Delta (both negative; (c, cbar) = (6, 3)), J = 0.3 Delta:")
row = []
for Delta in (20.3, 40.3, 80.3, 160.3):
    J = int(round(0.3*Delta)); rb, rs = 1 - 3/Delta, 1 - 6/Delta
    blk = Block3d(Delta, J, int(60*Delta)//2*2)
    ex = blk.h_abs(rb, rs)[0]
    r = np.sqrt(rb*rs); eta = (rb+rs)/(2*r)
    NJ = float(mp.factorial(J)/mp.rf(mp.mpf(1)/2, J))
    asym = NJ*float(mp.legendre(J, eta))/np.sqrt((1-r*r)*((1+r*r)**2 - 4*r*r*eta*eta))
    row.append(ex/asym)
print("   exact/asymptotic for Delta = 20.3, 40.3, 80.3, 160.3: " + "  ".join(f"{x:.4f}" for x in row) + "   <- does not tend to 1")

print("(b) collinear frame: charge content of the conserved vertices, l = 6 fixed")
mp.mp.dps = 40
x1 = [mp.mpf(0), 0, 0]; x2 = [mp.mpf('-1.7'), 0, 0]; x5 = [mp.mpf('0.8'), 0, 0]
D2 = mp.mpf('1.7'); l = 6
for l1 in (1, 2):
    for Delta in (mp.mpf(20), mp.mpf(200), mp.mpf(2000), mp.mpf(20000)):
        dims, spins = (mp.mpf(l1+1), D2, Delta), (l1, 0, l)
        monos, sv, V = conserved_coeffs_J(dims, spins, 'even')
        cvec = [V[len(monos)-1, a] for a in range(len(monos))]
        al = mp.mpf('0.37')
        z5 = [mp.mpf(1), mp.cos(al), mp.sin(al)]
        poly = lambda z: sum(cvec[a]*poly_even(monos[a], blocks(x1, z, x2, E3[0], x5, z5), spins) for a in range(len(monos)))
        s = tensor_at_1(poly, l1)
        if l1 == 1:
            t = [ETA[i]*s[i] for i in range(3)]
            long, trans = abs(t[0]), mp.sqrt(abs(t[1])**2 + abs(t[2])**2)
            print(f"   <J phi O>  Delta={mp.nstr(Delta, 6):>8s}:  |t.e_t| / |t_transverse| = {mp.nstr(long/trans, 6)}")
        else:
            # traceless part of the quadratic form; components with at least one time index vs purely spatial traceless part
            tr = -s[0][0] + s[1][1] + s[2][2]
            tl = [[s[i][j] - (ETA[i] if i == j else 0)*tr/3 for j in range(3)] for i in range(3)]
            time_part = mp.sqrt(abs(tl[0][0])**2 + 2*abs(tl[0][1])**2 + 2*abs(tl[0][2])**2)
            sp_tr = tl[1][1] + tl[2][2]
            spat = mp.sqrt(abs(tl[1][1] - sp_tr/2)**2 + abs(tl[2][2] - sp_tr/2)**2 + 2*abs(tl[1][2])**2)
            print(f"   <T phi O>  Delta={mp.nstr(Delta, 6):>8s}:  |components with a time index (q_1 = 0, +-1)| / |q_1 = +-2 part| = {mp.nstr(time_part/spat, 6)}")
