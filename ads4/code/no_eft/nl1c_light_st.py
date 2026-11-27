r"""no_eft/nl1c_light_st.py [NE; continuation of Eq. (fp-target) to light states is D-sketch]: the same ladder evaluation on
light SINGLE-trace spins J = 0, 2 (b = 0 ladder: ++ and the parity-odd structure), at every dimension from the unitarity bound to the
gap.  Side check only (the headline cubic bound assumes no light single trace); compare LS4/LS5 and the flat polynomials
F_D^(0) = 11.317027 + 20.771893 m^2, F_D^(2) of supplement Sec. 'Light single traces'.
Prints hhat normalised by Lam^2 * (its value convention), the J = 0 closed form, and the minimum over a fine Delta grid."""
import sys, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
import fp_poly as FP
from flint import arb, fmpq
import nl1_light_action as N
for Dg in (300, 1000):
    Z, Y = FP.load_moments(Dg, 300); lam = (Dg - 2)*(Dg - 1)
    A = FP.Action(Dg, Z, Y, FP.Q_exact(300, lam))
    c1 = A.Z[0]/3; c0 = 2*A.Y[0] + A.Lam2*A.Q[0]
    print(f'Delta_g={Dg}: J=0 closed form hhat = (Z_0/3) T + (2 Y_0 + Lam^2 Q_0) = {float(c1.mid()):.6f} T + {float(c0.mid()):.6e}'
          f'   -> per Lam^2: {float(c1.mid()):.6f} (T/Lam^2) + {float((c0/A.Lam2).mid()):.6f}')
    for J, D0 in ((0, 0.5), (2, 3.05)):
        mn = None
        D = D0
        while D <= Dg + J:
            if abs((D - 1)*(D - 2)) > 1e-9:
                h, T, r = N.hhat(A, J, 0, fmpq(int(round(D*1000)), 1000))
                v = float((h/A.Lam2).mid())
                if v != v: continue
                if mn is None or v < mn[0]:
                    mn = (v, D)
            D += 0.05 if D < 20 else 1.0
        print(f'   J={J}: min over Delta in [{D0}, {Dg + J}] of hhat/Lam^2 = {mn[0]:.6e} at Delta = {mn[1]:.2f}')
