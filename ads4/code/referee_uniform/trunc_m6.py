"""Referee: 1/m^6 truncation of the fixed-J remainder (r1/m^2 + r2/m^4) vs the order-8 fixed-eta series, near the gap."""
import sys, numpy as np
U = '/home/user/MAJOR/ads4/code/uniform'; sys.path.insert(0, U)
import remainder as RM, skelfast
for Dg in (20, 50):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    for mfl in (1.0, 1.1, 1.3):
        for J in (0, 2, 4):
            m = mfl*Lam; Bq = 2*J + 1.0
            al = (-Bq + np.sqrt(Bq*Bq + 4*m*m))/2
            if al < Dg - 2 - 1e-9: continue
            ab = al + 2*J + 1; m2 = al*ab; eta = (al**2 + ab**2)/(2*al*ab); Delta = al + J + 2
            x = np.array([0.2, 0.35, 0.5]); nus = x*m
            Pd = skelfast.Pders_vec(nus, eta)
            for k in (2, 3):
                r1, e1, c1 = RM.residual_fixedJ(k, x, J, m2)
                r8, e8, c8 = RM.residual_seriesC(k, nus, Delta, J, m2, eta, Pd)
                print(f"Dg={Dg} m/L={mfl} J={J} k={k}: rho_fixedJ={np.round(r1,4)} (err {np.round(e1,4)}, ctrl {c1})  rho_order8={np.round(r8,4)} (err {np.round(e8,4)})")
