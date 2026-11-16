"""Casimir-equation residual of a block evaluated from its jets E_{al,be} = r^{-Delta} th_s^al th_b^be [r^Delta Gt]
(th = rho d/drho), in the Dolan-Osborn form  D_z + D_zb + (d-2) z zb/(z-zb) [(1-z) d_z - (1-zb) d_zb]  = C/2,
D_z = z^2(1-z) d_z^2 - (1+a+b) z^2 d_z - a b z,  (a, b) = (c1, c2) of the reflected block, d = 3.
With z d_z = w th, w = (1+rho)/(1-rho), th w = 2 rho/(1-rho)^2."""
import numpy as np

def casimir_residual(E, rs, rb, Delta, J, a, b, d=3):
    """E[al, be, p] (al + be <= 2 used).  Returns |residual| / (largest individual term) per point: near rho_b -> 1 the
    individual terms are ~ (1-rho_b)^{-2} times g and cancel, so this is the relevant measure of the relative jet accuracy."""
    F, Fs, Fb, Fss, Fbb = E[0, 0], E[1, 0], E[0, 1], E[2, 0], E[0, 2]
    def one(rho, F1, F2):
        w = (1 + rho)/(1 - rho); thw = 2*rho/(1 - rho)**2
        z = 4*rho/(1 + rho)**2
        zd = w*F1                                   # z d_z g
        zd2 = w*w*F2 + w*thw*F1                     # (z d_z)^2 g
        Dz = (1 - z)*(zd2 - zd) - (1 + a + b)*z*zd - a*b*z*F
        sc = np.maximum.reduce([np.abs((1 - z)*zd2), np.abs((1 - z)*zd), np.abs((1 + a + b)*z*zd)])
        return Dz, z, zd, sc
    Ds, zs, zds, scs = one(rs, Fs, Fss)
    Db, zb, zdb, scb = one(rb, Fb, Fbb)
    m1 = (d - 2)/(zs - zb)*zb*(1 - zs)*zds; m2 = (d - 2)/(zs - zb)*zs*(1 - zb)*zdb
    C = Delta*(Delta - 3) + J*(J + 1)
    scale = np.maximum.reduce([scs, scb, np.abs(m1), np.abs(m2), np.abs(C/2*F)])
    return np.abs(Ds + Db + m1 - m2 - C/2*F)/scale
