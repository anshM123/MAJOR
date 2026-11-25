#!/usr/bin/env python3
"""C1_FINAL step A2 [DV]: the one-point matrices of E(n) and of the first moment E1(n) = int du u T(u,n) on the T family,
for the free boson (t4 = +4) and the free Dirac fermion (t4 = -4), which span the parity-even <TTT> at fixed c_T.

Checks (all at several P, rest frame and boosted, several directions phi):
 (0) sanity: slash(k) u(k) = 0; P^m t_mn = 0 (conservation) for both theories.
 (1) HM calibration: Hh_B, Hh_F have equal diagonals; flips opposite; |flip| = |diag|  (=> d=3 one-point matrix 1 +- flip,
     |t4| = 4, t4_B = -t4_F).  Einstein (t4 = 0) part (B+F)/2 is diagonal; t4 part (B-F)/2 is pure flip.
 (2) NEW: Mh (the first-moment matrix beyond the derivative term) has ZERO flip for both theories, and equal diagonals.
     Hence E1|_T = -i Hh nb.d - (i/2)(nb.dHh) + Mh with Mh_Einstein diagonal and Mh_t4 = 0  => premise F2's one-point pattern.
 (3) Ward identities: int dphi nb^mu Hh_hh = P^mu (momentum flux),  int dphi Mh_hh = 0 (dilatation flux with scale-invariant frames),
     int dphi nb^mu Hh_flip = 0.
 (4) the mechanism: d_w psi_h = (i Gamma_h + 1/w) psi_h with Gamma_{-h} = -Gamma_h (no admixture), both theories.
 (5) closed form M = (1/2){Gamma, H}: E1 acts on the T family as (1/2){E, -i nb.D}, D = frame-covariant derivative.
Output: cache/c1f_onepoint.txt
"""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c1f_freefields as ff

out = []
def P_(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

rng = np.random.default_rng(7)
Ps = [np.array([2.0, 0.0, 0.0]), np.array([3.0, 1.1, -0.7]), np.array([1.7, -0.9, 0.95]), np.array([5.0, 0.3, 4.1])]
phis = [0.0, 0.37, 1.9, -2.6, 3.0]
ok = True

# (0) sanity
k = np.array([1.3, 1.3*np.cos(0.8), 1.3*np.sin(0.8)])
P_("(0) |slash(k)u(k)| =", np.abs(ff.slash(k) @ ff.spinor(k)).max())
k1 = np.array([1.0, np.cos(0.3), np.sin(0.3)])*0.9; k2 = np.array([1.0, np.cos(2.2), np.sin(2.2)])*1.4
tB = np.outer(ff.lower(k1), ff.lower(k2)) + np.outer(ff.lower(k2), ff.lower(k1)) - ff.ETA*ff.dot(k1, k2)
Pt = k1 + k2
ub, v = ff.ubar(ff.spinor(k1)), ff.spinor(k2)
tF = np.array([[0.5*((ub @ (ff.ETA[m, m]*ff.GAM[m]) @ v)*ff.lower(k1 - k2)[n] + (ub @ (ff.ETA[n, n]*ff.GAM[n]) @ v)*ff.lower(k1 - k2)[m])
                for n in range(3)] for m in range(3)])
P_("(0) conservation |P^m t_mn|: boson", np.abs(Pt @ tB).max(), " fermion", np.abs(Pt @ tF).max(),
   "; fermion trace", abs(np.trace(np.linalg.inv(ff.ETA) @ tF)))

def offdiag(A): return np.array([A[0, 1], A[1, 0]])
P_("\n(1),(2): normalised matrices.  columns: P, phi | Hdiag_B, Hdiag_F | Hflip_B(+-), Hflip_F(+-) | |flip|/diag | Mflip_B, Mflip_F | Mdiag_B - Mdiag_F")
worst = dict(Hd=0, Hf=0, sat=0, Mf=0, Md=0, nrm=0)
for P in Ps:
    nB = ff.norm_density('B', P); nB2 = None
    for phi in phis:
        HB, MB = ff.HM('B', P, phi); HF, MF = ff.HM('F', P, phi)
        dHd = np.abs(np.diag(HB) - np.diag(HF)).max()/abs(HB[0, 0])
        dHf = np.abs(offdiag(HB) + offdiag(HF)).max()/abs(HB[0, 0])
        sat = abs(abs(HB[0, 1])/HB[0, 0].real - 1)
        mfl = max(np.abs(offdiag(MB)).max(), np.abs(offdiag(MF)).max())/abs(HB[0, 0])
        mdd = np.abs(np.diag(MB) - np.diag(MF)).max()/abs(HB[0, 0])
        worst['Hd'] = max(worst['Hd'], dHd); worst['Hf'] = max(worst['Hf'], dHf); worst['sat'] = max(worst['sat'], sat)
        worst['Mf'] = max(worst['Mf'], mfl); worst['Md'] = max(worst['Md'], mdd)
        P_(f"P={P}, phi={phi:+.2f} | {HB[0,0].real:.6e} {HF[0,0].real:.6e} | {HB[0,1]:.4e} {HF[0,1]:.4e} | "
           f"{abs(HB[0,1])/HB[0,0].real:.10f} | {abs(MB[0,1]):.1e} {abs(MF[0,1]):.1e} | {mdd:.1e}   (Mdiag_B = {MB[0,0].real:+.4e}, {MB[1,1].real:+.4e})")
P_("worst rel: Hdiag B-F", worst['Hd'], "; Hflip B+F", worst['Hf'], "; |flip|/diag - 1", worst['sat'],
   "; Mflip (B,F)", worst['Mf'], "; Mdiag B-F", worst['Md'])
ok &= worst['Hd'] < 1e-8 and worst['Hf'] < 1e-8 and worst['sat'] < 1e-8 and worst['Mf'] < 1e-7 and worst['Md'] < 1e-7
# helicity-independence of the norm density (h = -2 vs +2)
for th in ('B', 'F'):
    P = Ps[1]; M = np.sqrt(-ff.dot(P, P)); L = ff.boost(P)
    a = [abs(ff.PSI[th](L @ (M/2*np.array([1, np.cos(t), np.sin(t)])), L @ (M/2*np.array([1, -np.cos(t), -np.sin(t)])), h))**2
         for t in np.linspace(0, 2*np.pi, 400, endpoint=False) for h in (+2,)]
    b = [abs(ff.PSI[th](L @ (M/2*np.array([1, np.cos(t), np.sin(t)])), L @ (M/2*np.array([1, -np.cos(t), -np.sin(t)])), h))**2
         for t in np.linspace(0, 2*np.pi, 400, endpoint=False) for h in (-2,)]
    P_(f"norm density h=+2 vs h=-2 ({th}): rel diff {abs(np.mean(a)-np.mean(b))/np.mean(a):.1e}")

# (3) Ward identities
P_("\n(3) Ward: int dphi nb^mu Hh_hh vs P^mu; int dphi Mh_hh; int dphi nb^mu Hh_flip")
Nphi = 2000; grid = (np.arange(Nphi) + 0.5)*2*np.pi/Nphi
for th in ('B', 'F'):
    for P in Ps[1:3]:
        accH = np.zeros((3, 2, 2), complex); accM = np.zeros((2, 2), complex)
        for phi in grid:
            H, M = ff.HM(th, P, phi); nb = ff.nbar(phi)
            accH += nb[:, None, None]*H; accM += M
        accH *= 2*np.pi/Nphi; accM *= 2*np.pi/Nphi
        e1 = np.abs(accH[:, 0, 0] - P).max()/np.abs(P).max(); e2 = np.abs(accH[:, 1, 1] - P).max()/np.abs(P).max()
        e3 = np.abs(accH[:, 0, 1]).max()/np.abs(P).max(); e4 = np.abs(np.diag(accM)).max()/np.abs(P).max()
        P_(f"  {th} P={P}: momentum flux rel err (++,--) {e1:.1e} {e2:.1e}; flip moment {e3:.1e}; int Mh_diag {e4:.1e}")
        ok &= max(e1, e2, e3, e4) < 1e-6

# (4) mechanism: d_w log psi_h at fixed k2
P_("\n(4) d_w psi_h / psi_h = i Gamma_h + 1/w  (real part * w, Gamma_+ + Gamma_-)")
for th in ('B', 'F'):
    for P in Ps[1:]:
        for phi in phis[:3]:
            Hx, Mx, v, dv, w = ff._HM_one(th, P, phi, 1)
            r = dv/v
            P_(f"  {th} P={P} phi={phi:+.2f}: w*Re = {w*r[0].real:.10f} {w*r[1].real:.10f};  Im(+)+Im(-) = {r[0].imag + r[1].imag:.2e}")
            ok &= abs(w*r[0].real - 1) < 1e-7 and abs(w*r[1].real - 1) < 1e-7 and abs(r[0].imag + r[1].imag)*w < 1e-7
# (5) closed form: M = (1/2){Gamma, H}, Gamma = diag(Im d_w log psi_h)  =>  E1|_T = (1/2){E|_T, -i nb.D}, D = d + i Gamma
P_("\n(5) M = (1/2){Gamma,H} (time-of-arrival symmetrisation): max |M - (1/2){Gamma,H}| / |H_diag|")
w5 = 0.0
for th in ('B', 'F'):
    for P in Ps:
        for phi in phis:
            parts = [ff._HM_one(th, P, phi, 1)] + ([ff._HM_one(th, P, phi, 2)] if th == 'F' else [])
            for (Hx, Mx, v, dv, w) in parts:
                G = np.diag((dv/v).imag)
                w5 = max(w5, np.abs(Mx - 0.5*(G @ Hx + Hx @ G)).max()/abs(Hx[0, 0]))
P_(f"  worst = {w5:.2e}")
ok &= w5 < 1e-7
P_("\nALL CHECKS PASS:", ok)
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'), exist_ok=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'c1f_onepoint.txt'), 'w').write("\n".join(out) + "\n")
