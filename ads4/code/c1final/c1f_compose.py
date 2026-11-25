#!/usr/bin/env python3
"""C1_FINAL step A3 [DV]: end-to-end s-channel T-family term of the k = 2 functional at O(t4), with the EE (2,4) contraction.

By c1f_affine (a), at x0 = infinity the k = 2 weight f = v'-u' is affine in the retarded times with no u1*u3 term, so the
light-ray product in Psi^+ is a combination of the three structures  E3 E1,  E3 E1^(1),  E3^(1) E1  (coefficients: scalar functions
of n1, n3, x2, x4).  Its s-channel T-family term is  <Psi4| A3 P_T A1 |Psi2>  = int d^3P sum_h conj((A3 f4)_h) (A1 f2)_h, with the
normalised one-point operators (c1f_freefields):
      (E f)_h  = Hh f,        (E1 f)_h = -i Hh nb.d f - (i/2)(nb.d Hh) f + Mh f.
The (2,4) dual structure enters through the external helicity components f^(c)_h(P) = g(P) (e_c . eps_h^*(P)), summed over a
Cartesian tensor basis c with the identity (EE, c1f_affine (b)) or with an improper reflection R (control, helicity-flip sensitive).
Parity-even <TTT> at fixed c_T:  Hh = Hh_E + (t4/4) Hh_W,  Hh_E = (Hh_B + Hh_F)/2,  Hh_W = (Hh_B - Hh_F)/2  (t4_B = +4, t4_F = -4);
the same for Mh.  The O(t4) term is the cross term  X = I[E3, W1] + I[W3, E1]  (I bilinear in the operator sets of lines 3 and 1);
each of the two cross terms (W at the (12T) vertex, or at the (T34) vertex = the bulk energy-channel exchange) is tested separately.
Claim (premise F2 + EE flip-zero): X = 0 POINTWISE in P, for every structure, every n1, n3 and every external g2, g4.
Controls: I[E,E] != 0 (Einstein term, which carries the graviton pole), I[W,W] != 0 (the |g3|^2 term, flip x flip), and X != 0 under the reflected (flip-sensitive) contraction.
Output: cache/c1f_compose.txt
"""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c1f_freefields as ff

out = []
def P_(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

R = np.diag([1.0, 1.0, -1.0])   # improper Lorentz reflection y -> -y (control)
HEL = ff.HEL

def ops(P, phi):
    """returns dict part -> (Hh, Mh, dHh) for part in E (Einstein), W (t4 structure), B, F."""
    nb = ff.nbar(phi); d = 1e-5
    HB, MB = ff.HM('B', P, phi); HF, MF = ff.HM('F', P, phi)
    HBp, _ = ff.HM('B', P + d*nb, phi); HBm, _ = ff.HM('B', P - d*nb, phi)
    HFp, _ = ff.HM('F', P + d*nb, phi); HFm, _ = ff.HM('F', P - d*nb, phi)
    dB = (HBp - HBm)/(2*d); dF = (HFp - HFm)/(2*d)
    return {'E': ((HB + HF)/2, (MB + MF)/2, (dB + dF)/2), 'W': ((HB - HF)/2, (MB - MF)/2, (dB - dF)/2)}

def ext(P, g, c, side, refl):
    """helicity components of the external state for Cartesian component c=(a,b) at momentum P."""
    a, b = c
    vals = []
    for h in HEL:
        e = ff.eps(P, h).conj()
        if refl: e = R @ e @ R.T
        if side == 4: e = ff.ETA @ e @ ff.ETA      # lowered indices: sum_c f4^* f2 = g4^* g2 (eps_h4 . eps_h2^*)
        vals.append(g(P)*e[a, b])
    return np.array(vals)

def dext(P, g, c, side, refl, phi):
    nb = ff.nbar(phi); d = 1e-5
    return (ext(P + d*nb, g, c, side, refl) - ext(P - d*nb, g, c, side, refl))/(2*d)

def apply(kind, O, P, g, c, side, refl, phi):
    H, M, dH = O
    f = ext(P, g, c, side, refl)
    if kind == 'E': return H @ f
    df = dext(P, g, c, side, refl, phi)
    return -1j*(H @ df) - 0.5j*(dH @ f) + M @ f

def I(struct, O3, O1, P, g4, g2, phi3, phi1, refl):
    k3, k1 = struct
    tot = 0j
    for a in range(3):
        for b in range(3):
            v4 = apply(k3, O3, P, g4, (a, b), 4, False, phi3)
            v2 = apply(k1, O1, P, g2, (a, b), 2, refl, phi1)
            tot += np.vdot(v4, v2)
    return tot

rng = np.random.default_rng(11)
def make_g(seed):
    r = np.random.default_rng(seed); Q = np.array([2.5, 0.3, -0.2]) + 0.5*r.normal(size=3); a = r.normal(size=3)
    return lambda P: np.exp(-np.sum((P - Q)**2)/1.5 + 1j*(a @ P))

structs = [('E', 'E'), ('E', 'E1'), ('E1', 'E')]
okEE = True; okCtl = True; okEin = True
P_("struct | P | phi3 phi1 | |I[E,E]| | max(|X|,|I[E3,W1]|,|I[W3,E1]|)/|I[E,E]| (EE) | |X_R|/|I[E,E]| (control) | |X_EE| abs | |I[W,W]|/|I[E,E]| (EE)")
worst = 0.0; minctl = 1e9; minww = 1e9
for trial in range(8):
    P = np.array([2.4, 0.0, 0.0]) + np.array([0.6, 0.5, 0.5])*rng.normal(size=3)
    if -ff.dot(P, P) <= 0.3 or P[0] <= 0: continue
    phi3, phi1 = rng.uniform(-np.pi, np.pi, size=2)
    g4, g2 = make_g(100 + trial), make_g(200 + trial)
    O3 = ops(P, phi3); O1 = ops(P, phi1)
    for s in structs:
        IEE = I(s, O3['E'], O1['E'], P, g4, g2, phi3, phi1, False)
        Xa = I(s, O3['E'], O1['W'], P, g4, g2, phi3, phi1, False); Xb = I(s, O3['W'], O1['E'], P, g4, g2, phi3, phi1, False)
        X = Xa + Xb; sep = max(abs(Xa), abs(Xb))/abs(IEE)
        XR = I(s, O3['E'], O1['W'], P, g4, g2, phi3, phi1, True) + I(s, O3['W'], O1['E'], P, g4, g2, phi3, phi1, True)
        IWW = I(s, O3['W'], O1['W'], P, g4, g2, phi3, phi1, False); wwr = abs(IWW)/abs(IEE); minww = min(minww, wwr)
        rel = max(abs(X)/abs(IEE), sep); relR = abs(XR)/abs(IEE)
        worst = max(worst, rel); minctl = min(minctl, relR)
        P_(f"{s[0]:>2}.{s[1]:<2} | {np.round(P,3)} | {phi3:+.2f} {phi1:+.2f} | {abs(IEE):.4e} | {rel:.2e} | {relR:.3e} | {abs(X):.1e} | {wwr:.3e}")
        okEE &= rel < 1e-8; okCtl &= relR > 1e-3; okEin &= abs(IEE) > 0
P_(f"worst |X_EE|/|I[E,E]| = {worst:.2e};  smallest control ratio = {minctl:.3e}")
P_("CHECK O(t4) T-family term vanishes pointwise under EE:", okEE)
P_("CHECK control (reflected, flip-sensitive (2,4) contraction) is O(1):", okCtl)
P_("CHECK Einstein term non-zero:", okEin)
P_(f"smallest |I[W,W]|/|I[E,E]| under EE = {minww:.3e}  (the |g3|^2 = W x W term is flip x flip = preserved: EE keeps it, as it must)")
okEin &= minww > 1e-3
P_("ALL CHECKS PASS:", okEE and okCtl and okEin)
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'), exist_ok=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'c1f_compose.txt'), 'w').write("\n".join(out) + "\n")
