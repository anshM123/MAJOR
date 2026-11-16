"""
Stage 3: absence of odd powers of 1/m in the (1<->3)-symmetric part of the spinning vertex tensor.

Claim (REGGE_MOMENT_TTTT.md, section 6.1): let F(x, y; Delta, l) be the left vertex (operators 1, 2, O) of an analytic family of
conserved structures, in the symmetric frame.  The right vertex (O^S, 3, 4) is F(-x, -y; 3 - Delta, -1 - l).  Hence
     tau_L (x) tau_R (x, y; Delta, l)  ->  transpose   under (x, y; alpha, alphabar) -> (-x, -y; -alpha, -alphabar),
and, with x = X/m, y = Y/m, m^2 = alpha*alphabar, the symmetric part of the normalised vertex tensor is even in m at fixed
(X, Y, eta_AdS): its corrections start at 1/m^2, provided frames are built from (alpha, alphabar):
     pbarhat = (cosh(beta/2), -sinh(beta/2), 0),  tanh(beta/2) = (alphabar - alpha)/(alphabar + alpha) = (2l+1)/(2 Delta - 3).
Test: <J phi O_l> x <O^S J phi>, conserved parity-even structures; W := sym(t_L (x) t_R)/(|t_L||t_R|) in the frame
(pbarhat, n_par, e_perp).  Non-transverse symmetric components W_pn, W_pe, W_pp and the mixing W_ne must scale like 1/m^2
(they are O(1/m) for each vertex separately, and O(1/m) in the frame built from (Delta, -l)).
"""
import mpmath as mp
from threept import *
from s2_three_point_heavy import conserved_coeffs_J, EY, E0, N5

mp.mp.dps = 50
X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
D2 = mp.mpf('1.7')

def vertices(l1, Delta, l):
    al, alb = Delta - l - 2, Delta + l - 1
    m = mp.sqrt(al*alb)
    x = [c/m for c in X]; y = [c/m for c in Y]
    xx = dot(x, x)
    x1, x2, x3, x4 = [-c for c in y], [c/xx for c in x], y, [-c/xx for c in x]
    out = []
    for (xa, xb, D5, l5) in ((x1, x2, Delta, l), (x3, x4, 3 - Delta, -1 - l)):
        dims, spins = (mp.mpf(l1+1), D2, D5), (l1, 0, l5)
        monos, sv, V = conserved_coeffs_J(dims, spins, 'even')
        cvec = [V[len(monos)-1, a] for a in range(len(monos))]
        poly = lambda z: sum(cvec[a]*poly_even(monos[a], blocks(xa, z, xb, E3[0], E0, N5), spins) for a in range(len(monos)))
        out.append(tensor_at_1(poly, l1))
    th = mp.atanh((alb - al)/(alb + al))
    frames = {'Weyl-covariant frame (alpha, alphabar)': ([mp.cosh(th), -mp.sinh(th), 0], [mp.sinh(th), -mp.cosh(th), 0]),
              'naive frame (Delta, -l)': ([Delta/mp.sqrt(Delta**2-l**2), -l/mp.sqrt(Delta**2-l**2), 0], [l/mp.sqrt(Delta**2-l**2), -Delta/mp.sqrt(Delta**2-l**2), 0])}
    return m, out, frames

print("# <J phi O_l> x <O^S J phi>:  l/Delta ~ 0.4.  Symmetric part of t_L (x) t_R, components in (pbarhat, n_par, e_perp), normalised to W_nn")
for Delta in (mp.mpf('50.3'), mp.mpf('100.3'), mp.mpf('200.3'), mp.mpf('400.3'), mp.mpf('800.3')):
    l = mp.floor(mp.mpf(2)/5*Delta)
    m, (sL, sR), frames = vertices(1, Delta, l)
    tL = [ETA[i]*sL[i] for i in range(3)]; tR = [ETA[i]*sR[i] for i in range(3)]
    for name, (pb, npar) in frames.items():
        cL = [dot(tL, pb), dot(tL, npar), dot(tL, EY)]; cR = [dot(tR, pb), dot(tR, npar), dot(tR, EY)]
        W = [[(cL[i]*cR[j] + cL[j]*cR[i])/2 for j in range(3)] for i in range(3)]
        nn = W[1][1]
        print(f"  m = {mp.nstr(m, 6):>8s}  {name:40s}: m*W_pn = {mp.nstr(m*W[0][1]/nn, 6):>12s}  m*W_ne = {mp.nstr(m*W[1][2]/nn, 6):>12s}  m*W_pe = {mp.nstr(m*W[0][2]/nn, 6):>12s}   m^2*W_pn = {mp.nstr(m*m*W[0][1]/nn, 6):>10s}  m^2*W_ne = {mp.nstr(m*m*W[1][2]/nn, 6):>10s}")
