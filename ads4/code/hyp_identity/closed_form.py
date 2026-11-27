"""
Exact evaluation of V(T0) for free fields from the closed-form kernels of spectral_kernels.py:

  eta_Y(T) = (1/24) sum_species Int_0^oo dE rho(E) (-dn/dE) (1 + 4 E^2)      (R = 1)
     scalar : rho = E tanh(pi E)/(2 pi),  n = Bose
     Dirac  : rho = 2 * E coth(pi E)/(2 pi) (particles+antiparticles), n = Fermi
  s, c_H from the same densities of states.

V(T0) = 4 pi (d eta/dT)/(ds/dT) - 1 ; since ds/dT = c_H/T and c_H(T0) = 4 pi eta(T0),
V(T0) = dlog eta_Y/dlog T |_{T0} - 1.
Also W := dlog c_H/dlog T |_{T0} (third cumulant <HHH>, the Renyi S''_q datum).
"""
import mpmath as mp
mp.mp.dps = 40
T0 = 1/(2*mp.pi)
nB = lambda E, T: 1/mp.expm1(E/T)
nF = lambda E, T: 1/(mp.exp(E/T) + 1)
rB = lambda E: E*mp.tanh(mp.pi*E)/(2*mp.pi)
rF = lambda E: 2*E/mp.tanh(mp.pi*E)/(2*mp.pi)
CUT = [0, 1, 3, 10, 40]
def eta(T, n, r): return mp.quad(lambda E: r(E)*(-mp.diff(lambda e: n(e, T), E))*(1 + 4*E**2), CUT)/24
def cH(T, n, r):  return mp.quad(lambda E: r(E)*E*mp.diff(lambda t: n(E, t), T), CUT)
CT = {'scalar': mp.mpf(3)/(32*mp.pi**2), 'dirac': mp.mpf(3)/(16*mp.pi**2)}
res = {}
for k, (n, r) in {'scalar': (nB, rB), 'dirac': (nF, rF)}.items():
    e0 = eta(T0, n, r); c0 = cH(T0, n, r)
    g = T0*mp.diff(lambda T: eta(T, n, r), T0)/e0
    W = T0*mp.diff(lambda T: cH(T, n, r), T0)/c0
    res[k] = (g - 1, W)
    print(f"{k}: eta_Y(T0)/(pi^2C_T/48) = {mp.nstr(e0/(mp.pi**2*CT[k]/48), 30)}")
    print(f"   c_H(T0)/(pi^3C_T/12) = {mp.nstr(c0/(mp.pi**3*CT[k]/12), 30)}")
    print(f"   V(T0) = dlog eta/dlogT - 1 = {mp.nstr(g-1, 30)} = {mp.identify(g-1)}")
    print(f"   W = dlog c_H/dlogT = {mp.nstr(W, 30)} = {mp.identify(W)}")

# Einstein (hyperbolic AdS4 black hole): s ~ r^2, T = (3r - 1/r)/4pi, eta = s/4pi
r = mp.mpf(1)
T = lambda r: (3*r - 1/r)/(4*mp.pi); s = lambda r: r**2
c = lambda r: T(r)*mp.diff(s, r)/mp.diff(T, r)
dlogT = mp.diff(lambda r: mp.log(T(r)), 1)
VE = mp.diff(lambda r: mp.log(s(r)), 1)/dlogT - 1
WE = mp.diff(lambda r: mp.log(c(r)), 1)/dlogT
print(f"Einstein: V(T0) = {mp.nstr(VE,20)}, W = {mp.identify(WE)}")

# Affine-in-t4 tests (t4 = +4 scalar/boson structure, -4 fermion, 0 Einstein)
VF, WF = res['dirac']; VB, WB = res['scalar']
print("\nClaim V = -t4/30:  predicted scalar -2/15, fermion +2/15")
print(f"  fermion: {mp.identify(VF)}  (claim {'OK' if abs(VF-mp.mpf(2)/15)<1e-30 else 'FAILS'})")
print(f"  scalar : {mp.identify(VB)}  (claim {'OK' if abs(VB+mp.mpf(2)/15)<1e-30 else 'FAILS'};"
      f" offset from claim = {mp.identify(VB+mp.mpf(2)/15)})")
print(f"  line through (Einstein, fermion): V = {mp.identify((VF-VE)/(-4))} * t4  -> predicts boson {mp.identify(VE+4*(VF-VE)/(-4))}")
print(f"  scalar+fermion (if both were <TTT>-faithful) would force V(t4=0) = {mp.identify((VB+VF)/2)} != 0 = Einstein")
slopeW = (WF - WE)/(-4)
print(f"\nThird cumulant W: line through (Einstein, fermion) has slope {mp.identify(slopeW)},"
      f" predicts boson {mp.identify(WE+4*slopeW)}; thermal scalar gives {mp.identify(WB)}"
      f" (offset {mp.identify(WB-WE-4*slopeW)})")
