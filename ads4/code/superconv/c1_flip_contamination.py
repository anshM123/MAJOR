"""c1: size of the helicity-flip contamination of the even moments of Omega^U (M02) with the project kernel D1 (x) D1 (D1^2 per pair),
relative to the stripped component norm N_s, on near-gap states; and with the flip-free kernel S_ff = D1^2 + D1 D2 (exactly zero).

Model [M]: leading two-pair Regge formula (SPINNING_4PT 4.2) with exact l=2 M-tables at a = k + 5:
  eps_k m^{2k} C^+_k  =  2 eps_k { P_{-1/2+i nu}(cosh b)(|c_s|^2+|c_o|^2)
                                  + rho(nu) K_{4,0} X + rho(-nu) K_{0,4} X  + rho(nu) rho(-nu) K_{4,4} (|c_s|^2+|c_o|^2) }
  rho = Mflip/Mpres of D1^2 (F-normalised) = (2nu-i)(2nu-3i)(2nu-5i)(2nu-7i)/N1(nu,a);  X = interference (|X| <= (|c_s|^2+|c_o|^2)/2).
Bound per unit |c|^2 used here: |rho|(|K40|+|K04|)/2 + |rho|^2 |K44|  (O(1) factors of the coupling bookkeeping not fixed: estimate).
Normal form as omegau/common.py: part_k = u^{k-2} int phi_k(p) [eps_k m^{2k} C_k] dp, u = Lam^2/m^2, nu = Lam p.
Caveat: the Regge formula needs nu << m; here nu = Lam p reaches m/mfl -> used as a model up to x = p/mfl.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from common import load, state, KS, EPS, SIGMA
from rcommon import peval, flat_parts
from numpy.polynomial.legendre import leggauss

def N1(nu, a):
    return 64*a**4 - 768*a**3 + 64*a**2*nu**2 + 3280*a**2 - 320*a*nu**2 - 5840*a + 16*nu**4 + 488*nu**2 + 3705
def rho(nu, a):
    return (2*nu - 1j)*(2*nu - 3j)*(2*nu - 5j)*(2*nu - 7j)/N1(nu, a)

TH, WT = leggauss(2000); TH = np.pi*(TH + 1); WT = np.pi*WT
def K(nu, b, muL, muR):
    ch, sh = np.cosh(b), np.sinh(b)
    base = ch - sh*np.cos(TH)
    eR = (sh - ch*np.cos(TH) - 1j*np.sin(TH))/base
    f = base**(-0.5 + 1j*nu)*np.exp(1j*muL*TH)*eR**muR
    return np.sum(WT*f)/(2*np.pi)

xg, wg = leggauss(48); pn = (xg + 1)/2; pw = wg/2
if __name__ == "__main__":
    F, PH = load('M02')
    print("c1: M02 Omega^U, even k = 2, 4, 6 with kernel D1^2 (x) D1^2: flip contamination bound / N_s  (leading Regge model)")
    print("    columns: Dg, h, m/Lam, J, beta, N_s, |pres action|, flip bound (k=2,4,6 summed), ratio flip/N_s, [p <= 0.6 mfl only]")
    for Dg in (100, 1000):
        Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        for h in ('++', '+-'):
            Jlist = (0, 2, 4, 6) if h == '++' else (4, 6, 8)
            for mfl in (1.2, 1.4, 1.6, 2.0):
                u = 1/mfl**2
                Js, B2, B3, FW = flat_parts(F, h, u, 12)
                for J in Jlist:
                    st = state(Dg, Lam, mfl, J)
                    if st is None:
                        continue
                    b = st['beta']
                    Ns = abs(B2[J - Js[0]]) + abs(B3[J - Js[0]]) + abs(FW[J - Js[0]])
                    flip = 0.0; flip_r = 0.0; pres = 0.0
                    for i, p in enumerate(pn):
                        nu = Lam*p
                        P = K(nu, b, 0, 0).real
                        k40 = abs(K(nu, b, 4, 0)); k04 = abs(K(nu, b, 0, 4)); k44 = abs(K(nu, b, 4, 4))
                        for k in (2, 4, 6):
                            a = k + 5
                            r = abs(rho(nu, a))
                            wk = pw[i]*u**(k - 2)*abs(peval(PH[k], p))*2*EPS[k]
                            fb = wk*(r*(k40 + k04)/2 + r*r*k44)
                            flip += fb
                            if p <= 0.6*mfl: flip_r += fb
                            pres += pw[i]*u**(k - 2)*peval(PH[k], p)*2*EPS[k]*P
                    print(f"  {Dg:5d} {h} {mfl:4.2f} J={J:2d} beta={b:.4f}  N_s={Ns:.4g}  |pres|={abs(pres):.4g}  flip={flip:.4g}  flip/N_s={flip/Ns:.3g}  [x<=0.6: {flip_r/Ns:.3g}]", flush=True)
    # forward term: rho(0) at a = 9 (k = 4), d_{nu^2} of the bound at nu = 0 relative to d_{nu^2} P
    a = 9
    print(f"forward k=4: |rho(0)| = {abs(rho(0.0, a)):.4g} (D1^2), 0 (S_ff)")
    print("S_ff = D1^2 + D1 D2: flip response = (+1 - 1)/(4 dd dd) = 0 identically (Lemma L) -> all flip terms vanish; preserved response (N1-N2)/(64 dd dd) > 0.")
