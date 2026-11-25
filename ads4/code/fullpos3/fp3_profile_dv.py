"""
fullpos3/fp3_profile_dv.py (FULLPOS3.md §1) [DV]: the h-expansion profiles at c with the target weights psi2, psi3.
  I_j(c) = c^j G^{(j)}(c) for G = G_A, G_B, G_C, G_D (definitions in FULLPOS3.md §1.1), j <= JMAX;  composite Gauss-Legendre in p.
  P_1(theta) G = -(d/c)[I_3/3 + (b+1) I_2];  P_2(theta)G from P_2(k) = a1^2/2 - a2 (falling-factorial basis).
  Abs_j = 2 int |psi2| |c^j F_3^{(j)}(c)| dp  (absolute scale).
Usage: fp3_profile_dv.py c1,c2,...  [npanels] [nodes]     Output appended to cache/fp3_profile_dv.txt (json lines)
"""
import sys, os, json, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from flint import arb, acb, ctx
import fp3_common as F

JMAX = 9
p2, p3 = F.nf02()
Lam = math.sqrt(F.LAM2)


def profile(c, npan=64, nnod=64):
    t0 = time.time()
    xg, wg = np.polynomial.legendre.leggauss(nnod)
    ps = []; ws = []
    for i in range(npan):
        a, b = i/npan, (i + 1)/npan
        ps += list(a + (b - a)*(xg + 1)/2); ws += list((b - a)/2*wg)
    fac = [math.factorial(j) for j in range(JMAX + 1)]
    GA = np.zeros(JMAX + 1); GB = np.zeros(JMAX + 1); GC = np.zeros(JMAX + 1); GD = np.zeros(JMAX + 1); AbsA = np.zeros(JMAX + 1)
    cA = arb(c)
    for p, w in zip(ps, ws):
        nu = Lam*p
        ctx.prec = int(2.4*nu + 300)
        t3 = F.taylorF(3, nu, cA, JMAX); t4 = F.taylorF(4, nu, cA, JMAX)
        t3s = F.taylorF(3, acb(nu, 1), cA, JMAX)
        ell = acb(-0.5, nu)
        y2 = F.polyval_fr(p2, p); y3 = F.polyval_fr(p3, p)
        for j in range(JMAX + 1):
            cj = c**j*fac[j]
            v3 = float(t3[j].real.mid())*cj; v4 = float(t4[j].real.mid())*cj
            v3s = float((ell/2*t3s[j]).real.mid())*cj
            GA[j] += w*2*y2*v3
            AbsA[j] += w*2*abs(y2*v3)
            GB[j] += w*2*y2*(nu*nu + 6.25)/2*v4
            GC[j] += w*2*y2*(2.5*v4 - 2.75*v3 - v3s)
            GD[j] += w*F.LAM2*y3*v4
    return dict(c=c, GA=GA.tolist(), GB=GB.tolist(), GC=GC.tolist(), GD=GD.tolist(), AbsA=AbsA.tolist(), npan=npan, nnod=nnod,
                secs=time.time() - t0)


if __name__ == '__main__':
    cs = [float(x) for x in sys.argv[1].split(',')]
    npan = int(sys.argv[2]) if len(sys.argv) > 2 else 64; nnod = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    for c in cs:
        r = profile(c, npan, nnod)
        print(json.dumps(r), flush=True)
        open(os.path.join(F.CACHE, 'fp3_profile_dv.txt'), 'a').write(json.dumps(r) + '\n')
