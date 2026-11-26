# Route T' feasibility [DV-level scan]: termwise j-sum bound on the top leg (ideal psi2 weight), d^3-normalised, vs F_cert.
import sys, math, numpy as np, json
from fractions import Fraction
from an_common import nf02, LAM2
from scipy.special import loggamma
LAM = math.sqrt(LAM2)
p2, _ = nf02(); P2 = [float(x) for x in p2]
def psiabs(z):
    v = 0j
    for co in reversed(P2): v = v*z + co
    return abs(v)
rows = json.load(open('/home/user/MAJOR/ads4/code/lrtight/cache/lrt3_rcert.json'))['rows']
def fcert(beta):
    for r in rows:
        if float(Fraction(r['ba'])) <= beta <= float(Fraction(r['bb'])): return r['Fcert']
    return rows[-1]['Fcert']
def logB(nu):
    return (loggamma(2j*nu) - 2*loggamma(0.5 + 1j*nu)).real
def leg(L, Y, amin=998.0, nt=2000):
    d = math.exp(-L); c = 1 - d; beta = 2*math.atanh(math.sqrt(c))
    ts = np.linspace(0.25, LAM, nt); dt = ts[1]-ts[0]
    tot = 0.0
    for t in ts:
        nu = complex(t, Y)
        al = complex(2.5, 0) + 1j*nu
        A = 0.5 - 1j*nu; Bb = 1 - 2j*nu
        lc = 0.0; s = 0.0; lmax = -1e300; terms = []
        growth = 0.0
        for j in range(0, 20000):
            bmod = abs(al - j)
            # |S^(alpha-j)| envelope: 5 * prod (1 + |beta|^2/(3.9 a^2)) over the downward steps from sigma=1/4
            lt = lc + j*math.log(d) + (-al.real)*math.log(d)*(-1) * 0 + math.log(5) + growth
            terms.append(lt)
            if j > 50 and lt < max(terms) - 80: break
            growth += math.log1p(bmod**2/(3.9*amin**2))
            lc += math.log(abs((A + j)**2/((Bb + j)*(j + 1))))
        m = max(terms); s = m + math.log(sum(math.exp(x - m) for x in terms))
        # |d^{-alpha}| = d^{-5/2} e^{-Y L}
        lv = logB(nu) + s - 2.5*math.log(d) - Y*L + math.log(psiabs(nu/LAM))
        tot += math.exp(lv)*dt
    lead = 2*(2/LAM)*tot*d**3          # d^3-normalised, 2Re -> factor 2
    return lead, fcert(beta), beta
for L in [float(x) for x in sys.argv[1].split(',')]:
    out = []
    for Y in [float(x) for x in sys.argv[2].split(',')]:
        lg, fc, be = leg(L, Y)
        out.append(f'Y={Y}: {lg/fc:.2e}')
    print(f'L={L} beta={be:.3f} F_cert={fc:.3e}: leg/F_cert  ' + '  '.join(out), flush=True)
