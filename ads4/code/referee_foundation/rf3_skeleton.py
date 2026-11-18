"""
Referee check (a): independent evaluation of the exact scalar-skeleton moment m^{2k} C_{k,nu} for a SCALAR (J = 0) heavy exchange,
compared with the exactskel caches (e1_<Dg>_Df3.pkl, which come from skelexact.py (Dg=100) / skelexact3.py (Dg=300,1000)).

What is independent of the authors' code:
 * Mack data for J = 0 from the Dolan--Osborn u-expansion of the scalar block (not mack.py):
     Q_n  propto  ((b)_n)^2 / (n! (Delta-1/2)_n),   b = 1 - Df + Delta/2   (identical external Df, d = 3).
 * CHMRS (CfromBMellin) + (Bkt action), with the T-integral done by Barnes' first lemma, which I re-derived:
     int ds/(2 pi i) Gamma(s)^2 Gamma(c-s) Gamma(cbar-s) e^{2su} = |Gamma(c)|^4 2F1(c,cbar;A;1-e^{-2u})/Gamma(A),  c = (A+i nu)/2, A = a-1/2.
   and 1/(z)_kap, 1/(z (z)_kap) as Mellin transforms on t in (0,1) (Beta integrals), so that
     C_k propto 2^{1-2k} sum_n Q_n (-1)^kap/(Df - tau/2 - n)_kap * int_0^oo dy e^{-s_n y} phi_{n,kap}(e^{-y}) 2F1(c,cbar;A;1-e^{-y})/Gamma(A),
     s_n = tau/2 + n,  phi = [c_n (1-t)^{kap-1} - G_kap(t)/2]/Gamma(kap),  c_n = 1/(2Df - tau - 2n),  G_kap(t) = int_t^1 (1-s)^{kap-1}/s ds.
   The y-integral is done by mpmath quadrature with mpmath's hyp2f1 (not the w-series / beta-function sums of skelexact3),
   and the n-sum directly to N plus a fitted power-law tail.
 * One global normalisation constant (common to all k and nu) is fixed on the k = 2, smallest-nu node; everything else is a test.
Usage: python3 rf3_skeleton.py Dg mfl [node indices]   -> cache/rf3_<Dg>_<mfl>.txt
"""
import os, sys, pickle, time
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
ES = os.path.join(HERE, '..', 'exactskel', 'cache')
mp.mp.dps = 25

Dg = int(sys.argv[1]); mfl = float(sys.argv[2])
idx = [int(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0, 20, 32, 44, 54, 60, 63]
Df = mp.mpf(3)
db = pickle.load(open(os.path.join(ES, 'e1_%d_Df3.pkl' % Dg), 'rb'))
rec = db[(round(mfl, 6), 0)]
Delta = mp.mpf(float(rec['Delta'])); tau = Delta
m2 = (Delta - 2)*(Delta - 1)
nus = rec['nus']; ex = rec['ex']
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
pr(f"# rf3: Dg={Dg} m/Lam={mfl} J=0 Delta={float(Delta):.6f} Df=3; nodes {idx}; x = {[round(float(rec['x'][i]),4) for i in idx]}")

b = 1 - Df + Delta/2
def logQ(n):
    return 2*mp.log(mp.rf(b, n)) - mp.log(mp.factorial(n)) - mp.log(mp.rf(Delta - mp.mpf(0.5), n))

def G(kap, t, y):
    # int_t^1 (1-s)^{kap-1}/s ds = -log t + sum_{j>=1} C(kap-1,j)(-1)^j (1-t^j)/j
    s = y
    for j in range(1, kap):
        s += mp.binomial(kap - 1, j)*(-1)**j*(1 - t**j)/j
    return s

def In(n, kap, Fcache):
    sn = tau/2 + n
    cn = 1/(2*Df - tau - 2*n)
    def f(y):
        t = mp.e**(-y)
        phi = (cn*(1 - t)**(kap - 1) - G(kap, t, y)/2)/mp.gamma(kap)
        return mp.e**(-sn*y)*phi*Fcache(y)
    # scale: e^{-s_n y}; integrate on [0, 60/s_n] split, the rest is negligible (F grows at most like e^{pi nu/2})
    Y = mp.mpf(80)/sn
    return mp.quad(f, mp.linspace(0, Y, 9))

from scipy.special import gammaln, roots_genlaguerre
import numpy.polynomial.chebyshev as C
YL, WL = roots_genlaguerre(80, 0.0)

def Fgrid(k, nu, zmax=0.995, npts=400):
    """log of 2F1(c,cbar;A;z)/Gamma(A) on Chebyshev nodes in z, from mpmath, for barycentric interpolation"""
    a = 2*Df + k - 1; A = a - mp.mpf(0.5)
    c = (A + 1j*mp.mpf(nu))/2
    j = np.arange(npts); zc = 0.5*zmax*(1 - np.cos(np.pi*(j + 0.5)/npts))
    lv = np.array([float(mp.log(mp.re(mp.hyp2f1(c, mp.conj(c), A, mp.mpf(z))))) for z in zc]) - float(mp.loggamma(A))
    cf = C.chebfit(2*zc/zmax - 1, lv, npts - 1)
    return lambda z: C.chebval(2*np.asarray(z)/zmax - 1, cf), zmax

def Ck(k, nu, N=6000):
    kap = k//2
    lF, zmax = Fgrid(k, nu)
    tf = float(tau); Dff = float(Df); bb = float(b); Dl = float(Delta)
    n = np.arange(N, dtype=float)
    lQ = 2*(gammaln(bb + n) - gammaln(bb)) - gammaln(n + 1) - (gammaln(Dl - 0.5 + n) - gammaln(Dl - 0.5))
    # (-1)^kap/(Df - tau/2 - n)_kap : product of (Df - tau/2 - n + i), all negative for n large; sign handled explicitly
    P = np.ones(N)
    for i in range(kap):
        P = P*(Dff - tf/2 - n + i)
    pre = (-1)**kap/P
    sn = tf/2 + n
    cn = 1/(2*Dff - tf - 2*n)
    I = np.zeros(N)
    for jn in range(N):
        y = YL/sn[jn]
        t = np.exp(-y); z = 1 - t
        ok = z < zmax
        Gk = y.copy()
        for jj in range(1, kap):
            Gk += float(mp.binomial(kap - 1, jj))*(-1)**jj*(1 - t**jj)/jj
        phi = (cn[jn]*(1 - t)**(kap - 1) - Gk/2)/float(mp.gamma(kap))
        val = np.where(ok, phi*np.exp(np.where(ok, lF(np.minimum(z, zmax)), 0.0)), 0.0)
        I[jn] = np.sum(WL*val)/sn[jn]
    lt = lQ + np.log(np.abs(pre*I)); sg = np.sign(pre*I)
    mx = lt.max()
    terms = sg*np.exp(lt - mx)
    tot = terms.sum()
    ns = np.arange(N//2, N)
    p = -np.polyfit(np.log(ns), np.log(np.abs(terms[N//2:])), 1)[0]
    tail = terms[-1]*(N - 1)**p/(p - 1)/(N - 0.5)**(p - 1)
    lscale = mx + (1 - 2*k)*np.log(2) + k*np.log(float(m2))
    return mp.e**lscale*mp.mpf(tot), mp.e**lscale*mp.mpf(tail), p, tot

t0 = time.time()
res = {}
for k in (2, 4, 6):
    for i in idx:
        v, tail, p, _ = Ck(k, nus[i])
        res[(k, i)] = (v + tail, tail)
        pr(f"  k={k} node {i:2d} nu={nus[i]:8.4f}: raw={mp.nstr(v+tail, 15)}  tail/total={float(tail/(v+tail)):.2e} (decay p={p:.2f})  [{time.time()-t0:.0f} s]")
norm = res[(2, idx[0])][0]/ex[2][idx[0]]
pr("## comparison with exactskel cache after ONE global normalisation (k=2, first node)")
for k in (2, 4, 6):
    for i in idx:
        mine = res[(k, i)][0]/norm
        pr(f"  k={k} nu={nus[i]:8.4f} x={rec['x'][i]:.4f}: referee {float(mine):+.12e}  cache {ex[k][i]:+.12e}  rel.diff {float(mine/ex[k][i]-1):+.2e}"
           f"   (bulk-point E_k = {float((2-rec['x'][i]**2)/(1-rec['x'][i]**2)**(1+k/2)):.4e})")
open(os.path.join(HERE, 'cache', 'rf3_%d_%g.txt' % (Dg, mfl)), 'w').write("\n".join(out) + "\n")
