"""
fullpos4/fp4_soft.py (FULLPOS4.md §3.3) [DV]: soft-signed / hard-absolute norm of the frequency-integrated outgoing §55 contour.
For each alpha = A - i eps the soft integral is done SIGNED along the whole gamma path:
   Psi(al) = int_path e^{-i a ga} xi^{s-1}(1-xi^2)^{-p} P+_h(beta) d ga,
   A >= 0: ga = C - i eps, C: A -> pi/2, then ga = pi/2 - i t, t >= eps   (= I_U - i I_T, the corner is paired automatically)
   A <  0: ga = al - i v, v >= 0                                          (= -i I_L)
and  L(G) = (C/pi) Re int dA e^{i S al} Psi(al)/sin(al).
Reported: N_soft = (C/pi) int |e^{i S al} Psi(al)/sin al| dA (this is what an alpha-cell certificate pays if the soft integral is
evaluated exactly per alpha), split into |A| bands, against |G(c)| and N_s.  Weight psi2 (main part), p = 3.
Usage: fp4_soft.py a J b [nA_per_panel]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import fp4_int as I
from fp4_pf import Cconst
from fp4_abs import gl, graded, geom


def soft_path(al, A, a, eps, lam):
    """nodes, d gamma weights for the soft path of al (A >= 0: arc + tail; A < 0: lower ray)."""
    if A >= 0:
        L = math.pi/2 - A
        wmax = math.pi/a                           # half-wavelength panels, 8 GL nodes
        e = [0.0]; h = max(abs(al)*1e-6, 1e-13)
        while e[-1] + h < L:
            e.append(e[-1] + h); h = min(h*1.3, wmax)
        e.append(L)
        X, Wt = zip(*[gl(lo, hi, 8) for lo, hi in zip(e[:-1], e[1:])])
        d = np.concatenate(X); wd = np.concatenate(Wt)
        g1 = A + d - 1j*eps; w1 = wd + 0j
        tmax = min(3.0, 60.0/max(a - lam, 20.0))
        t, wt = graded(0.0, tmax, 1e-5, 1.3)
        g2 = math.pi/2 - 1j*(eps + t); w2 = -1j*wt
        return np.concatenate([g1, g2]), np.concatenate([w1, w2])
    v, wv = graded(0.0, 60.0/a, max(abs(al)*1e-6, 1e-13), 1.3)
    return al - 1j*v, -1j*wv


def run(a, J, b, nA=8, p=3):
    t0 = time.time()
    amp = I.Amp(2, nq=12)
    n = J - b; N = n + a; s = 2*b + 1; S = 2*J + 1 + a; eps = 1.0/S
    Delta = N + b + 2; T = (Delta - 1)*(Delta - 2); C = J*(J + 1); c = (C - b*(b + 1))/(T - b*(b + 1)); w = 1 - C/T
    Cc = Cconst(n, N, s)
    Ap, wAp = graded(0.0, math.pi/2, eps/20, 1.6, nA)
    As = np.concatenate([-Ap[::-1], Ap]); wA = np.concatenate([wAp[::-1], wAp])
    rows = []
    for A, wa in zip(As, wA):
        al = A - 1j*eps
        g, wg = soft_path(al, A, a, eps, amp.lam)
        xi, beta = geom(np.full(len(g), al), g)
        ok = np.abs(xi) > 1e-13
        Pr = np.zeros(len(g), complex); Et = np.zeros(len(g), complex)
        Pr[ok], Et[ok] = amp(beta[ok], factored=True)
        pre = xi**(s - 1)*(1 - xi*xi)**(-p)
        car = -1j*a*g
        f = pre*(Pr*np.exp(car) + Et*np.exp(car + 1j*amp.lam*beta))
        Psi = complex((f*wg).sum()); Pabs = float((np.abs(f)*np.abs(wg)).sum())
        hard = np.exp(1j*S*al)/np.sin(al)
        rows.append((A, wa, Psi*hard, Pabs*abs(hard)))
    rows = np.array(rows, dtype=complex)
    Aa = rows[:, 0].real; wa = rows[:, 1].real; v = rows[:, 2]; va = rows[:, 3].real
    Nsoft = Cc/math.pi*(np.abs(v)*wa).sum(); Nraw = Cc/math.pi*(va*wa).sum()
    bands = {}
    for lo, hi in [(0, 1), (1, 10), (10, 100), (100, 1000), (1000, 1e9)]:
        m = (np.abs(Aa)*S >= lo) & (np.abs(Aa)*S < hi)
        bands[f'|A|S in [{lo},{hi})'] = dict(soft=Cc/math.pi*(np.abs(v[m])*wa[m]).sum(), raw=Cc/math.pi*(va[m]*wa[m]).sum())
    import ls_core as LSC
    beta_c = math.acosh((1 + c)/(1 - c))
    Gc = (1 - c)**(-3)*2*LSC.OmLR(2, [beta_c])[0]
    m_ = math.sqrt(T - C); Ns = LSC.hyb_fast(1000, J, m_/I.LAM, '++' if b == 0 else '+-')['Ns']
    return dict(a=a, J=J, b=b, S=S, c=c, L=-math.log(1 - c), C=Cc, Gc=Gc, Ns=Ns, w3=w**3, N_soft=Nsoft, N_raw=Nraw,
                soft_over_Gc=Nsoft/abs(Gc), raw_over_Gc=Nraw/abs(Gc), w3soft_over_Ns=w**3*Nsoft/Ns, bands=bands,
                nalpha=len(As), secs=time.time() - t0)


if __name__ == '__main__':
    a, J, b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    nA = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    r = run(a, J, b, nA)
    print(json.dumps(r), flush=True)
    open(os.path.join(HERE, 'cache', 'fp4_soft.jsonl'), 'a').write(json.dumps(r) + '\n')
