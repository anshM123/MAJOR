"""
fullpos4/fp4_diag.py (FULLPOS4.md §4) [DV]: the author's §74 diagonal extraction + the remaining outgoing contour.

Soft variable xi: Y = e^{i ga} = (X - xi)/(1 - xi X), X = e^{i al} = e^{eps + iA};  d ga/sin al = 2 d xi/((1 - xi X)(1 - xi/X)),
carrier e^{i S al - i a ga} = X^B ((1 - xi X)/(1 - xi/X))^a.   With K(A, xi) = X^B (1 - xi X)^{a-1} (1 - xi/X)^{-a-1}:
   L(G) = (C/pi) int dA int_{xi-path(A)} K xi^{s-1} F(xi^2) d xi        (full F)
        = 2 Re (C/pi) int dA int K xi^{s-1} F+(xi^2) d xi               (outgoing; xi-path = image of the §55 gamma path)
Split each soft path at its LAST crossing xi_c(A) of |xi| = xi0:  [0, xi0] real  +  arc |xi| = xi0 from xi0 to xi_c  +  rest of path.
   D = (C/pi) int_0^{xi0} xi^{s-1} G(xi^2) h_B(xi) d xi,   h_B(xi) = int_{-pi/2}^{pi/2} K(A, xi) dA   (real G: outgoing logs cancel)
   R = 2 Re (C/pi) int dA [ int_arc K xi^{s-1} F+ d xi + int_{rest of gamma path} e^{i S al - i a ga} xi^{s-1} F+ / sin al d ga ]
   L(G) = D + R.
Weight psi2 (p = 3) via fp4_int.Amp (lam = Lam or a toy bandwidth).  Reported: D, its norm, N_rem (soft-signed per A), N_rem_raw.
Usage: fp4_diag.py a J b [nA] [lam] [xi0] [signed_check]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import fp4_int as I
from fp4_pf import Cconst
from fp4_abs import gl, graded, geom


def logK(A, xi, B, a, eps):
    X = np.exp(eps + 1j*A)
    return B*(eps + 1j*A) + (a - 1)*np.log(1 - xi*X) - (a + 1)*np.log(1 - xi/X)


def hB(xi, B, a, eps, nper=10):
    """h_B(xi) = int_{-pi/2}^{pi/2} K dA by Gauss-Legendre panels of a quarter wavelength of e^{iBA} (direct)."""
    npan = int(4*B*math.pi/(2*math.pi)*1.0) + 8
    e = np.linspace(-math.pi/2, math.pi/2, npan + 1)
    x, w = np.polynomial.legendre.leggauss(nper)
    A = (e[:-1, None] + (e[1:] - e[:-1])[:, None]*(x[None, :] + 1)/2).ravel()
    wA = ((e[1:] - e[:-1])[:, None]/2*w[None, :]).ravel()
    out = np.zeros(len(xi), complex)
    for i, z in enumerate(xi):
        out[i] = (np.exp(logK(A, z, B, a, eps))*wA).sum()
    return out


def hB_closed(xi, B, a, eps, J, jmax=40):
    """author's closed form (§74): (2(-1)^J e^{B eps}/B) Re sum_j (-1)^j Z^{(j)}(eps)/B^j, Z(v) = exp[(a-1)log(1 - i xi e^v) - (a+1)log(1 + i xi e^{-v})]
    (Taylor coefficients by mpmath)."""
    import mpmath as mp
    mp.mp.dps = 40
    out = []
    for z in xi:
        Zf = lambda v: mp.exp((a - 1)*mp.log(1 - 1j*z*mp.exp(v)) - (a + 1)*mp.log(1 + 1j*z*mp.exp(-v)))
        tc = mp.taylor(Zf, eps, jmax)
        s = sum((-1)**j*mp.factorial(j)*tc[j]/mp.mpf(B)**j for j in range(jmax + 1))
        out.append(float(2*(-1)**J*mp.exp(B*eps)/B*mp.re(s)))
    return np.array(out)


def path_param(A, a, eps, lam):
    """soft gamma path of al = A - i eps as a function of a real parameter tau in [0, tau_end]; returns (gamma(tau), dgamma/dtau, tau_end, L)."""
    al = A - 1j*eps
    if A >= 0:
        L = math.pi/2 - A
        tmax = min(3.0, 60.0/max(a - lam, 20.0))
        def g(t):
            t = np.asarray(t, float)
            return np.where(t <= L, A + t - 1j*eps, math.pi/2 - 1j*(eps + (t - L)))
        def dg(t):
            t = np.asarray(t, float)
            return np.where(t <= L, 1.0 + 0j, -1j + 0*t)
        return g, dg, L + tmax, L
    vmax = max(60.0/a, 0.2)                   # far enough that the path leaves |xi| < xi0 (crossing must exist)
    return (lambda t: al - 1j*np.asarray(t, float)), (lambda t: -1j + 0*np.asarray(t, float)), vmax, None


def nodes_from(t0, t1, a, Lsplit, h0):
    """GL nodes on [t0, t1]: graded from t0 (h0, ratio 1.3) up to half-wavelength panels pi/a; split at Lsplit (arc/tail corner)."""
    wmax = math.pi/a
    e = [t0]; h = h0
    while e[-1] + h < t1:
        nxt = e[-1] + h
        if Lsplit is not None and e[-1] < Lsplit < nxt:
            nxt = Lsplit
        e.append(nxt); h = min(h*1.3, wmax)
    e.append(t1)
    pan = [gl(lo, hi, 8) for lo, hi in zip(e[:-1], e[1:]) if hi > lo]
    if not pan:
        return np.array([]), np.array([])
    X, Wt = zip(*pan)
    return np.concatenate(X), np.concatenate(Wt)


def crossing(A, a, eps, lam, xi0):
    """last crossing tau_c of |xi(gamma(tau))| = xi0 along the soft path."""
    g, dg, tend, L = path_param(A, a, eps, lam)
    al = A - 1j*eps
    ts = np.concatenate([np.geomspace(1e-14, tend, 4000)])
    xi, _ = geom(np.full(len(ts), al), g(ts))
    inside = np.abs(xi) < xi0
    k = np.nonzero(inside)[0]
    k = k[-1] if len(k) else 0
    lo, hi = (ts[k], ts[k + 1]) if k + 1 < len(ts) else (ts[k], tend)
    for _ in range(80):
        mid = 0.5*(lo + hi)
        x, _ = geom(np.array([al]), g(np.array([mid])))
        if abs(x[0]) < xi0:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def run(a, J, b, nA=4, lam=None, xi0=0.02, p=3, signedA=None):
    t0 = time.time()
    amp = I.Amp(2, lam=lam, nq=12)
    LAMw = amp.lam
    n = J - b; N = n + a; s = 2*b + 1; S = 2*J + 1 + a; B = 2*J + 1; eps = 1.0/S
    Delta = N + b + 2; T = (Delta - 1)*(Delta - 2); C = J*(J + 1); c = (C - b*(b + 1))/(T - b*(b + 1)); w = 1 - C/T
    Cc = Cconst(n, N, s)
    # ---- diagonal block D
    xd, wd = graded(0.0, xi0, xi0*1e-4, 1.3, 16)
    hb = hB(xd, B, a, eps)
    betad = np.log((1 + xd)/(1 - xd)) + 0j
    Gd = (1 - xd**2)**(-p)*2*np.real(amp(betad))
    integ = xd**(s - 1)*Gd*hb.real
    D = Cc/math.pi*(integ*wd).sum(); Dnorm = Cc/math.pi*(np.abs(integ)*wd).sum()
    him = float(np.max(np.abs(hb.imag))/max(np.max(np.abs(hb.real)), 1e-300))
    # closed-form check of h_B at 3 points
    try:
        hc = hB_closed(np.array([xi0*0.1, xi0*0.5, xi0]), B, a, eps, J)
        hd = hB(np.array([xi0*0.1, xi0*0.5, xi0]), B, a, eps).real
        hcheck = float(np.max(np.abs(hc - hd))/np.max(np.abs(hd)))
    except Exception as ex:
        hcheck = str(ex)
    # ---- remaining contour per A
    if signedA is None:
        Ap, wAp = graded(0.0, math.pi/2, eps/20, 1.6, nA)
    else:
        Ap, wAp = signedA(S, eps)
    As = np.concatenate([-Ap[::-1], Ap]); wA = np.concatenate([wAp[::-1], wAp])
    vals = []; raws = []
    for A in As:
        al = A - 1j*eps; X = np.exp(eps + 1j*A)
        tc = crossing(A, a, eps, LAMw, xi0)
        g, dg, tend, L = path_param(A, a, eps, LAMw)
        xc, _ = geom(np.array([al]), g(np.array([tc])))
        th_c = float(np.angle(xc[0]))
        # arc xi = xi0 e^{i th}, th: 0 -> th_c
        th, wth = gl(0.0, th_c, 24) if abs(th_c) > 0 else (np.array([]), np.array([]))
        xa = xi0*np.exp(1j*th); dxa = 1j*xa*wth
        Pr, Et = amp(np.log((1 + xa)/(1 - xa)), factored=True)
        ba = np.log((1 + xa)/(1 - xa))
        lk = logK(A, xa, B, a, eps)
        fa = 2*xa**(s - 1)*(1 - xa*xa)**(-p)*(Pr*np.exp(lk) + Et*np.exp(lk + 1j*LAMw*ba))/((1 - xa*X)*(1 - xa/X))*0
        # (K already contains (1-xi X)^{a-1}(1-xi/X)^{-a-1} = carrier/((1-xi X)(1-xi/X)); the 2 is the d ga/sin al factor)
        fa = 2*xa**(s - 1)*(1 - xa*xa)**(-p)*(Pr*np.exp(lk) + Et*np.exp(lk + 1j*LAMw*ba))
        Iarc = (fa*dxa).sum(); Rarc = (np.abs(fa)*np.abs(dxa)).sum()
        # rest of the gamma path
        tt, wt = nodes_from(tc, tend, a, L, max(tc*1e-3, 1e-12))
        if len(tt) == 0:
            vals.append(Iarc); raws.append(Rarc); continue
        gg = g(tt); dgg = dg(tt)*wt
        xi, beta = geom(np.full(len(gg), al), gg)
        Pr, Et = amp(beta, factored=True)
        car = 1j*S*al - 1j*a*gg
        f = xi**(s - 1)*(1 - xi*xi)**(-p)*(Pr*np.exp(car) + Et*np.exp(car + 1j*LAMw*beta))/np.sin(al)
        Irest = (f*dgg).sum(); Rrest = (np.abs(f)*np.abs(dgg)).sum()
        vals.append(Iarc + Irest); raws.append(Rarc + Rrest)
    vals = np.array(vals); raws = np.array(raws)
    Rsigned = 2*Cc/math.pi*(vals*wA).sum().real/2         # (C/2pi)*2Re: see header (int d ga/sin al -> 2 int d xi ...)
    Nrem = Cc/math.pi*(np.abs(vals)*wA).sum(); Nraw = Cc/math.pi*(raws*wA).sum()
    bands = {}
    for lo, hi in [(0, 1), (1, 10), (10, 100), (100, 1000), (1000, 1e9)]:
        m = (np.abs(As)*S >= lo) & (np.abs(As)*S < hi)
        bands[f'[{lo},{hi})'] = Cc/math.pi*(np.abs(vals[m])*wA[m]).sum()
    out = dict(a=a, J=J, b=b, S=S, c=c, L=-math.log(1 - c), lam=LAMw, xi0=xi0, C=Cc, D=D, Dnorm=Dnorm, hB_imag_rel=him, hB_closed_check=hcheck,
               R_signed=Rsigned, total=D + Rsigned, N_rem_soft=Nrem, N_rem_raw=Nraw, bands_rem_soft=bands, nalpha=len(As), w3=w**3)
    if LAMw > 100:
        import ls_core as LSC
        beta_c = math.acosh((1 + c)/(1 - c))
        Gc = (1 - c)**(-3)*2*LSC.OmLR(2, [beta_c])[0]
        out.update(Gc=Gc, D_over_Gc=D/Gc, Dnorm_over_Gc=Dnorm/abs(Gc), Nrem_over_Gc=Nrem/abs(Gc), Nraw_over_Gc=Nraw/abs(Gc))
        if J <= 4000:
            m_ = math.sqrt(T - C); Ns = LSC.hyb_fast(1000, J, m_/I.LAM, '++' if b == 0 else '+-')['Ns']
            out.update(Ns=Ns, w3Nrem_over_Ns=w**3*Nrem/Ns, w3D_over_Ns=w**3*D/Ns)
    out['secs'] = time.time() - t0
    return out


if __name__ == '__main__':
    a, J, b = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    a = int(a) if a == int(a) else a
    nA = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    lam = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != '-' else None
    xi0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.02
    r = run(a, J, b, nA, lam, xi0)
    print(json.dumps(r), flush=True)
    open(os.path.join(HERE, 'cache', 'fp4_diag.jsonl'), 'a').write(json.dumps(r) + '\n')
