"""agproof/ag_lock.py (AG_PROOF.md §3): is the exact scalar's beyond-hybrid correction G = exact - hybrid, near the edge nu = Lam,
a pair of travelling waves PHASE-LOCKED to the flat waves of Lemma F?  If so, its + wave on the right side p = 1 of Lemma P4W's rectangle
can be extracted and its decay compared with the AG' majorant.

Model (tested, not assumed): on a disc D = {|nu - nu_c| <= R} around the edge,
    G(nu) = cbar(nu) flat_-(nu) + c(nu) flat_+(nu),   cbar(nu) = sum_j b_j s^j,  c(nu) = sum_j d_j s^j,  s = (nu - nu_c)/R,
(flat_pm = E_2 W_pm(theta(nu))/2, Lemma F).  Real-analyticity of G predicts d_j = conj(b_j); this is NOT imposed in the fit and is
reported as a consistency check.  (b, d) by linear least squares on exact values of G at complex nu in D (real axis, upper and lower
half; the exact scalar k = 2 moment at complex nu is the entire Laplace/Barnes series of ag_exact.py).  The fit residual measures the
model error.  Then
    G_+ := c flat_+   on the right side 1 + i y,  0 <= y <= y_D (the part of the side inside D),
and its decay -log|G_+(1+iy)/G_+(1)| is compared with Lemma F's flat profile and P4_WINDOW's AG' majorant (local rate
0.97 Lam beta/sqrt(1 - x1^2) up to y = 15/kappa).  The edge integral V is also formed (inside D; beyond y_D the flat-profile is used,
scaled by |c(1+iy_D)|).
-> cache/lock.txt, cache/lock.json      Usage: python3 ag_lock.py [J ...]   (default 40 48 64 100)
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
from scipy.special import jv
import ag_exact as AE
import ag_flat as AF
import p4w_common as C            # read-only (imported through ag_flat's path)
import p4w_contour as K           # read-only

mp.mp.dps = 30
L = AE.LAM
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def flat_pm(J, m, nu):
    """(flat_+, flat_-) at complex nu: the convergent wave by its series, the other as P_J - (it)."""
    x = mp.mpc(nu)/m
    th = 2*mp.asin(x)
    E2 = (2 - x*x)/(1 - x*x)**2
    P = mp.legendre(J, 1 - 2*x*x)
    if mp.im(th) >= 0:
        Wp = AF.Wplus(J, th); Wm = P - Wp
    else:
        Wm = AF.Wminus(J, th); Wp = P - Wm
    return complex(E2*Wp/2), complex(E2*Wm/2)


def curv(st, nu):
    """P_{-1/2+i nu}(cosh beta) - J0(bhat p) at complex nu (the hybrid's leading-Regge curvature, per unit leading image)."""
    nu = np.atleast_1d(np.asarray(nu, complex))
    P = C.Pd(st['beta'], nu)[0]
    bh = (2*st['J'] + 1)/st['mfl']
    return P - jv(0, bh*nu/L)


def analyse(J, deg=10, Rfac=None):
    st = AE.state(J); m = st['m']; beta = st['beta']; x1 = L/m
    k0 = 0.97*beta; kl = k0/np.sqrt(1 - x1*x1)             # nu-units
    T = 25.0/k0
    dist = m - L                                           # distance from p = 1 to the flat branch point nu = m
    off = float(os.environ.get('AG_OFF', '1.5'))           # disc centre at Lam - off*dist; R = 0.8 (m - centre): D avoids nu = m
    nuc = L - off*dist
    R = 0.8*(m - nuc)
    t0 = time.time()
    L3 = AE.build(st, numax=float(np.hypot(L, T))*1.02, prec=400, conf=(320, 960, 0.5))
    # fit points: polar grid in D (both half-planes and the real axis)
    pts = []
    for rr in np.linspace(0.05, 1.0, 14):
        nang = max(8, int(40*rr))
        for a in np.linspace(-np.pi, np.pi, nang, endpoint=False):
            pts.append(nuc + rr*R*np.exp(1j*a))
    pts += list(nuc + np.linspace(-R, R, 41))              # real axis
    pts = np.array(pts)
    pts = pts[np.abs(pts - m) > 0.15*dist]                  # stay off the flat pole at nu = m
    ex = np.array([complex(float(v.real.mid()), float(v.imag.mid())) for v in (AE.ex2_complex(L3, z) for z in pts)])/2
    fp = np.array([flat_pm(J, m, z) for z in pts])
    fl = fp[:, 0] + fp[:, 1]
    G = ex - fl                                            # g = exact - flat (AG' is a statement about g)
    if os.environ.get('AG_SUBCURV'):
        G = G - curv(st, pts)                               # G = exact - hybrid
    s = (pts - nuc)/R
    out = dict(J=J, mfl=st['mfl'], x1=x1, beta=beta, R=R, nuc=nuc, dist=dist, kl=kl, fits={})
    for dg in sorted(set([deg - 2, deg, deg + 2])):
        A = np.hstack([fp[:, 1:2]*s[:, None]**np.arange(dg + 1), fp[:, 0:1]*s[:, None]**np.arange(dg + 1)])
        env = np.abs(fp[:, 0]) + np.abs(fp[:, 1])            # residual relative to the local flat-wave envelope
        w = 1/env
        co, *_ = np.linalg.lstsq(A*w[:, None], G*w, rcond=None)
        b = co[:dg + 1]; d = co[dg + 1:]
        cb_ = A[:, :dg + 1] @ b; c_ = A[:, dg + 1:] @ d
        res = np.abs(A @ co - G)/(np.abs(cb_) + np.abs(c_))
        onax = np.abs(pts.imag) < 1e-9
        cons = float(np.max(np.abs(d - np.conj(b)))/np.max(np.abs(b)))
        # right side p = 1 inside D
        yD = float(np.sqrt(max(R*R - (L - nuc)**2, 0.0)))
        ys = np.linspace(0.0, yD, 121)
        zz = L + 1j*ys
        sz = (zz - nuc)/R
        cz = np.polyval(d[::-1], sz); cbz = np.polyval(b[::-1], sz)
        fpz = np.array([flat_pm(J, m, z) for z in zz])
        Gp = cz*fpz[:, 0]
        decay = -np.log(np.abs(Gp/Gp[0]))
        flatd = -np.log(np.abs(fpz[:, 0]/fpz[0, 0]))
        maj = np.where(ys <= 15/kl, kl*ys, 15 + k0*(ys - 15/kl))
        # real-axis amplitude of G at the edge: |G_+(Lam)| vs the full |G| envelope near p = 1
        out['fits'][dg] = dict(res_max=float(res.max()), res_axis=float(res[onax].max()), cons=cons, yD=yD, yD_kl=yD*kl,
                               cabs=[float(v) for v in np.abs(cz)], ys=[float(v) for v in ys], decay=[float(v) for v in decay],
                               flatd=[float(v) for v in flatd], maj=[float(v) for v in maj],
                               gamma=float(np.max(np.abs(cz))/np.abs(cz[0])), Gp0=float(abs(Gp[0])),
                               min_excess=float(np.min(decay - maj)), min_excess_flat=float(np.min(decay - flatd)))
        # edge integral (k = 2, p-units): inside D with G_+, beyond y_D the flat profile scaled at y_D (|c| frozen)
        yp = ys/L
        ps = np.abs(K.psi_c(2, 1 + 1j*yp))
        VG = np.trapezoid(ps*np.abs(Gp/Gp[0]), yp)
        VAGin = np.trapezoid(ps*np.exp(-maj), yp)
        out['fits'][dg].update(VG_in=float(VG), VAG_in=float(VAGin))
    out['t'] = time.time() - t0
    return out


def main():
    Js = [int(a) for a in sys.argv[1:]] or [40, 48, 64, 100]
    pr(__doc__.split('->')[0])
    res = {}
    for J in Js:
        o = analyse(J)
        res[J] = o
        pr(f"\n== J = {J} (gap), m/L = {o['mfl']:.4f}, x1 = {o['x1']:.4f}, beta = {o['beta']:.4f}; distance p=1 -> nu=m: {o['dist']:.1f}; "
           f"disc centre {o['nuc']:.1f}, R = {o['R']:.1f}; kappa_loc = {o['kl']:.4f}  [{o['t']:.0f}s]")
        for dg, f in o['fits'].items():
            pr(f"   deg {dg:2d}: rel. residual max {f['res_max']:.1e} (real axis {f['res_axis']:.1e}); |d - conj b|/|b| = {f['cons']:.1e}; "
               f"side p=1 inside D up to y = {f['yD']:.1f} = {f['yD_kl']:.1f}/kappa_loc; gamma = max|c|/|c(1)| = {f['gamma']:.3f}; "
               f"min_y [decay(G_+) - AG' majorant] = {f['min_excess']:+.2f}, - flat profile: {f['min_excess_flat']:+.2f}; "
               f"V(G_+)/V(AG') inside D = {f['VG_in']/f['VAG_in']:.3f}")
        f = o['fits'][max(o['fits'])]
        ys = np.array(f['ys']); kl = o['kl']
        sel = [int(np.argmin(np.abs(ys - c/kl))) for c in (1, 2, 5, 8, 10, 12, 15) if c/kl <= ys[-1]]
        pr('   profile at y*kappa_loc = ' + ', '.join(f"{ys[i]*kl:.1f}: G_+ {f['decay'][i]:.2f} | flat {f['flatd'][i]:.2f} | AG' {f['maj'][i]:.2f} | |c|/|c(1)| {f['cabs'][i]/f['cabs'][0]:.2f}" for i in sel))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(AE.CACHE, 'lock.json'), 'w'), default=float)
    open(os.path.join(AE.CACHE, 'lock.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
