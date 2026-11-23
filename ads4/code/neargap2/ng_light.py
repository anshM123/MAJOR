"""
ng_light (NEARGAP2.md §3): images of LIGHT single traces (J = 0 and J = 2, either parity; ++ structure only, no +- below J = 4) under
Omega^EE_theta (Result D) and under A1PRIME's repaired Q functional Omega^Q_up + eps* Omega^EE_theta (light J = 0; the k = 4, 5 parts of
Omega^Q_up drop out, A1PRIME §3), at finite Delta_L < Delta_g, with the per-state machinery of the heavy states:
  hybrid image  = exact flat image (NF02 / Q3 weights, exact Wigner d^J_00 = P_J(1 - 2y), y = p^2 Lam^2/m^2)     [P: polynomial]
                + leading-Regge AdS curvature  2 psi2 (P_{-1/2+i nu}(cosh b) - J0(nu B)) + u psi3 (same),        [model of the AdS image]
                  cosh b = 1 + (2J+1)^2/(2 m^2),  B = 2 sinh(b/2) = (2J+1)/m,  nu = Lam p   (ee1000_common.hyb / q3, read-only)
  Lemma C (POS_ERRORS §2.1): |P - J0| <= (1 - sech(b/2)) + min(2, nu (B - b))  -> sign-free bound on the curvature.
Light-state mass: m^2 = (Delta_L - J - 2)(Delta_L + J - 1) (the heavy-state convention, alpha (alpha + 2J + 1)); light states with
Delta_L <= J + 2 have m^2 <= 0 and are reported through the m -> 0 limit only.
Normalisation (CHLPS/A1PRIME units): F = m^2 (B2 + B3) per |g0|^2 (J = 0) and per |g2|^2/m^4 (J = 2) [m in Lam units];
N = m^2(|B2| + |B3|).  Checks: F_D^0 = 11.317027 + 20.771893 m^2, F_D^2 = 3 - 6.766265 m^2 - 6.905463 m^4 + 20.771893 m^6,
F_Q^0 = -0.430430 + 2.032073 m^2 (A1PRIME §4.2-4.3).
Usage: python3 ng_light.py   -> cache/ng_light.txt, cache/ng_light.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import ng_guard  # noqa: F401
os.chdir(SF)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import mpmath as mp
import ee1000_common as E
Q = E.Q
NGC = os.path.join(HERE, 'cache')
QJ = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
PSI2Q = np.array([float(v) for v in QJ['psi2']]); PSI3Q = np.array([float(v) for v in QJ['psi3']])
EPS = 0.4304302/11.3170270          # eps* of A1PRIME


def Pleg(nus, eta):
    """P_{-1/2+i nu}(eta), mpmath hyp2f1 (exact; same formula as remk.Pders_vec n = 0)."""
    z = (1 - mp.mpf(eta))/2
    return np.array([float(mp.hyp2f1(mp.mpf(1)/2 + 1j*mp.mpf(v), mp.mpf(1)/2 - 1j*mp.mpf(v), 1, z).real) for v in nus])


def Pleg_md(nus, beta, nq=None):
    """P_{-1/2+i nu}(cosh beta) by the Mehler-Dirichlet integral (sqrt2/pi) int_0^beta cos(nu t) dt/sqrt(cosh b - cosh t), with
    t = beta (1 - v^2) (removes the endpoint singularity), Gauss-Legendre panels resolving the phase nu beta (vectorised)."""
    nus = np.asarray(nus, float)
    nq = nq or int(min(40000, 64*np.ceil((1.5*nus.max()*beta + 1000)/64)))
    xg, wg = np.polynomial.legendre.leggauss(64)
    k = nq//64; e = np.linspace(0, 1, k + 1)
    v = np.concatenate([(e[i] + e[i + 1])/2 + (e[i + 1] - e[i])/2*xg for i in range(k)])
    w = np.concatenate([(e[i + 1] - e[i])/2*wg for i in range(k)])
    t = beta*(1 - v*v)
    den = np.sqrt(np.maximum(np.cosh(beta) - np.cosh(t), 1e-300))
    g = w*2*beta*v/den
    out = np.empty_like(nus)
    for i0 in range(0, len(nus), 64):
        out[i0:i0 + 64] = np.cos(np.outer(nus[i0:i0 + 64], t)) @ g
    return np.sqrt(2)/np.pi*out


def pieces(Dg, J, m2abs, npn=None):
    """flat and LR-curvature pieces of a light (or heavy) ++ state of absolute mass^2 m2abs at Delta_g, for both weight sets."""
    Lm = E.lam_of(Dg); m = np.sqrt(m2abs); mfl = m/Lm; u = 1/mfl**2
    F = Q.NearGap(Dg).F
    # p-grid: Gauss-Legendre panels (the light integrand has P_J(1 - 2 u p^2), degree 2J in p, times smooth weights)
    if npn is None:                     # resolve the Bessel/Legendre oscillation: nu B <= Lam (2J+1)/m (about 0.8 nodes per radian + 256)
        npn = int(min(16384, 64*np.ceil((0.8*Lm*(2*J + 1)/m + 256)/64)))
    xg, wg = np.polynomial.legendre.leggauss(64)
    edges = np.linspace(0, 1, npn//64 + 1)
    pn = np.concatenate([(edges[i] + edges[i + 1])/2 + (edges[i + 1] - edges[i])/2*xg for i in range(len(edges) - 1)])
    pw = np.concatenate([(edges[i + 1] - edges[i])/2*wg for i in range(len(edges) - 1)])
    y = u*pn*pn; x = np.sqrt(y); nu = Lm*pn
    D = np.polynomial.legendre.legval(1 - 2*y, [0]*J + [1])
    beta = 2*np.arcsinh((J + 0.5)/m); B = 2*np.sinh(beta/2)
    eta = np.cosh(beta)
    P = Pleg_md(nu, beta); from scipy.special import j0 as _j0; jz = _j0(nu*B)
    dP = P - jz
    lemC = (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, nu*(B - beta))
    out = dict(J=J, m2=m2abs, mfl=mfl, beta=beta, Lam=Lm)
    for tag, p2, p3 in (('D', F.psi2(pn), F.psi3(pn)),
                        ('Q', np.polynomial.polynomial.polyval(pn, PSI2Q), np.polynomial.polynomial.polyval(pn, PSI3Q))):
        B2 = (pw*p2*(2 - y)*D).sum(); B3 = u*(pw*p3*D).sum()
        cv2 = (pw*p2*2*dP).sum(); cv3 = u*(pw*p3*dP).sum()
        cb = (pw*np.abs(p2)*2*lemC).sum() + u*(pw*np.abs(p3)*lemC).sum()
        ca = (pw*np.abs(p2)*2*np.abs(dP)).sum() + u*(pw*np.abs(p3)*np.abs(dP)).sum()
        f = mfl**2*(mfl**4 if J == 2 else 1.0)        # CHLPS units: per |g0|^2, per |g2|^2/m^4 (m in Lam units)
        out[tag] = dict(B2=B2*f, B3=B3*f, flat=(B2 + B3)*f, N=(abs(B2) + abs(B3))*f, cv=(cv2 + cv3)*f, cvabs=ca*f, cvbound=cb*f)
    return out


def m2_of(Dl, J):
    return (Dl - J - 2.0)*(Dl + J - 1.0)


if __name__ == '__main__':
    t0 = time.time()
    lines = []
    # check of the Mehler-Dirichlet evaluator against mpmath hyp2f1 where the latter converges
    chk = []
    for beta in (0.01, 0.3, 1.0):
        nus = np.array([0.5, 7.3, 55.0, 300.0, 998.0])
        try:
            ref = Pleg(nus, np.cosh(beta)); chk.append(float(np.max(np.abs(Pleg_md(nus, beta) - ref))))
        except Exception:
            pass

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    res = []
    # flat checks against A1PRIME's exact polynomials (Dg only sets Lam; flat pieces are Lam-free in these units)
    pr("# ng_light: light J = 0, 2 (++) images. Units: per |g0|^2 (J=0), per |g2|^2/m^4 (J=2); m in units of Lam.")
    pr(f"# Mehler-Dirichlet P vs mpmath hyp2f1 (beta = 0.01, 0.3, 1; nu = 0.5 ... 998): max abs diff {max(chk):.1e}")
    pr("# flat check vs A1PRIME polynomials (m/Lam = 0.1, 0.4602, 0.675, 0.99):")
    for mf in (0.1, 0.460237, 0.675, 0.99):
        Lm = E.lam_of(1000.0); r0 = pieces(1000.0, 0, (mf*Lm)**2); r2 = pieces(1000.0, 2, (mf*Lm)**2)
        q = mf*mf
        pr(f"   m={mf:.4f}: F_D^0 {r0['D']['flat']:+.6f} (poly {11.317027 + 20.771893*q:+.6f}); F_Q^0 {r0['Q']['flat']:+.6f} "
           f"(poly {-0.430430 + 2.032073*q:+.6f}); F_D^2 {r2['D']['flat']:+.6f} (poly {3 - 6.766265*q - 6.905463*q*q + 20.771893*q**3:+.6f})")
    for Dg in (1000.0, 300.0):
        Lm = E.lam_of(Dg)
        pr(f"\n## Delta_g = {Dg:g} (Lam = {Lm:.3f}); columns: Delta_L, m/Lam, beta; per functional: flat F, F/N flat, AdS(LR) F, "
           f"(LR curvature)/N, Lemma-C bound/N, worst-case (flat - bound)/N")
        DLs = sorted(set([2.6, 3.0, 3.5, 4.5, 5, 6, 8, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300] +
                         [round(Dg*f, 3) for f in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999)]))
        for J in (0, 2):
            pr(f"  J = {J}:")
            for Dl in DLs:
                if Dl >= Dg or m2_of(Dl, J) <= 0:
                    continue
                r = pieces(Dg, J, m2_of(Dl, J))
                row = dict(Dg=Dg, J=J, Dl=Dl, mfl=r['mfl'], beta=r['beta'])
                s = f"   D_L={Dl:8.3f} m/L={r['mfl']:.5f} b={r['beta']:.4f} |"
                tags = ('D', 'Q', 'R') if J == 0 else ('D',)
                for tag in tags:
                    if tag == 'R':                     # repaired Q: Q + eps* D
                        a = {k: r['Q'][k] + EPS*r['D'][k] for k in ('flat', 'cv', 'cvbound')}
                        a['N'] = r['Q']['N'] + EPS*r['D']['N']
                    else:
                        a = r[tag]
                    ads = a['flat'] + a['cv']
                    row[tag] = dict(flat=a['flat'], N=a['N'], ads=ads, cvN=a['cv']/a['N'], cbN=a['cvbound']/a['N'],
                                    worst=(a['flat'] - a['cvbound'])/a['N'])
                    s += f" {tag}: {a['flat']:+.5f} {a['flat']/a['N']:+.4f} AdS {ads:+.5f} cv/N {a['cv']/a['N']:+.2e} C/N {a['cvbound']/a['N']:.2e} w {(a['flat'] - a['cvbound'])/a['N']:+.4f} |"
                pr(s)
                # quadrature check: same state with twice the nodes
                if Dl in (2.6, 3.0, 5, 20, 100) or Dl > 0.98*Dg:
                    r2x = pieces(Dg, J, m2_of(Dl, J), npn=2*int(min(16384, 64*np.ceil((0.8*Lm*(2*J + 1)/np.sqrt(m2_of(Dl, J)) + 256)/64))))
                    row['quadcheck'] = {t: (r2x[t]['cv'] - r[t]['cv'])/r[t]['N'] for t in ('D', 'Q')}
                    pr(f"      quadrature check (2x nodes): d(cv)/N  D {row['quadcheck']['D']:+.1e}  Q {row['quadcheck']['Q']:+.1e}")
                res.append(row)
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(NGC, 'ng_light.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(NGC, 'ng_light.json'), 'w'))
