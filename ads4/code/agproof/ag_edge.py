"""agproof/ag_edge.py (AG_PROOF.md §3.2): the EDGE CONTRIBUTION of the exact scalar's correction, measured on the real axis, against
the edge term of Lemma P4W.

For a function f on [0, 1] (p = nu/Lam) the edge (p = 1) contribution to int psi_2 f is isolated with an entire cutoff
    chi(p) = erfc((p_e - p)/sigma)/2   (0 below p_e, 1 above),     E[f] := int_0^1 psi_2(p) chi(p) f(p) dp,
on Gauss-Legendre nodes (exact values of f).  For f = f_+ + f_- with f_pm analytic travelling waves (decaying into the upper/lower
half-plane), E[f] = 2 Re[-i int_0^inf psi_2(1+iy) f_+(1+iy) dy] + O(e^{-(kappa sigma)^2/4}) + (cutoff region), i.e. exactly what the
right side of Lemma P4W's rectangle bounds.  Calibration: f = flat (Lemma F), where the right-side integral is computable exactly.
Objects (per unit leading image, Delta_phi = 3 exact scalar, gap states J = 40, 48, 64, 100):
   g = exact - flat,  G = exact - hybrid = g - (P - J0).
Compared with: P4W's edge term for G at Lam_AG = 1 (2 [V(1; kappa_1, T1 = 15) (e_sc(1) + A) + V(1; 0.97 Lam beta) A], k = 2, without the
ring weight 2), and with |g(p)| near p = 1.  Stability: several (p_e, sigma).
-> cache/edge.txt, cache/edge.json        Usage: python3 ag_edge.py [J ...]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
from scipy.special import erfc, j0, eval_legendre
import ag_exact as AE
import ag_flat as AF
import p4w_common as C          # read-only
import p4w_contour as K         # read-only

L = AE.LAM
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def psi2(p):
    return K.psi_c(2, np.asarray(p, float)).real


def main():
    Js = [int(a) for a in sys.argv[1:]] or [40, 48, 64, 100]
    pr(__doc__.split('->')[0])
    res = {}
    for J in Js:
        t0 = time.time()
        st = AE.state(J); m = st['m']; beta = st['beta']; mfl = st['mfl']; x1 = L/m
        L3 = AE.build(st, numax=L*1.02, prec=400, conf=(320, 960, 0.5))
        k0 = 0.97*beta*L; kl = k0/np.sqrt(1 - x1*x1)            # p-units
        T = 25.0/k0
        # nodes on [0.55, 1]: composite GL, dense enough for local frequency <= 1.2 kl
        nseg = int(np.ceil(0.45*1.2*kl/np.pi*2)) + 20
        xs, ws = np.polynomial.legendre.leggauss(16)
        br = np.linspace(0.55, 1.0, nseg + 1)
        p = np.concatenate([a + (b - a)*(xs + 1)/2 for a, b in zip(br[:-1], br[1:])])
        wp = np.concatenate([ws*(b - a)/2 for a, b in zip(br[:-1], br[1:])])
        nu = L*p; x = nu/m
        ex = np.array([float(L3.value_fast(2, float(v)).mid()) for v in nu])/2
        flat = (2 - x*x)/(1 - x*x)**2*eval_legendre(J, 1 - 2*x*x)/2
        Pn = C.Pd(beta, nu)[0].real
        cur = Pn - j0((2*J + 1)/mfl*p)
        g = ex - flat; G = g - cur
        ps = psi2(p)
        # exact right-side integral for the flat wave (Lemma F): 2 Re[-i int_0^Y psi2(1+iy) flat_+(1+iy) dy], Y = 60/kl (p-units)
        yy, wy = np.polynomial.legendre.leggauss(200); Y = 60/kl; yy = (yy + 1)/2*Y; wy = wy*Y/2
        fpl = np.array([complex(AF.flatplus(J, m, L*(1 + 1j*y))) for y in yy])
        Eflat_contour = float(2*np.real(-1j*np.sum(K.psi_c(2, 1 + 1j*yy)*fpl*wy)))
        rows = []
        for pe, sg in ((0.80, 0.03), (0.85, 0.03), (0.85, 0.05), (0.90, 0.02), (0.75, 0.05)):
            chi = erfc((pe - p)/sg)/2
            if pe - 5*sg < 0.55:
                continue
            Ef = float(np.sum(wp*ps*chi*flat)); Eg = float(np.sum(wp*ps*chi*g)); EG = float(np.sum(wp*ps*chi*G))
            Ec = float(np.sum(wp*ps*chi*cur))
            rows.append(dict(pe=pe, sigma=sg, kappa_sigma=kl*sg*np.sqrt(1 - (pe*L/m)**2)/np.sqrt(1 - x1*x1), Eflat=Ef, Eg=Eg, EG=EG, Ecur=Ec))
        # P4W's edge term (k = 2, no ring weight): 2[V(1;kl,T1=15)(e_sc(1) + A) + V(1;k0) A]
        esc1 = float(K.e_sc(beta, mfl)(np.array([1.0]))[0]); A = K.A_curv(beta, J, mfl)
        V1 = K.V(2, 1.0, kl, T, kap0=k0, T1=min(T, 15/kl)); V0 = K.V(2, 1.0, k0, T)
        edge_P4W = 2*(V1*(esc1 + A) + V0*A)
        near = p > 1 - 2*np.pi/kl
        o = dict(J=J, mfl=mfl, x1=x1, beta=beta, kl=kl, Eflat_contour=Eflat_contour, rows=rows, esc1=esc1, A=A, V1=V1, V0=V0,
                 edge_P4W=edge_P4W, gmax_near=float(np.max(np.abs(g[near]))), Gmax_near=float(np.max(np.abs(G[near]))),
                 flat1=float(abs(AF.flatplus(J, m, L))), t=time.time() - t0)
        res[J] = o
        pr(f"\n== J = {J}, m/L = {mfl:.4f}, x1 = {x1:.4f}, beta = {beta:.4f}, kappa_loc = {kl:.1f} (p-units); {len(p)} nodes  [{o['t']:.0f}s]")
        pr(f"   flat: right-side contour integral (Lemma F, exact) = {Eflat_contour:+.4e};  |flat_+(1)| = {o['flat1']:.3e}")
        pr(f"   near p = 1 (last period): max|g| = {o['gmax_near']:.3e}, max|G| = {o['Gmax_near']:.3e};  e_sc(1) = {esc1:.3e}, A = {A:.2e}")
        pr(f"   P4W edge term (Lam_AG = 1, k = 2, per wave pair, no ring weight) = {edge_P4W:.3e}  [V(1;kl) = {V1:.3e}, V(1;k0) = {V0:.3e}]")
        for r in rows:
            pr(f"   cutoff p_e = {r['pe']:.2f}, sigma = {r['sigma']:.2f} (kappa sigma ~ {r['kappa_sigma']:.0f}): E[flat] = {r['Eflat']:+.4e}, "
               f"E[g] = {r['Eg']:+.4e}, E[G] = {r['EG']:+.4e}, E[P-J0] = {r['Ecur']:+.2e};  |E[G]|/P4W-edge = {abs(r['EG'])/edge_P4W:.2e}")
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(AE.CACHE, 'edge.json'), 'w'), default=float)
    open(os.path.join(AE.CACHE, 'edge.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
