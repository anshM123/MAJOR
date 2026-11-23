"""agproof/ag_flat.py (AG_PROOF.md §1–2): Lemma F (the flat partial wave as an exact pair of analytic travelling waves, with an explicit
decay bound) and Proposition N (the AG' majorant against the flat wave's own decay at the edge nu = Lam).

Lemma F.  For integer n >= 0 and 0 < Re theta < pi:
    a_k = (1/2)_k (n+1)_k / (k! (n+3/2)_k) > 0,   F_n(q) = sum_k a_k q^k = 2F1(1/2, n+1; n+3/2; q),
    C_n = sqrt(pi) Gamma(n+1) / (2 Gamma(n+3/2)),
    W+(theta) = (2/(pi i)) C_n e^{i(n+1)theta} F_n(e^{2 i theta}),   W-(theta) = -(2/(pi i)) C_n e^{-i(n+1)theta} F_n(e^{-2 i theta}),
    P_n(cos theta) = W+(theta) + W-(theta)      (DLMF 14.13.1, mu = 0),
    |W+(theta)| <= (2/pi) C_n e^{-(n+1) eta} (1 - e^{-2 eta})^{-1/2},   eta = Im theta > 0      (a_k <= (1/2)_k/k!).
With theta(nu) = 2 arcsin(nu/m) the flat image flat(nu) = E_2(x) P_J(1 - 2x^2)/2 (x = nu/m) splits as flat_+ + flat_-,
flat_pm = E_2 W_pm(theta(nu))/2, analytic in the upper / lower half nu-plane for |Re nu| < m.

This script checks the identity and the bound (mpmath, 30 digits), then compares, on the right side p = 1 of Lemma P4W's rectangle, the
flat wave's own decay |flat_+(Lam + i y)| / |flat_+(Lam)| with the AG' majorant exp(-kappa_AG'(y)) (local rate 0.97 Lam beta/sqrt(1-x1^2)
up to y = 15/kappa, then 0.97 Lam beta), state by state (Theorem W's 47 states, the 4 exact gap states, Theorem L''s states):
  * max_y of the ratio (violation factor), and the edge integrals
    V_F  = int_0^T |psi_2(1 + i y)| |flat_+(1+iy)|/|flat_+(1)| dy      (true flat profile, amplitude included)
    V_Fp = int_0^T |psi_2(1 + i y)| exp(-(J+1/2) Im[theta(1+iy) - theta(1)]) dy   (phase-only flat profile)
    V_AG = int_0^T |psi_2(1 + i y)| exp(-kappa_AG'(y)) dy                (P4_WINDOW's majorant, T1 = 15)
  (k = 2; k = 3 in brackets).  -> cache/flat.txt, cache/flat.json
"""
import os, sys, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
CACHE = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(CODE, 'p4win'))
import numpy as np
import mpmath as mp
import p4w_common as C            # read-only
import p4w_contour as K           # read-only

mp.mp.dps = 30
L = C.L
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def Cn(n):
    return mp.sqrt(mp.pi)*mp.gamma(n + 1)/(2*mp.gamma(n + mp.mpf(1.5)))


def Wplus(n, th):
    th = mp.mpc(th)
    q = mp.exp(2j*th)
    return 2/(mp.pi*1j)*Cn(n)*mp.exp(1j*(n + 1)*th)*mp.hyp2f1(0.5, n + 1, n + 1.5, q)


def Wminus(n, th):
    th = mp.mpc(th)
    q = mp.exp(-2j*th)
    return -2/(mp.pi*1j)*Cn(n)*mp.exp(-1j*(n + 1)*th)*mp.hyp2f1(0.5, n + 1, n + 1.5, q)


def check_lemmaF():
    pr('Lemma F checks (mpmath, 30 digits):')
    worst_id = 0; worst_bd = 0
    for n in (0, 1, 7, 26, 40, 100, 344):
        for th in (0.3, 1.1, 2.0, 2.9, mp.mpc(1.3, 0.2), mp.mpc(2.6, 0.05), mp.mpc(0.4, 1.0)):
            P = mp.legendre(n, mp.cos(th))
            s = Wplus(n, th) + Wminus(n, th)
            worst_id = max(worst_id, float(abs(s - P)/max(abs(P), mp.mpf(1e-30))))
            eta = mp.im(th)
            if eta > 0:
                bd = 2/mp.pi*Cn(n)*mp.exp(-(n + 1)*eta)/mp.sqrt(1 - mp.exp(-2*eta))
                worst_bd = max(worst_bd, float(abs(Wplus(n, th))/bd))
    pr(f'  identity P_n(cos th) = W+ + W-: max rel. error {worst_id:.1e} (n = 0..344, real and complex th)')
    pr(f'  bound |W+| / [(2/pi) C_n e^(-(n+1) eta) (1-e^(-2 eta))^(-1/2)] <= {worst_bd:.3f} (must be <= 1)')
    return worst_id, worst_bd


def flatplus(J, m, nu):
    x = mp.mpc(nu)/m
    th = 2*mp.asin(x)
    E2 = (2 - x*x)/(1 - x*x)**2
    return E2*Wplus(J, th)/2


def edge_compare(J, mfl, ny=161):
    m = mfl*L; beta = C.beta_of(J, mfl); x1 = 1/mfl
    k0 = 0.97*beta*L; kl = k0/np.sqrt(1 - x1*x1)            # p-units
    T = 25.0/k0
    T1 = min(T, 15.0/kl)
    # Gauss-Legendre on [0, T] split at T1 and at a few multiples of 1/kl (integrand peaks near y ~ 5/kl)
    br = sorted(set([0.0] + [min(T, c/kl) for c in (2, 5, 10, 15, 30, 60)] + [T]))
    xs, ws = np.polynomial.legendre.leggauss(24)
    ys = []; wy = []
    for a, b in zip(br[:-1], br[1:]):
        if b > a:
            ys += list(a + (b - a)*(xs + 1)/2); wy += list(ws*(b - a)/2)
    ys = np.array(ys); wy = np.array(wy)
    f0 = flatplus(J, m, L)
    fl = np.array([float(abs(flatplus(J, m, L*(1 + 1j*y)))) for y in ys])/float(abs(f0))
    th = lambda y: complex(2*mp.asin(mp.mpc(1, y)*L/m))
    ph = np.array([(J + 0.5)*(th(y).imag - th(0).imag) for y in ys])
    maj = np.where(ys <= T1, kl*ys, kl*T1 + k0*(ys - T1))
    out = dict(J=J, mfl=mfl, x1=x1, beta=beta, T=T, kl=kl)
    for k in (2, 3):
        ps = np.abs(K.psi_c(k, 1 + 1j*ys))
        out[f'VF{k}'] = float((ps*fl*wy).sum())
        out[f'VFp{k}'] = float((ps*np.exp(-ph)*wy).sum())
        out[f'VAG{k}'] = float((ps*np.exp(-maj)*wy).sum())
    rat = fl/np.exp(-maj)
    out['maxratio'] = float(rat.max()); out['y_at_max_kl'] = float(ys[rat.argmax()]*kl)
    rat_p = np.exp(-ph)/np.exp(-maj)
    out['maxratio_phase'] = float(rat_p.max())
    # profile samples at y = c/kl
    out['prof'] = {c: (float(np.interp(c/kl, ys, -np.log(fl))), float(np.interp(c/kl, ys, ph)), float(np.interp(c/kl, ys, maj)))
                   for c in (2, 5, 10, 15) if c/kl <= T}
    return out


def states_W():
    st = json.load(open(os.path.join(CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']
    return sorted([(s['J'], s['mfl']) for s in st], key=lambda v: (v[0] % 2, v[0], v[1]))


def states_exact():
    out = []
    for J in (40, 48, 64, 100):
        al = 998.0; ab = al + 2*J + 1
        out.append((J, np.sqrt(al*ab)/L))
    return out


def states_L():
    sys.path.insert(0, os.path.join(CODE, 'p4win'))
    import p4w_L as PL                                             # read-only
    return [(J, PL.mfl_exact(lab, J, mt)) for lab, J, mt, _, _ in PL.ROWS]


def main():
    check_lemmaF()
    res = {}
    for name, sts in (('W47', states_W()), ('exact4', states_exact()), ('Lprime', states_L())):
        pr(f'\n== {name}: right side p = 1; profiles at y = c/kappa_loc (c = 2, 5, 10, 15): [true flat -log ratio | phase-only | AG\' majorant]; '
           f'edge integrals V (k = 2) true-flat / phase-only / AG\' and their ratios; max_y |flat_+|-ratio / AG\' majorant')
        rows = []
        for J, mfl in sts:
            o = edge_compare(J, mfl)
            rows.append(o)
            prof = ' '.join(f'{c}:[{v[0]:.1f}|{v[1]:.1f}|{v[2]:.1f}]' for c, v in o['prof'].items())
            pr(f"  J={J:7d} m/L={mfl:.3f} x1={o['x1']:.3f} kl/k0={1/np.sqrt(1-o['x1']**2):.2f} | {prof} | V2: {o['VF2']:.2e} / {o['VFp2']:.2e} / "
               f"{o['VAG2']:.2e}  VF/VAG = {o['VF2']/o['VAG2']:.2f} (k=3: {o['VF3']/o['VAG3']:.2f}), VFp/VAG = {o['VFp2']/o['VAG2']:.2f} | "
               f"max ratio {o['maxratio']:.2g} at y = {o['y_at_max_kl']:.1f}/kappa (phase-only {o['maxratio_phase']:.2g})")
        res[name] = rows
    json.dump(res, open(os.path.join(CACHE, 'flat.json'), 'w'), default=float)
    open(os.path.join(CACHE, 'flat.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
