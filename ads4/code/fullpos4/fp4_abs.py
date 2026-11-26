"""
fullpos4/fp4_abs.py (FULLPOS4.md §3) [DV]: norm of the frequency-integrated outgoing §55 contour at a physical state.
  Weight psi2 (the main part G_A^psi of H/w^3), bare p = 3 kernel:  L(G^psi) = (C/pi) Re(I_U - i I_T - i I_L),
  W = xi^{s-1} (1-xi^2)^{-3} P+_h(beta) e^{i S al - i a ga}/sin(al),   P+_h from fp4_int.Amp.
Measured: N_abs = (C/pi) (int|W_U| + int|W_T| + int|W_L|), split by |beta| bands, against
  |G(c)| = (1-c)^{-3} |2 OmLR(beta_c)| (= L(G) to ~1e-6, FULLPOS3 §3.1) and N_s (T0 hyb_fast);  reported as w^3 N_abs / N_s.
Usage: fp4_abs.py a J b [kappa]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import fp4_int as I
from fp4_pf import Cconst


def gl(a, b, n=8):
    x, w = np.polynomial.legendre.leggauss(n)
    return a + (b - a)*(x + 1)/2, (b - a)/2*w


def graded(a, b, h0, r=1.3, n=8):
    e = [a]; h = h0
    while e[-1] + h < b:
        e.append(e[-1] + h); h *= r
    e.append(b)
    X, Wt = zip(*[gl(lo, hi, n) for lo, hi in zip(e[:-1], e[1:])])
    return np.concatenate(X), np.concatenate(Wt)


def geom(al, ga):
    x = np.exp(1j*al); y = np.exp(1j*ga)
    xi = (y - x)/(x*y - 1)
    beta = np.log((1 + xi)/(1 - xi))
    return xi, beta


MAJ = os.environ.get('FP4_MAJ', '0') == '1'
BANDS = [0, 0.01, 0.1, 0.3, 1, 3, 1e9]


def run(a, J, b, amp=None, p=3, kappa=None, toy=False):
    """kappa: retain only |A| <= A0 = kappa*Lam/S (hard localisation; the lifted remainder is not included)."""
    t0 = time.time()
    amp = amp or I.Amp(2)
    n = J - b; N = n + a; s = 2*b + 1; S = 2*J + 1 + a; eps = 1.0/S
    Delta = N + b + 2; T = (Delta - 1)*(Delta - 2); C = J*(J + 1); c = (C - b*(b + 1))/(T - b*(b + 1)); w = 1 - C/T
    Cc = Cconst(n, N, s)
    res = {}
    band_abs = np.zeros(len(BANDS) - 1)
    A0 = math.pi/2 if kappa is None else min(math.pi/2, kappa*I.LAM/S)
    A, wA = graded(0.0, A0, eps/20)
    tot = {}; sig = {}; totm = {}
    for piece in ('U', 'T', 'L'):
        ALs, GAs, WTs = [], [], []
        for Ai, wa in zip(A, wA):
            if piece == 'U':
                al = Ai - 1j*eps
                d, wd = graded(0.0, math.pi/2 - Ai, max(abs(al)*1e-5, 1e-12))
                ga = Ai + d - 1j*eps
            elif piece == 'T':
                al = Ai - 1j*eps
                d, wd = graded(0.0, 60.0/a + 3.0, 1e-4)
                ga = math.pi/2 - 1j*(eps + d)
            else:
                al = -Ai - 1j*eps
                d, wd = graded(0.0, 60.0/a, max(abs(al)*1e-5, 1e-12))
                ga = al - 1j*d
            ALs.append(np.full(len(d), al)); GAs.append(ga); WTs.append(wa*wd)
        al = np.concatenate(ALs); ga = np.concatenate(GAs); wt = np.concatenate(WTs)
        xi, beta = geom(al, ga)
        ok = np.abs(xi) > 1e-13
        Pr = np.zeros(len(al), complex); Et = np.zeros(len(al), complex)
        Pr[ok], Et[ok] = amp(beta[ok], factored=True)
        car = 1j*S*al - 1j*a*ga
        pre = xi**(s - 1)*(1 - xi*xi)**(-p)/np.sin(al)
        Wv = pre*(Pr*np.exp(car) + Et*np.exp(car + 1j*amp.lam*beta))      # e^{i Lam beta} combined with the carrier (no overflow)
        ab = np.abs(Wv)*wt
        if MAJ:
            Mr = np.zeros(len(al)); ME = np.zeros(len(al))
            Mr[ok], ME[ok] = amp.majorant(beta[ok], factored=True)
            abm = np.abs(pre)*(np.exp(car.real)*Mr + np.exp(car.real - amp.lam*beta.imag)*ME)*wt
            totm[piece] = float(abm.sum())
        tot[piece] = float(ab.sum()); sig[piece] = complex((Wv*wt).sum())
        bb = np.abs(beta)
        for i in range(len(BANDS) - 1):
            band_abs[i] += ab[(bb >= BANDS[i]) & (bb < BANDS[i + 1])].sum()
    Nabs = Cc/math.pi*sum(tot.values())
    Hsig = Cc/math.pi*(sig['U'] - 1j*sig['T'] - 1j*sig['L']).real
    if toy:
        return dict(Nabs=Nabs, signed=Hsig, pieces=sig)
    import ls_core as LSC
    beta_c = math.acosh((1 + c)/(1 - c))
    Gc = (1 - c)**(-3)*2*LSC.OmLR(2, [beta_c])[0]
    m = math.sqrt(T - C); hy = LSC.hyb_fast(1000, J, m/I.LAM, '++' if b == 0 else '+-')
    Ns = hy['Ns']
    out = dict(a=a, J=J, b=b, S=S, c=c, L=-math.log(1 - c), beta_c=beta_c, C=Cc, Gc=Gc, Ns=Ns, w3=w**3,
               Nabs=Nabs, Nmaj=(Cc/math.pi*sum(totm.values()) if MAJ else None), signed_quad=Hsig, signed_over_Gc=Hsig/Gc, pieces={k: Cc/math.pi*v for k, v in tot.items()},
               bands={f'{BANDS[i]}-{BANDS[i+1]}': Cc/math.pi*band_abs[i] for i in range(len(BANDS) - 1)},
               kappa=kappa, A0=A0, ratio_G=Nabs/abs(Gc), ratio_Ns=w**3*Nabs/Ns, w3Gc_over_Ns=w**3*Gc/Ns, npts=len(A), secs=time.time() - t0)
    return out


if __name__ == '__main__':
    a, J, b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    kap = float(sys.argv[4]) if len(sys.argv) > 4 else None
    r = run(a, J, b, kappa=kap)
    print(json.dumps(r), flush=True)
    open(os.path.join(HERE, 'cache', 'fp4_abs.jsonl'), 'a').write(json.dumps(r) + '\n')
