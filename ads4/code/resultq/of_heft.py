"""
of_heft (RESULT_Q_FERMION_HEFT.md, task 2): scalar-Mellin toy check of H_EFT for a single heavy exchange (Delta_h, J_h).

Setting: identical external scalars Df = 3 (= Delta_T), d = 3, k = 4; CHMRS 2106.10274 Mellin sum rule (eq. BMellinDefinition)
    Bhat_{4,T}[M] = oint_{|S'| = oo} dS'/(2 pi i) M(S', T') Kr(S'),   T' = T + 2Df - S',
    Kr(S') = 16 / [(2Df - S')^2 (2Df + 2 - S') (S' - T) (S' - T + 2)],
and its pure-moment (de-shifted) form  C_{4,nu}[X] ~ int_{Re T = c0} dT/(4 pi i) Gamma(T/2)^2 gamma_{A-T}(nu) X(T),  A = 2Df + 4 - 3/2 (eq.
CfromBMellin; at nu = 0 this is the de-shifted value functional of Result Q in the scalar proxy, RESULT_Q_STAGE4 §0 item 6).

Exchange Witten diagram (three channels; residues of the block, mack.py conventions, common constants dropped):
    M_h(S, T) = sum_n [ rho_n(T)/(S - s_n) + rho_n(S)/(T - s_n) + rho_n(T)/(U - s_n) ],  s_n = tau + 2n,  U = 4Df - S - T,
    rho_n(X) = K_n Qpoly^n(y = X - 2Df + tau/2 + n)   (exactmom.BlockMellin nodes/weights: direct n < N0 + Euler-Maclaurin tail).
Contour bookkeeping at fixed T (Cauchy; all residues of M_h Kr):
    H(T)   = residues at the s- and t-channel poles of M_h        = Bhat_{4,T}[heavy block]   (checked against exactmom.Bhat)
    L(T)   = residues at the kernel poles S' = 2Df (double), 2Df + 2, T, T - 2 = the s-channel DOUBLE-TRACE part (OPE side, exact)
    Arc(T) = H + L = the contour at infinity = the exchange's actual total contribution (CHMRS definition of Bhat).
    Lemma A analogue: Arc = 0 for J_h < 4 (M_h K ~ S'^{J_h - 5}).
Term-wise (P): the low-energy (1/s_n) expansion of M_h truncated at order R,
    M^(R) = sum_n [ rho_n(T)(1 - (S/s_n)^{R+1})/(S - s_n) + (t) + (u) ]   (a crossing-symmetric POLYNOMIAL for every R),
    each polynomial's light value is its kernel-pole sum L[M^(R)] (= its residue at infinity = CHMRS (P)); TP_R = C_{4,nu}[L[M^(R)]].
Usage: OMP_NUM_THREADS=1 python3 of_heft.py  -> cache/of_heft.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
from scipy.special import loggamma
import exactmom as E

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
DF = 3.0
K = 4
A = 2*DF + K - 1.5
OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


class Exchange:
    def __init__(self, Delta, J, N0=400):
        self.B = E.BlockMellin(Delta, J, DF, N0=N0)
        self.Delta, self.J, self.tau = float(Delta), int(J), float(Delta) - J
        lr = self.B.logK0 + 2*self.B.g1
        self.Lr = lr.max()
        self.wr = np.exp(lr - self.Lr)*self.B.wn          # K_n e^{-Lr} x (sum weights)
        self.sn = self.tau + 2*self.B.ns

    def rho(self, X):
        """[len(X), nn] array rho_n(X) (without the sum weight)."""
        X = np.atleast_1d(np.asarray(X, complex)).ravel()
        y = X[:, None] - 2*DF + self.tau/2 + self.B.ns[None, :]
        return self.B.Qeval(y)

    def M(self, S, T, R=None):
        """M_h(S, T) (R = None) or its truncation M^(R) (R: int or 1-d int array -> extra trailing axis)."""
        S = np.atleast_1d(np.asarray(S, complex)).ravel(); T = np.atleast_1d(np.asarray(T, complex)).ravel()
        U = 4*DF - S - T
        sn = self.sn[None, :]
        out = 0
        for (x, other) in ((S, T), (T, S), (U, T)):
            rq = self.rho(other)*self.wr[None, :]
            den = x[:, None] - sn
            if R is None:
                out = out + np.sum(rq/den, axis=1)
            else:
                Rs = np.atleast_1d(R)
                q = (x[:, None]/sn)                                   # [np, nn]
                fac = 1 - q[:, :, None]**(Rs[None, None, :] + 1)
                out = out + np.einsum('pn,pnr->pr', rq/den, fac)
        return out

    def rho_at(self, X):
        """rho_n(X_n) for a per-node argument X[n] (diagonal evaluation without the [N, N] matrix)."""
        X = np.asarray(X, complex)
        y = (X - 2*DF + self.tau/2 + self.B.ns)[None, :]
        return self.B.Qeval(y)[0]


def kr(Sp, T):
    return 16.0/((2*DF - Sp)**2*(2*DF + 2 - Sp)*(Sp - T)*(Sp - T + 2))


def Hfast(X, T):
    """H(T) for an array T without N^2 work."""
    out = np.zeros(T.size, complex)
    for i, t in enumerate(T):
        Sp = X.sn; Tp = t + 2*DF - Sp
        out[i] = np.sum(X.rho_at(Tp)*X.wr*kr(Sp, t))
        St = t + 2*DF - X.sn
        out[i] += np.sum(-X.rho_at(St)*X.wr*kr(St, t))
    return out


NC = int(os.environ.get('OF_NC', '40'))      # trapezoid nodes per kernel circle
TH = 2*np.pi*np.arange(NC)/NC


def Lkernel(X, T, R=None, rad=0.8):
    """kernel-pole residues of M (or M^(R)) Kr at S' = 2Df, 2Df+2, T, T-2: trapezoid on circles of radius rad (exponentially accurate)."""
    T = np.atleast_1d(np.asarray(T, complex)).ravel()
    tot = 0
    for c in (lambda t: 2*DF + 0*t, lambda t: 2*DF + 2 + 0*t, lambda t: t, lambda t: t - 2):
        cen = c(T)
        Sp = cen[:, None] + rad*np.exp(1j*TH)[None, :]               # [nT, NC]
        Tp = T[:, None] + 2*DF - Sp
        Mv = X.M(Sp.ravel(), Tp.ravel(), R=R)
        Kv = kr(Sp, T[:, None]).ravel()
        dS = (1j*rad*np.exp(1j*TH))[None, :]*np.ones_like(Sp)
        if R is None:
            vals = (Mv*Kv*dS.ravel()).reshape(Sp.shape)
            tot = tot + np.sum(vals, axis=1)/NC/(1j)            # (1/2pi i) oint = (1/NC) sum f dS/dtheta / i ... dtheta = 2pi/NC
        else:
            vals = (Mv*(Kv*dS.ravel())[:, None]).reshape(Sp.shape + (np.size(R),))
            tot = tot + np.sum(vals, axis=1)/NC/(1j)
    return tot


def Lchunk(X, T, R=None, chunk=6):
    """Lkernel over a long T array in chunks (memory: [chunk*4*NC, nn, len(R)] complex)."""
    outs = [Lkernel(X, T[i:i + chunk], R=R) for i in range(0, T.size, chunk)]
    return np.concatenate(outs, axis=0)


def tgrid(nu=0.0, c0=2.0, h=0.1, ymax=30.0):
    """nodes T = c0 + i y (y >= 0) and weights for C_{4,nu}[X] = 2 Re sum w X(T) (conjugate symmetry), normalised so that the weight
    integrates Gamma(T/2)^2 gamma_{A-T}(nu)/gamma_{A}(nu)^2 dT/(4 pi i)."""
    y = np.arange(0.0, ymax + h/2, h)
    T = c0 + 1j*y
    lw = (2*loggamma(T/2) + loggamma((A - T + 1j*nu)/2) + loggamma((A - T - 1j*nu)/2) - 2*loggamma((A + 1j*nu)/2)
          - 2*loggamma((A - 1j*nu)/2))
    w = np.exp(lw)*h/(4*np.pi)
    w[0] *= 0.5
    return T, w


def cint(Tw, vals):
    T, w = Tw
    return 2*np.real(np.sum(w*vals)) if np.ndim(vals) == 1 else 2*np.real(np.sum(w[:, None]*vals, axis=0))


def main():
    pr('# of_heft: scalar-Mellin toy check of H_EFT (single heavy exchange, Df = 3, d = 3, k = 4); see docstring for definitions')
    # ---------------------------------------------------------------- (0) validation
    pr('(0) VALIDATION')
    X = Exchange(20.0, 2)
    Ts = np.array([2.0 + 0.5j, 2.0 + 3j, 1.0 + 7j, 4.0 - 2j])
    H = Hfast(X, Ts)
    Lb, Bv = X.B.Bhat(4, Ts)
    rat = H/(Bv*np.exp(Lb - X.Lr))
    pr(f'  (0a) H(T) (my residue sum) / exactmom.Bhat_4,T (block), (Delta, J) = (20, 2), 4 complex T: {np.round(rat, 12)}  (expected -2: Bhat has the -1/2 density factor)')
    for N0 in (200, 400, 800):
        Y = Exchange(20.0, 4, N0=N0)
        pr(f'  (0b) n-sum convergence, (20, 4), N0 = {N0}: M(6, 2+3i) = {complex(Y.M(6.0, 2 + 3j)[0]):.12e},  H(2+3i) = {complex(Hfast(Y, np.array([2 + 3j]))[0]):.12e}')
    for rad in (0.5, 0.8):
        pr(f'  (0c) kernel-circle radius {rad}: L(2+3i) for (20, 4) = {complex(Lkernel(Y, np.array([2 + 3j]), rad=rad)[0]):.13e}')
    # C_{4,nu} normalisation vs exactmom (nu -> 0 extrapolated exact heavy image), same states, ratio must be state-independent up to
    # the known (-1)^J / q_{Delta,J} m^8 conversion; we compare the ratio C[H]/exactmom * q / m^8 across states.
    pr('  (0d) C_{4,0}[H] on Re T = 2 vs exactmom.HeavyC (nu -> 0, cubic fit in nu^2 through nu = 0.1..0.4); ratio x (-1)^J q m^-8 e^{-Lr}:')
    Tw = tgrid(0.0)
    for (D, J) in ((20.0, 0), (20.0, 2), (50.0, 4), (30.0, 6)):
        Y = Exchange(D, J)
        CH = cint(Tw, Hfast(Y, Tw[0]))
        hc = E.HeavyC(D, J, DF)
        nus = np.array([0.1, 0.2, 0.3, 0.4]); v = hc.value(4, nus)
        c0 = np.linalg.solve(np.vander(nus**2, 4, increasing=True), v)[0]
        conv = (-1)**J*np.exp(hc.logq - Y.Lr)/hc.m2**4
        pr(f'      (Delta, J) = ({D}, {J}): C[H] = {CH:+.10e}, exactmom m^8 C_4,0 = {c0:+.10e}, normalised ratio = {CH/(c0*conv):+.10f}')
    pr(f'  [{time.time() - T0:.0f}s]')
    # ---------------------------------------------------------------- (1) fixed-T identity
    pr('(1) FIXED-T CONTOUR BOOKKEEPING  Arc(T) = H(T) + L(T)  (Lemma-A analogue: Arc = 0 iff J_h < 4); T = 2 + 3i, 1 + 7i:')
    for D in (20.0, 50.0):
        for J in (0, 2, 3, 4, 6):
            if J == 3:
                continue
            Y = Exchange(D, J)
            T = np.array([2 + 3j, 1 + 7j])
            Hh = Hfast(Y, T); Ll = Lkernel(Y, T)
            pr(f'    (Delta_h, J_h) = ({D:5.1f}, {J}): |Arc|/|H| = {abs(Hh[0] + Ll[0])/abs(Hh[0]):.1e}, {abs(Hh[1] + Ll[1])/abs(Hh[1]):.1e}   '
               f'(L/H = {(Ll[0]/Hh[0]).real:+.6f}{(Ll[0]/Hh[0]).imag:+.1e}i)')
    pr(f'  [{time.time() - T0:.0f}s]')
    # ---------------------------------------------------------------- (2) integrated, term-wise
    pr('(2) DE-SHIFTED VALUE FUNCTIONAL (pure moment C_{4,nu=0}; also nu = 2), Re T = 2:  H = heavy block,  D = actual double-trace value')
    pr('    (OPE side, kernel poles),  Arc = H + D (actual total),  TP_R = term-wise (P) of the order-R low-energy expansion of M_h;')
    pr('    "best" = optimal truncation (min over R of |TP_R - D|), R* its order; all relative to |H|.')
    Rs = np.array([0, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 60])
    for nu in (0.0, 2.0, 5.0):
        Tw = tgrid(nu)
        for D in (20.0, 50.0, 100.0):
            for J in (0, 2, 4, 6):
                Y = Exchange(D, J)
                Hh = cint(Tw, Hfast(Y, Tw[0]))
                Ll = cint(Tw, Lchunk(Y, Tw[0]))
                if nu != 0.0:
                    pr(f'    nu = {nu}  (Delta_h, J_h) = ({D:5.1f}, {J}):  D/H = {Ll/Hh:+.10f}   Arc/H = {(Hh + Ll)/Hh:+.3e}')
                    continue
                TP = cint(Tw, Lchunk(Y, Tw[0], R=Rs))
                err = np.abs(TP - Ll)/abs(Hh)
                rb = int(np.argmin(err))
                pr(f'    nu = {nu}  (Delta_h, J_h) = ({D:5.1f}, {J}):  H = {Hh:+.6e}  D/H = {Ll/Hh:+.10f}   Arc/H = {(Hh + Ll)/Hh:+.4e}   '
                   f'|TP_R - D|/|H| at R = ' + ', '.join(f'{r}: {e:.1e}' for r, e in zip(Rs, err)) + f';  best at R = {Rs[rb]}')
        pr(f'  [{time.time() - T0:.0f}s]')
    pr(f'# total {time.time() - T0:.0f} s')
    open(os.path.join(CACHE, 'of_heft.txt'), 'w').write('\n'.join(OUT) + '\n')


def moments():
    """(3) exact and absolute moments of the stored sc3 de-shift weights (functionals_q4_deshift.json): the term-wise (P) value of a
    light-side polynomial p(nu^2) = sum_n p_n nu^{2n} under the de-shifted functional is sum_n p_n m_n."""
    import json
    import sympy as sp
    d = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['sc3']
    nu_s = sp.Symbol('nu')
    Rj = {int(j): sp.lambdify(nu_s, sp.sympify(v), 'numpy') for j, v in d['R_j'].items()}
    dnu, N = d['dnu'], d['N']; nu = (np.arange(N) + 0.5)*dnu
    pr(f'(3) sc3 DE-SHIFT MOMENTS m_n = sum_j int chi R_j (nu - i j)^(2n) (grid dnu = {dnu}, N = {N}, nu_max = {nu[-1]:.2f}); |.|: absolute moments')
    for key in ('chi_value', 'chi_deriv'):
        chi = np.array(d[key]); row = []
        for n in range(10):
            m = sum(np.sum(dnu*chi*Rj[j](nu + 0j)*(nu - 1j*j)**(2*n)) for j in Rj)
            a = sum(np.sum(dnu*np.abs(chi*Rj[j](nu + 0j))*np.abs(nu - 1j*j)**(2*n)) for j in Rj)
            row.append(f'n={n}: {m.real:+.3e} (|.| {a:.2e})')
        pr(f'    {key}: ' + ';  '.join(row))
    open(os.path.join(CACHE, 'of_heft_moments.txt'), 'w').write('\n'.join(OUT) + '\n')


def diag():
    """(4) cancellation diagnostic of the T-integral: kappa = int |w| (|H| + |D|) / |C[H]|; float64 floor ~ 1e-16 kappa."""
    pr(f'    [NC = {NC} nodes per kernel circle]')
    pr('(4) T-INTEGRAL CANCELLATION (nu = 0): kappa = int|w|(|H|+|D|) / |C[H]|  -> expected floor of D/H + 1 and of |TP_R - D|/|H| ~ 1e-16 kappa')
    Tw = tgrid(0.0)
    for (D_, J) in ((20.0, 0), (20.0, 2), (50.0, 0), (50.0, 2), (100.0, 0), (100.0, 2)):
        Y = Exchange(D_, J)
        Hv = Hfast(Y, Tw[0]); Lv = Lchunk(Y, Tw[0])
        CH = cint(Tw, Hv)
        kap = 2*np.sum(np.abs(Tw[1])*(np.abs(Hv) + np.abs(Lv)))/abs(CH)
        pr(f'    (Delta_h, J_h) = ({D_}, {J}): kappa = {kap:.2e}  (floor ~ {1e-16*kap:.0e});  (H + D)/H = {(CH + cint(Tw, Lv))/CH:+.2e}')
    open(os.path.join(CACHE, f'of_heft_diag_nc{NC}.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    if sys.argv[1:] == ['diag']:
        diag()
    elif sys.argv[1:] == ['mom']:
        moments()
    else:
        main()
