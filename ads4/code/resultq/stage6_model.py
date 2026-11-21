"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §1): ONE heavy-side model for the whole physical region, and the pieces to cross-validate it.

Heavy state (Delta, J), twist >= Delta_g:  m^2 = (Delta-J-2)(Delta+J-1),  Lam^2 = (Dg-2)(Dg-1),  u = Lam^2/m^2 <= u_gap(J) = (Dg-1)/(Dg+2J-1),
bhat = (2J+1) sqrt(u),  cosh beta = 1 + (2J+1)^2/(2m^2) = 1 + bhat^2/(2 Lam^2)   <- beta is a function of bhat ALONE (given Lam).
Channels: '++' (J even), '+-' (J >= 4); the parity-odd even-J block has the ++ image at this order (Stage 3 §2.1).

Models (all with BOTH forward pieces de-shifted, i.e. pure forward moments; Stage 4/5):
 * LR  (Stage 12/5 leading Regge):  B2 = 2 I20(beta) - u I22(beta),  B3 = s u I3(beta),
        B4 = 2u^2 P0,  dB4 = -2u^2 Lam^2 D + c4 u^3 P0,  B5 = s u^3 P0,
        I20 = int psi2 P_{-1/2+i Lam p}(cosh beta) dp, I22 = int psi2 p^2 P.., I3 = int psi3 P..; P0 = P_{-1/2}, D = -d_{nu^2}P|_0.
 * HYB (this stage; = Stage 3's hybrid for the smeared parts, extended to the forward parts):  HYB = FLAT_EXACT + (LR - LR_flat),
        FLAT_EXACT: Stage-0 action at (u, J) with exact Wigner functions (common_q.parts) and flat forward normal forms
                    B4 = 2u^2, dB4 = (c4 - 2JJ) u^3, B5 = s u^3;
        LR_flat:    the beta -> 0 (fixed bhat) limit of LR: P_{-1/2+i Lam p} -> J0(bhat p), P0 -> 1, Lam^2 D -> bhat^2/4.
        So HYB = LR + [FLAT_EXACT - LR_flat]: the leading-Regge AdS curvature plus the finite-J (O(u bhat), O(u)) flat corrections.
        The curvature corrections Delta X(bhat) = X_AdS(beta(bhat)) - X_flat(bhat) depend on bhat only -> 1D tables + splines.
 * EXA (Stage 3, where exact CHMRS skeleton moments exist): HYB smeared parts + B4 = u^2 ex4(0), dB4 = (c4-2JJ)u^3 + u^2[Lam^2 d ex4 - u(5-2JJ)],
        B5 = s u^3 ex5(0)/2.
Heavy action FF = B2 + B3 - a (dB4 - B5) - c B4;  Stage-0 normalisation N = |B2| + |B3| + |a(dB4 - B5)| + c|B4|  (THE normalisation of the
LP margin delta); Stage-5 normalisation Nrm = |B2| + |B3| + a(2u^2 Lam^2 D + 42 u^3 P0 + u^3 P0) + 2c u^2 P0 (>= N), kept for comparison.
"""
import os, sys, json, glob, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.interpolate import CubicSpline
import common_q as CQ
import stage12_forward_regge as FR
import stage5_positivity as S5P

CACHE = os.path.join(HERE, 'cache')
C4 = {'++': 1.0, '+-': 41.0}; SIG = {'++': 1.0, '+-': -1.0}
FUNCS = {0.05: 'Q0-full-robust', 0.2: 'Q3-full-delta0.2'}


def lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def u_gap(J, Dg):
    return (Dg - 1.0)/(Dg + 2.0*np.asarray(J, float) - 1.0)


def jmax_gap(u, Dg):
    return np.floor(((Dg - 1.0)/u - Dg + 1.0)/2.0 + 1e-9)


def load(name):
    F, c2, c3 = S5P.load_func(name)
    return F, c2, c3


class Tables:
    """1D tables in bhat of the AdS (beta(bhat)) and flat (beta -> 0) leading-Regge quantities, for one functional and one Lam."""

    def __init__(self, name, Lam, bmax, verbose=False):
        fn = os.path.join(CACHE, 'stage6_tab_%s_%d_%d.pkl' % (name.split()[0], round(Lam), round(bmax)))
        if os.path.exists(fn):
            self.__dict__.update(pickle.load(open(fn, 'rb'))); self._spl(); return
        F, c2, c3 = load(name)
        T = [S5P.CosT(c2), S5P.CosT(np.convolve(c2, [0, 0, 1])), S5P.CosT(c3)]
        b = np.unique(np.concatenate([np.linspace(0, 60, 3001), np.linspace(60, 400, 3401), np.arange(400, bmax + 1.0, 0.5)]))
        out = np.zeros((10, len(b)))
        for i, bh in enumerate(b):
            if bh == 0:
                f = [T[0](np.array([0.]))[0], T[1](np.array([0.]))[0], T[2](np.array([0.]))[0]]
                out[:, i] = f + [1., 0.] + f + [0., 1.]
                continue
            beta = 2*np.arcsinh(bh/(2*Lam))
            I20, I22, I3, P0, D, H5 = S5P.beta_quantities(beta, Lam, T)
            n = int(min(2.4e5, 800 + 5*bh)); n = max(800, (n//20)*20)
            ph, wph = FR.composite_gl(0., np.pi/2, n//20, 20)
            om = bh*np.sin(ph); W = (2/np.pi)*wph
            out[:, i] = [I20, I22, I3, P0, Lam**2*D, (W*T[0](om)).sum(), (W*T[1](om)).sum(), (W*T[2](om)).sum(), bh*bh/4, H5]
            if verbose and i % 2000 == 0:
                print('  tab', name, Lam, bh, flush=True)
        self.b = b; self.v = out
        pickle.dump(dict(b=b, v=out), open(fn, 'wb'))
        self._spl()

    def _spl(self):
        I20, I22, I3, P0, LD, J20, J22, J3, LDf, H5 = self.v
        self.s = {'dI20': CubicSpline(self.b, I20 - J20), 'dI22': CubicSpline(self.b, I22 - J22), 'dI3': CubicSpline(self.b, I3 - J3),
                  'dP0': CubicSpline(self.b, P0 - 1), 'dLD': CubicSpline(self.b, LD - LDf),
                  'J20': CubicSpline(self.b, J20), 'J22': CubicSpline(self.b, J22), 'J3': CubicSpline(self.b, J3)}

    def __call__(self, key, bh):
        return self.s[key](bh)


def hyb(F, h, u, J, B2f, B3f, tab, Lam, fw=None):
    """HYB (or EXA if fw = (ex4, Lam^2 d ex4, ex5) given) heavy action. Arrays broadcast. Returns FF, N, Nrm, parts dict."""
    a, c = F['a'], F['c']; s = SIG[h]; c4 = C4[h]
    bh = (2*J + 1)*np.sqrt(u); JJ = J*(J + 1.0)
    dI20, dI22, dI3 = tab('dI20', bh), tab('dI22', bh), tab('dI3', bh)
    dP0, dLD = tab('dP0', bh), tab('dLD', bh)
    P0 = 1 + dP0; LD = bh*bh/4 + dLD
    B2 = B2f + 2*dI20 - u*dI22
    B3 = B3f + s*u*dI3
    if fw is None:
        B4 = 2*u*u*P0
        dB4 = (c4 - 2*JJ)*u**3 - 2*u*u*dLD + c4*u**3*dP0
        B5 = s*u**3*P0
    else:
        ex4, d4L, ex5 = fw
        B4 = u*u*ex4
        dB4 = (c4 - 2*JJ)*u**3 + u*u*(d4L - u*(5 - 2*JJ))
        B5 = s*u**3*ex5/2
    FW = -a*(dB4 - B5)
    FF = B2 + B3 + FW - c*B4
    N = np.abs(B2) + np.abs(B3) + np.abs(FW) + c*np.abs(B4)
    Nrm = np.abs(B2) + np.abs(B3) + a*(2*u*u*LD + 42*u**3*P0 + u**3*P0) + 2*c*u*u*P0
    return FF, N, Nrm, dict(B2=B2, B3=B3, FW=FW, B4=B4, bh=bh)


def lr(F, h, u, J, tab, Lam):
    """LR model on the same states (J may be non-integer). Uses the AdS tables directly (I = J_flat + dI)."""
    a, c = F['a'], F['c']; s = SIG[h]; c4 = C4[h]
    bh = (2*J + 1)*np.sqrt(u)
    I20 = tab('J20', bh) + tab('dI20', bh); I22 = tab('J22', bh) + tab('dI22', bh); I3 = tab('J3', bh) + tab('dI3', bh)
    P0 = 1 + tab('dP0', bh); LD = bh*bh/4 + tab('dLD', bh)
    B2 = 2*I20 - u*I22; B3 = s*u*I3
    B4 = 2*u*u*P0; dB4 = -2*u*u*LD + c4*u**3*P0; B5 = s*u**3*P0
    FW = -a*(dB4 - B5)
    FF = B2 + B3 + FW - c*B4
    N = np.abs(B2) + np.abs(B3) + np.abs(FW) + c*np.abs(B4)
    Nrm = np.abs(B2) + np.abs(B3) + a*(2*u*u*LD + 42*u**3*P0 + u**3*P0) + 2*c*u*u*P0
    return FF, N, Nrm


def flat_table(name, xs, Jcap, Dgmax=1000):
    """flat exact B2f, B3f for every u = 1/x^2 in xs and J <= min(Jcap, jmax_gap(u, Dgmax)); cached. Returns dict h -> list of (u, Js, B2, B3)."""
    fn = os.path.join(CACHE, 'stage6_flat_%s_%d_%d.pkl' % (name.split()[0], len(xs), Jcap))
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    F, c2, c3 = load(name)
    out = {'++': [], '+-': []}
    for x in xs:
        u = 1.0/x**2
        Jm = int(min(Jcap, jmax_gap(u, Dgmax)))
        for h in ('++', '+-'):
            if Jm < (0 if h == '++' else 4):
                out[h].append((u, np.zeros(0, int), np.zeros(0), np.zeros(0))); continue
            Js, B2, B3 = CQ.parts(h, u, max(Jm, 4), [c2], [c3])
            sel = Js <= Jm
            if h == '++':
                sel &= Js % 2 == 0
            out[h].append((u, Js[sel], B2[sel, 0], B3[sel, 0]))
    pickle.dump(out, open(fn, 'wb'))
    return out


def gap_states(name, Dg, Jtop):
    """flat exact parts at the gap state u = u_gap(J) of every J <= Jtop (both channels)."""
    F, c2, c3 = load(name)
    out = {'++': [], '+-': []}
    for J in range(0, Jtop + 1):
        u = float(u_gap(J, Dg))
        for h in ('++', '+-'):
            if (h == '++' and J % 2) or (h == '+-' and J < 4):
                continue
            Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3])
            k = list(Js).index(J)
            out[h].append((u, J, B2[k, 0], B3[k, 0]))
    return out


def load_exact(Dg, Df=3):
    db = {}
    for fn in glob.glob(os.path.join(CACHE, 'stage3_fw_%d_Df%g_w*.pkl' % (Dg, Df))) + \
            glob.glob(os.path.join(CACHE, 'stage6_fw_%d_Df%g_w*.pkl' % (Dg, Df))):
        db.update(pickle.load(open(fn, 'rb')))
    return db
