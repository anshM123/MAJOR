"""
tail.py (NU_TAIL_CHECK.md): how much of the published k = 2 de-shift weights' action on a near-gap heavy state comes from nu > Lambda
(nu > m), when the per-nu heavy image is the EXACT one (which grows above the state's mass) rather than the leading-Regge one.

Proxy for the EE raw row Gamma(nu) (= raw image/(R_- + R_+), ~ P_{-1/2+i nu}(cosh beta) ~ 1 at small nu): the exact scalar CHMRS action
Gamma(nu) := m^4 C_{2,nu}[Delta, J]/2 (Df = 3), computed contour-free in Arb (exactskel/skelexact3.py, read-only; N0 = 240, nq = 320,
converged to <= 1e-9 at nu = 1.31 m).  Its normalisation (-> 1 at nu -> 0) matches the rows'.  [The EE integrand is the scalar one
times a rational factor with no new divisor (RD_CLOSE §1); the e^{c nu} growth comes from the common 1/N_com ~ gamma_a^-2 ~ e^{pi nu}
normalisation, so the scalar growth profile is the natural proxy -- a model statement, [M] for EE.]

Weights (prl/data, read-only), each times the published cutoff theta(nu) = (erf((nu+nu_c)/a) - erf((nu-nu_c)/a))/2, nu_c = Lam+100, a = 10:
  Dgrid / Qgrid : d_abel_X_{Dg}.txt / k2/xq_abel_X_{Dg}.txt  (functional Re sum_k dnu X_k Gamma_k)
  Dcont / Qcont : k2/k2_continuum_{D,Q}_{Dg}.json, X = chi (R_- + R_+), chi evaluated from the JSON formula in 40-digit mpmath
                  on [Lam-5, Lam+45] (beyond Lam+45 the Gaussians are < e^-1000); below Lam-5 the published table.
  Dexact        : the grid weight up to Lam, continued beyond Lam by its own measured exponential tail (rate fitted on [Lam-2, Lam+3]),
                  i.e. the exact min-norm X without the double-precision floor.
For every (weight, state): A_psi = int_0^Lam psi Gamma (psi-form action, the form all per-state evidence uses),
A_X = Re sum_{nu<=Lam} dnu X theta Gamma, and the tails beyond Lam and beyond m (signed, and with |X|).
usage: python3 -B tail.py Dg [hi]   -> cache/tail_<Dg>[_hi].txt
  'hi': the higher-mass J = 0 states at Dg = 1000 (m = Lam+150, Lam+250), coarser tail sampling
"""
import os, sys, json, time, warnings
sys.dont_write_bytecode = True
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import mpmath as mp
from scipy.special import erf
import skelexact3 as S3

DATA = os.path.join(HERE, '..', '..', 'prl', 'data')
Dg = int(sys.argv[1])
HI = len(sys.argv) > 2 and sys.argv[2] == 'hi'
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
nuc, aw = Lam + 100.0, 10.0
theta = lambda nu: 0.5*(erf((nu + nuc)/aw) - erf((nu - nuc)/aw))
out = []
def log(s=''):
    print(s, flush=True); out.append(s)

# ---------------- kernels R_-+ (prl/data/kernels.json, D_k2_EE) ----------------
K = json.load(open(os.path.join(DATA, 'kernels.json')))['D_k2_EE']
def Rsum(nu):
    I = 1j
    ns = {'nu': nu, 'I': I}
    return eval(K['R_minus'].replace('nu', '(nu)'), {}, ns) + eval(K['R_plus'].replace('nu', '(nu)'), {}, ns)

# ---------------- weights ----------------
def load_grid(fn):
    a = np.loadtxt(fn, comments='#')
    return a[:, 0], a[:, 1] + 1j*a[:, 2], a[:, 5]

def cont_X(fn, nus):
    d = json.load(open(fn))
    mp.mp.dps = 40
    cs = d['centres_sigma']; xr = [mp.mpf(v) for v in d['x_re']]; xi = [mp.mpf(v) for v in d['x_im']]
    yr = [mp.mpf(v) for v in d['y_re']]; yi = [mp.mpf(v) for v in d['y_im']]
    res = np.zeros(len(nus), complex)
    for i, nu0 in enumerate(nus):
        nu = mp.mpf(float(nu0))
        Q = nu**2/(nu**2 + mp.mpf(1)/4)**2
        sr = mp.mpf(0); si = mp.mpf(0)
        for (c, s), a, b in zip(cs, xr, xi):
            if abs(float(nu0) - c) > 40*s and abs(float(nu0) + c) > 40*s:
                continue
            c = mp.mpf(c); s = mp.mpf(s)
            g = mp.exp(-(nu - c)**2/(2*s*s)) + (mp.exp(-(nu + c)**2/(2*s*s)) if c != 0 else 0)
            sr += a*g; si += b*g
        l = (nu**2 + 1)*mp.log(nu**2 + 1)*mp.exp(-nu**2/2)
        sr = Q*sr + sum(yr[j]*nu**(2*j) for j in range(3))*l
        si = Q*si + sum(yi[j]*nu**(2*j) for j in range(3))*l
        res[i] = complex(sr, si)*complex(Rsum(complex(nu0)))
    return res

W = {}
nuD, XD, psiD = load_grid(os.path.join(DATA, f'd_abel_X_{Dg}.txt'))
nuQ, XQ, psiQ = load_grid(os.path.join(DATA, 'k2', f'xq_abel_X_{Dg}.txt'))
dnu = nuD[1] - nuD[0]
W['Dgrid'] = (nuD, XD, psiD)
W['Qgrid'] = (nuQ, XQ, psiQ)
for tag in ('D', 'Q'):
    fn = os.path.join(DATA, 'k2', f'k2_continuum_{tag}_{Dg}.json')
    tb = np.loadtxt(os.path.join(DATA, 'k2', f'k2_continuum_{tag}_{Dg}_table.txt'), comments='#')
    nu_t, X_t, psi_t = tb[:, 0], tb[:, 3] + 1j*tb[:, 4], tb[:, 5]
    lo = nu_t < Lam - 5
    nu_hi = np.arange(nu_t[lo][-1] + dnu, Lam + 45, dnu)
    t0 = time.time()
    X_hi = cont_X(fn, nu_hi)
    log(f'# {tag}cont: evaluated chi in mpmath on [{nu_hi[0]:.2f}, {nu_hi[-1]:.2f}] ({len(nu_hi)} pts, {time.time()-t0:.0f}s)')
    W[tag + 'cont'] = (np.concatenate([nu_t[lo], nu_hi]), np.concatenate([X_t[lo], X_hi]),
                       np.concatenate([psi_t[lo], np.zeros(len(nu_hi))]))
# exact-tail continuation of D's grid weight: fit log|X| on [Lam-2, Lam+3] (before the floor)
sel = (nuD > Lam - 2) & (nuD < Lam + 3)
rate, c0 = np.polyfit(nuD[sel], np.log(np.abs(XD[sel])), 1)
Xex = XD.copy()
beyond = nuD > Lam
Xex[beyond] = np.exp(c0 + rate*nuD[beyond])        # magnitude only; sign irrelevant for the |.| bound
W['Dexact'] = (nuD, Xex, psiD)
log(f'# Dg = {Dg}, Lam = {Lam:.4f}; cutoff theta: nu_c = Lam+100 = {nuc:.2f}, a = {aw}; grid dnu = {dnu}')
log(f'# D grid weight: exact-tail decay rate beyond Lam (fit on [Lam-2, Lam+3]) = {-rate:.4f} (pi/2 = {np.pi/2:.4f}); '
    f'|X(Lam)| = {np.exp(c0 + rate*Lam):.3e}; max|X| = {np.abs(XD).max():.3f}')
for name, (nu, X, psi) in W.items():
    for lo_, hi_ in ((Lam, Lam + 10), (Lam + 10, Lam + 30), (Lam + 30, Lam + 45), (Lam + 45, 1e9)):
        s = (nu > lo_) & (nu <= hi_)
        if s.any():
            log(f'#   {name:6s} |X| on ({lo_ - Lam:+.0f}, {min(hi_, nu[-1]) - Lam:+.0f}] beyond Lam: max {np.abs(X[s]).max():.2e}, '
                f'median {np.median(np.abs(X[s])):.2e}  (grid end Lam{nu[-1]-Lam:+.1f})')

# ---------------- exact heavy images ----------------
def states(Dg):
    if HI:
        return [(1150.0, 0), (1250.0, 0)]      # J = 60 was run in this mode first (cache/tail_1000_J60.txt; numerically unvalidated, see md)
    if Dg == 1000:
        return [(1000.0, 0), (1002.0, 2), (1010.0, 10), (1030.0, 30)]
    return [(300.0, 0), (302.0, 2), (310.0, 10), (330.0, 30), (350.0, 0)]

class Img:
    def __init__(self, D, J, numax):
        self.m = np.sqrt((D - J - 2)*(D + J - 1))
        M = int(3*numax**2/(2*D) + 6*np.sqrt(numax + 1) + 60)
        self.L = S3.Laplace3(D, J, 3.0, ks=(2,), prec=400 if J < 40 else 600, N0=240, nq=320, M=M)
    def __call__(self, nu):
        v = self.L.value(2, float(nu))
        return float(v.mid())/2.0

log('')
log('# per state: Gamma = exact scalar m^4 C_{2,nu}/2 (-> 1 at nu -> 0).  A_psi = int_0^Lam psi Gamma dnu (psi-form; the published')
log('#   evidence form); A_X = Re sum_{nu<=Lam} dnu X theta Gamma; T_L / T_m = Re sum over nu > Lam / nu > m of dnu X theta Gamma (signed);')
log('#   |T_L| := sum_{nu>Lam} dnu |X| theta |Gamma|.  Norm: S_psi = int |psi| dnu (action of a state with Gamma = 1).')
for (D, J) in states(Dg):
    t0 = time.time()
    numax = max(nuD[-1], Lam + 50)
    img = Img(D, J, numax)
    m = img.m
    # coarse samples: below Lam (for A_psi), nonuniform; tail window step 0.5
    nb = np.unique(np.concatenate([np.linspace(0.02, 0.9*Lam, 181), np.linspace(0.9*Lam, Lam, 201)]))
    if J > 0 and not HI:
        nb = np.unique(np.concatenate([nb, np.linspace(0.02, Lam, max(400, 12*J))]))
    nt = np.arange(np.floor(Lam) - 2, numax + 0.5, 0.5)
    if HI:
        nt = np.unique(np.concatenate([np.arange(np.floor(Lam) - 2, m + 5, 0.5), np.arange(np.floor(m + 5), numax + 1.0, 1.0)]))
    Gb = np.array([img(v) for v in nb]); Gt = np.array([img(v) for v in nt])
    def G_on(nu):
        g = np.zeros(nu.size)
        a = nu <= nt[0]
        g[a] = np.interp(nu[a], nb, Gb)        # below the tail window: linear (fine grid)
        b = ~a
        sg = np.sign(Gt)
        ok = np.all(sg[1:] == sg[:-1])
        if ok:
            g[b] = sg[0]*np.exp(np.interp(nu[b], nt, np.log(np.abs(Gt))))
        else:
            g[b] = np.interp(nu[b], nt, Gt)
        return g
    log(f'## Delta = {D}, J = {J}: m = {m:.3f} = Lam{m - Lam:+.2f} (x1 = Lam/m = {Lam/m:.4f}); Gamma(m) = {img(m):.3e}, '
        f'Gamma(Lam+30) = {img(Lam+30):.3e}, Gamma(nu_c) = {img(nuc):.3e}, Gamma(nu_c+50) = {img(nuc+50):.3e}'
        + (f', Gamma(grid end {nuD[-1]:.1f}) = {img(nuD[-1]):.3e}' if nuD[-1] > nuc + 50 else ''))
    lg = np.log(np.abs(Gt)); sl = np.gradient(lg, nt)
    for xx in (1.0, 1.02, 1.05, 1.1, 1.2, 1.3):
        nn = xx*m
        if nt[0] <= nn <= nt[-1]:
            log(f'     d ln|Gamma|/d nu at x = {xx}: {np.interp(nn, nt, sl):.3f}')
    for name, (nu, X, psi) in W.items():
        th = theta(nu)
        g = G_on(nu)
        dn = nu[1] - nu[0]
        inL = nu <= Lam
        Apsi = np.sum(psi[inL]*g[inL])*dn
        Spsi = np.sum(np.abs(psi))*dn
        AX = np.real(np.sum(X[inL]*th[inL]*g[inL]))*dn
        TL = np.real(np.sum(X[~inL]*th[~inL]*g[~inL]))*dn
        TLa = np.sum(np.abs(X[~inL])*th[~inL]*np.abs(g[~inL]))*dn
        im = nu > m
        Tm = np.real(np.sum(X[im]*th[im]*g[im]))*dn
        Tma = np.sum(np.abs(X[im])*th[im]*np.abs(g[im]))*dn
        if im.any():
            pk = np.argmax(np.abs(X[im])*th[im]*np.abs(g[im]))
            pks = f'peak of |X theta Gamma| at nu = Lam{nu[im][pk]-Lam:+.1f}'
        else:
            pks = 'no grid point above m'
        log(f'   {name:6s}: A_psi = {Apsi:+.5e}  A_X = {AX:+.5e}  S_psi = {Spsi:.4e} | T_L = {TL:+.3e} (|T_L| {TLa:.3e}, '
            f'/|A_psi| {TLa/abs(Apsi):.1e}) | T_m = {Tm:+.3e} (|T_m| {Tma:.3e}, /|A_psi| {Tma/abs(Apsi):.1e}); {pks}')
    log(f'   ({time.time()-t0:.0f}s)')
open(os.path.join(HERE, 'cache', f'tail_{Dg}' + ('_hi' if HI else '') + '.txt'), 'w').write('\n'.join(out) + '\n')
