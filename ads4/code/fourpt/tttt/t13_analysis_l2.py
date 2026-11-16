"""
End-to-end <TTTT> check (coordinator step 2 of the third request): exact tensor Regge moments of the Taylor-mode waves W[X_hL, X_hR]
(tmoments_l2 caches) versus the factorised leading prediction with the l = 2 heavy vertices (t14):
   (-1)^J Pihat^q_{k,nu}[s13^3 s24^3 W[X_hL, X_hR]] -> Pref * sum_h T^{hL hR}_h M^{q1q3}_{h1h3}(nu) M^{q2q4}_{h2h4}(-nu) K_{muL,muR}(beta),
   Pref = (pi/4) 4^6/b^(0) gamma_a^2 m^{-2a}, a = 2l+k+1 = 11 (k = 6),
   T^{hL hR}_h = c_T sqrt(lambda_hL lambda_hR) m^-8 on the helicity classes of (hL, hR) [same: h2 = -h1; opp: h2 = h1], 0 otherwise
   (leading heavy vertices of t14: T^{hh} = lambda_h E, E = c_T/m^8 Weyl-even, off-diagonal = geometric mean; c_T = 8 <=> 1/4 per quadruple).
Comparisons (no fitted constant in the prediction):
  * class fits (8 parity classes, 625 modes, predicted weights) -> R_c = T_c(fit)/T_c(pred) for each wave;
  * one-scale fit of all 625 modes per wave;
  * Sigma_+- for the mixed wave W[e, e], e = X_same + w X_opp (w chosen so that both helicities contribute equally):
    Sigma_+ = T_oo + T_ss (multiplies R_L R_R), Sigma_- = T_oo - T_ss (multiplies the opposite-sign term g1^2 S_L S_R);
  * prediction-free factorisation  T^{so} T^{os}/(T^{ss} T^{oo}) = 1 from the fitted class amplitudes.
Usage: python3 t13_analysis_l2.py cache/modes_l2_....npz c_T [k nu]
"""
import os, sys, itertools
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from predict4 import Mval, prefactor, Kcached, weyl
from moments4 import kernel_np
from tmoments_l2 import COMBOS
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
L = 2
HS = list(itertools.product((1, -1), repeat=4))
CLASSES = {}
for h in HS:
    key = h if h[0] == 1 else tuple(-x for x in h)
    CLASSES.setdefault(key, []).append(h)
CNAME = {(1, 1, -1, -1): 'oo', (1, -1, -1, 1): 'ss', (1, 1, 1, 1): 'ff_oo', (1, -1, 1, -1): 'ff_ss',
         (1, 1, -1, 1): 'pf_os', (1, -1, -1, -1): 'pf_so', (1, 1, 1, -1): 'fp_os', (1, -1, 1, 1): 'fp_so'}
QS = list(itertools.product(range(-L, L + 1), repeat=4))

class Modes:
    def __init__(self, fn):
        d = np.load(fn)
        self.Delta, self.J = float(d['Delta']), int(d['J'])
        self.g = dict(us=d['us'], ts=d['ts'], U=d['U'], S=d['S'], W=d['W'], ok=d['ok'])
        self.conv, self.logscale = d['conv'], d['logscale']
        self.modes = {}
        for key, v in zip(d['keys'], d['modes']):
            self.modes.setdefault(COMBOS[int(key[0])], {})[tuple(int(x) for x in key[1:])] = v
        self._K = {}

    def moment(self, k, nu, combo, q):
        g = self.g; U, S = g['U'], g['S']
        xi = -np.log(g['us'])
        q1, q3, q2, q4 = q
        key = (nu, q1 + q3, q2 + q4)
        if key not in self._K:
            self._K[key] = kernel_np(nu, [(q1 + q3, q2 + q4)], xi)[(q1 + q3, q2 + q4)]
        F = self._K[key]
        rho, rhob = S**2, (S*U)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + self.logscale
        logw[~(g['ok'] & self.conv)] = -np.inf
        sh = np.max(logw)
        G = self.modes[combo][q] if not isinstance(combo, dict) else sum(w*self.modes[c][q] for c, w in combo.items())
        val = np.sum(np.where(np.isfinite(logw), np.exp(logw - sh), 0.0)*F[:, None]*G*g['W'])
        return val, sh

def design(Delta, J, k, nu):
    Ap, Am, m, beta = weyl(Delta, J)
    a = 2*L + k + 1
    pre = prefactor(L, Delta, J, k, nu)
    A = {}
    for q in QS:
        q1, q3, q2, q4 = q
        row = {}
        for c, hs in CLASSES.items():
            acc = mp.mpc(0)
            for h in hs:
                ML, muL = Mval(L, q1, q3, h[0], h[2], nu, a)
                if ML == 0: continue
                MR, muR = Mval(L, q2, q4, h[1], h[3], -nu, a)
                if MR == 0: continue
                acc += ML*MR*Kcached(nu, muL, muR, beta)
            row[c] = complex(pre*acc)
        A[q] = row
    return A

def lam(h, D, J):
    P = lambda D_, J_: mp.rf((2 - D_ + J_)/2, 4)*mp.rf((D_ + J_ - 7)/2, 4)
    A_ = mp.rf(-J, 4)/mp.rf(J + 1, 4) if h == 'opp' else mp.rf(D - 1, 4)/mp.rf(2 - D, 4)
    return 2**(2*J + 1)*P(3 - D, -1 - J)/(A_*P(D, J))

def Tpred(hL, hR, Delta, J, cT):
    Ap, Am, m, beta = weyl(Delta, J)
    D = mp.mpf(Delta)
    val = mp.sqrt(lam(hL, D, J)*lam(hR, D, J))*cT*m**-8
    T = {}
    for cl in CLASSES:
        okL = (cl[1] == -cl[0]) if hL == 'same' else (cl[1] == cl[0])
        okR = (cl[3] == -cl[2]) if hR == 'same' else (cl[3] == cl[2])
        T[cl] = float(val) if (okL and okR) else 0.0
    return T

def fit(A, meas, Tp, logref):
    cl = list(CLASSES)
    pred = {q: sum(A[q][c]*Tp[c] for c in cl) for q in QS}
    scale = max(abs(v) for v in pred.values())
    rows, rhs = [], []
    for q in QS:
        w = 1.0/max(abs(pred[q]), 1e-3*scale)
        rows.append([A[q][c]*np.exp(-logref)*w for c in cl]); rhs.append(meas[q]*w)
    M = np.array(rows); y = np.array(rhs)
    Mr = np.vstack([M.real, M.imag]); yr = np.concatenate([y.real, y.imag])
    sol, *_ = np.linalg.lstsq(Mr, yr, rcond=None)
    resid = np.linalg.norm(Mr @ sol - yr)/np.linalg.norm(yr)
    # one-scale fit with uniform weights
    pv = np.array([pred[q]*np.exp(-logref) for q in QS]); mv = np.array([meas[q] for q in QS])
    s1 = np.vdot(pv, mv)/np.vdot(pv, pv)
    r1 = np.linalg.norm(mv - s1*pv)/np.linalg.norm(mv)
    return dict(zip(cl, sol)), resid, s1, r1

if __name__ == "__main__":
    fn = sys.argv[1]; cT = float(sys.argv[2])
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    nu = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    mo = Modes(fn)
    D, J = mo.Delta, mo.J
    Ap, Am, m, beta = weyl(D, J)
    pr(f"# <TTTT> end to end (fourpt/tttt/t13_analysis_l2.py): {os.path.basename(fn)}: Delta={D:.6f} J={J} m={mp.nstr(m, 6)} beta={mp.nstr(beta, 4)}; k={k} nu={nu}; c_T={cT}")
    A = design(D, J, k, nu)
    logref = None
    fits = {}
    waves = {c: c for c in COMBOS}
    lo, ls = float(lam('opp', mp.mpf(D), J)), float(lam('same', mp.mpf(D), J))
    w_mix = np.sqrt(ls/lo)                      # equal predicted weight of the two helicities in W[e,e] (T^{hh} ~ lambda_h)
    mix = {('same', 'same'): 1.0, ('same', 'opp'): w_mix, ('opp', 'same'): w_mix, ('opp', 'opp'): w_mix**2}
    waves['mix'] = mix
    for name, cmb in waves.items():
        meas = {}
        for q in QS:
            v, sh = mo.moment(k, nu, cmb, q)
            if logref is None: logref = sh
            meas[q] = (-1)**J*complex(v)*np.exp(sh - logref)
        if name == 'mix':
            Tp = {c: sum(wt*Tpred(hL, hR, D, J, cT)[c] for (hL, hR), wt in mix.items()) for c in CLASSES}
        else:
            Tp = Tpred(name[0], name[1], D, J, cT)
        sol, resid, s1, r1 = fit(A, meas, Tp, logref)
        Tmax = max(abs(v) for v in Tp.values())
        R = {CNAME[c]: sol[c]/Tp[c] for c in CLASSES if abs(Tp[c]) > 1e-6*Tmax}
        fits[name] = (sol, Tp)
        extra = ""
        if name == 'mix':
            Sp = (sol[(1, 1, -1, -1)] + sol[(1, -1, -1, 1)])/(Tp[(1, 1, -1, -1)] + Tp[(1, -1, -1, 1)])
            dm = Tp[(1, 1, -1, -1)] - Tp[(1, -1, -1, 1)]
            Sm = (sol[(1, 1, -1, -1)] - sol[(1, -1, -1, 1)])
            extra = f";  Sigma_+ fit/pred = {Sp:.4f};  Sigma_- = T_oo - T_ss: fit {Sm:.4e} vs pred {dm:.4e} (relative to Sigma_+ pred: {Sm/(Tp[(1,1,-1,-1)] + Tp[(1,-1,-1,1)]):.2e})"
        offs = max(abs(sol[c]) for c in CLASSES if abs(Tp[c]) <= 1e-6*Tmax) if any(abs(Tp[c]) <= 1e-6*Tmax for c in CLASSES) else 0
        pr(f"  wave {name}: class ratios fit/pred " + ", ".join(f"{n}: {v:.4f}" for n, v in R.items()) +
           f";  |classes predicted 0|/T_pred max: {offs/Tmax:.1e};  LSQ resid {resid:.1e};  one-scale fit {s1.real:.4f}{s1.imag:+.1e}i (resid {r1:.1e})" + extra)
    # prediction-free factorisation: T^{so} T^{os} / (T^{ss} T^{oo}) from fitted classes (average over the two classes of each wave)
    def amp(name, classes):
        return np.mean([fits[name][0][c] for c in classes])
    Tss = amp(('same', 'same'), [(1, -1, -1, 1), (1, -1, 1, -1)]); Too = amp(('opp', 'opp'), [(1, 1, -1, -1), (1, 1, 1, 1)])
    Tso = amp(('same', 'opp'), [(1, -1, 1, 1), (1, -1, -1, -1)]); Tos = amp(('opp', 'same'), [(1, 1, 1, -1), (1, 1, -1, 1)])
    pr(f"  prediction-free factorisation T^so T^os/(T^ss T^oo) = {Tso*Tos/(Tss*Too):.4f};  T^so/T^os = {Tso/Tos:.4f} (leading prediction 1);  T^oo/T^ss = {Too/Tss:.4f} vs lambda_o/lambda_s = {lo/ls:.4f}")
    open(os.path.join(HERE, 'out_t13_analysis_l2_%s_k%d_nu%g.txt' % (os.path.basename(fn)[:-4], k, nu)), 'w').write("\n".join(out) + "\n")
