"""Referee check of t13_analysis_l2: reproduce fits from the saved mode caches, design conditioning, mutation tests.
Usage: python3 ref_t13.py <npz>   (read-only use of fourpt/tttt code and caches)"""
import os, sys, itertools, pickle
import numpy as np, mpmath as mp
T = '/home/user/MAJOR/ads4/code/fourpt/tttt'
sys.path.insert(0, T)
import t13_analysis_l2 as A13
from t13_analysis_l2 import Modes, CLASSES, CNAME, QS, Tpred, lam, L
from predict4 import Mfun, prefactor, Kcached, weyl
from formulas import gamma_a
from tmoments_l2 import COMBOS

fn = sys.argv[1]; k, nu, cT = 6, 0.0, 8.0
mo = Modes(fn); D, J = mo.Delta, mo.J
Ap, Am, m, beta = weyl(D, J)
print(f"# {os.path.basename(fn)} D={D} J={J} m={float(m):.4f} beta={float(beta):.4f}")
print("combos in cache:", COMBOS, " grid points:", mo.g['U'].shape, " conv&ok:", int((mo.g['ok'] & mo.conv).sum()))

# --- design split in R/S pieces so that the g1^2 S_L S_R ('opposite-sign') term can be mutated
a = 2*L + k + 1
pre = prefactor(L, D, J, k, nu)
cache = f'/home/user/MAJOR/ads4/code/referee_tttt/cache/design_{J}.pkl'
if os.path.exists(cache):
    parts = pickle.load(open(cache, 'rb'))
else:
    MF = Mfun(L)
    def RS(q1, q3, h1, h3, nu_):
        R, S, mu = MF[(q1, q3, h1, h3)]
        g1 = gamma_a(a - 1, nu_)/gamma_a(a, nu_)
        return R(nu_, a), g1*S(nu_, a), mu
    parts = {}
    for q in QS:
        q1, q3, q2, q4 = q
        for c, hs in CLASSES.items():
            acc = np.zeros(4, complex)
            for h in hs:
                RL, SL, muL = RS(q1, q3, h[0], h[2], nu)
                RR, SR, muR = RS(q2, q4, h[1], h[3], -nu)
                K = Kcached(nu, muL, muR, beta)
                acc += np.array([complex(pre*RL*RR*K), complex(pre*RL*SR*K), complex(pre*SL*RR*K), complex(pre*SL*SR*K)])
            parts[(q, c)] = acc
    pickle.dump(parts, open(cache, 'wb'))

def design(sgn=(1, 1, 1, 1)):
    s = np.array(sgn)
    return {q: {c: complex(np.sum(s*parts[(q, c)])) for c in CLASSES} for q in QS}

A0 = design()
# sanity: compare with t13 design for a few q
Aref = A13.design(D, J, k, nu)
err = max(abs(A0[q][c] - Aref[q][c])/max(1e-300, max(abs(Aref[q][cc]) for cc in CLASSES)) for q in QS[:40] for c in CLASSES)
print("design split vs t13.design (40 q):", f"{err:.1e}")

# measurements
lo, ls = float(lam('opp', mp.mpf(D), J)), float(lam('same', mp.mpf(D), J))
w_mix = np.sqrt(2*ls/lo)
waves = {c: c for c in COMBOS}
waves['mix'] = {('same', 'same'): 1.0, ('same', 'opp'): w_mix, ('opp', 'same'): w_mix, ('opp', 'opp'): w_mix**2}
logref = None; meas = {}
for name, cmb in waves.items():
    mm = {}
    for q in QS:
        v, sh = mo.moment(k, nu, cmb, q)
        if logref is None: logref = sh
        mm[q] = (-1)**J*complex(v)*np.exp(sh - logref)
    meas[name] = mm

def Tp_of(name, lamswap=False):
    def tp(hL, hR):
        if lamswap:
            sw = {'same': 'opp', 'opp': 'same'}
            t = Tpred(sw[hL], sw[hR], D, J, cT)  # value with swapped lambdas ...
            cls = Tpred(hL, hR, D, J, cT)        # ... on the correct helicity classes
            v = max(t.values())
            return {c: (v if cls[c] != 0 else 0.0) for c in CLASSES}
        return Tpred(hL, hR, D, J, cT)
    if name == 'mix':
        return {c: sum(wt*tp(hL, hR)[c] for (hL, hR), wt in waves['mix'].items()) for c in CLASSES}
    return tp(*name)

def run(A, tag, lamswap=False):
    print(f"--- {tag}")
    cl = list(CLASSES)
    M = np.array([[A[q][c]*np.exp(-logref) for c in cl] for q in QS])
    Mr = np.vstack([M.real, M.imag])
    if tag == 'baseline':
        cn = np.linalg.norm(Mr, axis=0)
        print("   design column norms (rel):", {CNAME[c]: f"{v/cn.max():.1e}" for c, v in zip(cl, cn)}, " cond(Mr) =", f"{np.linalg.cond(Mr):.1e}",
              " cond(cols-normalised) =", f"{np.linalg.cond(Mr/cn):.1e}", " rank-params:", Mr.shape)
    for name in waves:
        Tp = Tp_of(name, lamswap)
        sol, resid, s1, r1 = A13.fit(A, meas[name], Tp, logref)
        Tm = max(abs(v) for v in Tp.values())
        R = {CNAME[c]: round((sol[c]/Tp[c]).real, 4) for c in CLASSES if abs(Tp[c]) > 1e-6*Tm}
        ex = ''
        if name == 'mix':
            oo, ss = (1, 1, -1, -1), (1, -1, -1, 1)
            ex = f" Sig+ {((sol[oo]+sol[ss])/(Tp[oo]+Tp[ss])).real:.4f} Sig- {((sol[oo]-sol[ss])/(Tp[oo]-Tp[ss])).real:.4f}"
        print(f"   {str(name):18s} one-scale {s1.real:+.4f} (resid {r1:.3f}); class LSQ resid {resid:.3f}; classes {R}{ex}")

run(A0, 'baseline')
run(A0, 'mutation 1: swap lambda_same <-> lambda_opp in T_pred', lamswap=True)
run(design((1, 1, 1, -1)), 'mutation 2: flip sign of g1^2 S_L S_R (opposite-sign) term in design')
run(design((1, -1, -1, 1)), 'mutation 3: g1 -> -g1 (flip R S cross terms)')
# mutation 4: swap the helicity class assignment same <-> opp at the left vertex in the prediction
print("--- mutation 4: predicted helicity classes swapped (same <-> opp) on both vertices, values kept")
sw = {'same': 'opp', 'opp': 'same'}
for name in [c for c in COMBOS]:
    t_true = Tpred(*name, D, J, cT); t_cls = Tpred(sw[name[0]], sw[name[1]], D, J, cT)
    v = max(t_true.values()); Tp = {c: (v if t_cls[c] != 0 else 0.0) for c in CLASSES}
    sol, resid, s1, r1 = A13.fit(A0, meas[name], Tp, logref)
    print(f"   {str(name):18s} one-scale {s1.real:+.4f} (resid {r1:.3f})")

# --- diagnostics: per-class contribution norms, norm of measured vs predicted, pair-exchange identity so <-> os
print("--- diagnostics")
cl = list(CLASSES)
for name in waves:
    Tp = Tp_of(name)
    sol, resid, s1, r1 = A13.fit(A0, meas[name], Tp, logref)
    y = np.array([meas[name][q] for q in QS])
    pv = np.array([sum(A0[q][c]*Tp[c] for c in cl) for q in QS])*np.exp(-logref)
    contrib = {CNAME[c]: np.linalg.norm(np.array([A0[q][c] for q in QS])*np.exp(-logref)*sol[c])/np.linalg.norm(y) for c in cl}
    print(f"   {str(name):18s} |meas|/|pred| = {np.linalg.norm(y)/np.linalg.norm(pv):.4f}; fitted class contributions/|meas|:",
          {k: f"{v:.1e}" for k, v in contrib.items()})
so, os_ = meas[('same', 'opp')], meas[('opp', 'same')]
d1 = max(abs(so[q] - os_[(q[2], q[3], q[0], q[1])]) for q in QS)/max(abs(v) for v in so.values())
d2 = max(abs(so[q] - os_[q]) for q in QS)/max(abs(v) for v in so.values())
print(f"   so vs os under (q1,q3)<->(q2,q4): {d1:.1e};  so vs os same q: {d2:.1e}")
# deviation coefficient m^2 (ratio - 1) for the one-scale fits
for name in waves:
    Tp = Tp_of(name); s1 = A13.fit(A0, meas[name], Tp, logref)[2].real
    print(f"   {str(name):18s} m^2 (one-scale - 1) = {float(m)**2*(s1-1):+.1f}")
