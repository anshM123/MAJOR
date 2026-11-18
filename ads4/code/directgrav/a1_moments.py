"""
a1: exact <TTTT> kernel moments of the directgrav waves (or any tttt-format mode file) at the low Gauss nodes of the functionals,
against the factorised leading prediction (t13/t16 machinery, imported read-only), plus the exact scalar skeleton at the same (Delta, J, nu).

Per (k, nu, wave, kernel):  R_grav = measured kernel moment / leading prediction (raw, no corner estimate; as t16 'kernel_moment').
Waves: ss = W[X_s,X_s] (++ state), oo = W[X_o,X_o] (+- state), so = mixed (interference; reported as |K(so)|/sqrt|K(ss)K(oo)|).
Kernels: FF = S_ff x S_ff (even k, D-safe / Omega^{NF,SC} and Omega^SC), E = D1^2 x D1^2 (even k, Omega^U), O = A x A (odd k, both).
Skeleton (Df = 3, 2): R_ex = m^{2k} C_k(exact, exactskel) / (2 P_{-1/2+i nu}(cosh beta));  R_pos = position-space moment / prediction on the
same grid (calibration of the position-space route, c1).
Usage: python3 a1_moments.py <modes.npz> Dg mfl [numax]   -> cache/a1_<tag>.pkl (+ printed table)
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
for p in (TT, os.path.join(TT, '..'), os.path.join(TT, '..', 'oddk'), os.path.join(TT, '..', '..', 'e2e'),
          os.path.join(HERE, '..', 'exactskel'), os.path.join(HERE, '..', 'omegau')):
    sys.path.insert(0, p)
import numpy as np, mpmath as mp
from t13_analysis_l2 import Modes, design, Tpred
from t16_kernels_ext import rel

# Circle-mode weights of the TENSOR kernels (corrected normalisation n_|q| = (1/36, 1/16, 1), fourpt/tttt/t18_kstar.py and
# SPINNING_TTTT_EXT.md sect. 6.2; the naive weights of t16_kernels_ext.KER over-weight the isotropic mode and are NOT used here).
KER = {
    'E': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25},                                            # D1^2 (Omega^U even kernel)
    'FF': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/32, (-1, 1): -1/32},           # S_ff = D1^2 + D1 D2 (flip-free)
    'KS': {(2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/64, (-1, 1): -1/64},                         # K* = (D1 + D2)(D1 - D2/2)
    'O': {(2, -2): 1/(4j), (-2, 2): -1/(4j)},                                                     # A = D5 D1 (odd kernel)
}


def kernel_moment(mo, k, nu, wave, kern, A, Tp):
    """as t16_kernels_ext.kernel_moment (raw measured and predicted two-pair kernel moments), with the corrected KER above"""
    from t13_analysis_l2 import CLASSES
    K = KER[kern]
    terms = []; pred = 0j
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4)
            v, sh = mo.moment(k, nu, wave, q)
            terms.append((w13*w24, v, sh))
            pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    shm = max(t[2] for t in terms)
    meas = sum(w*v*np.exp(sh - shm) for w, v, sh in terms)*(-1)**mo.J
    return complex(meas), complex(pred), shm
from moments import ScalarMoments
from t3_scalar_mixed import prediction as sprediction
from skelexact import ExactState
from common import gl01
import remk
mp.mp.dps = 20

fn = sys.argv[1]; Dg = int(sys.argv[2]); mfl = float(sys.argv[3])
NUMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 16.0
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
pn, pw = gl01(64)
mo = Modes(fn); D, J = mo.Delta, mo.J
al, alb = D - J - 2, D + J - 1
m = np.sqrt(al*alb); eta = (al*al + alb*alb)/(2*al*alb)
inodes = [i for i in range(64) if Lam*pn[i] <= NUMAX]
NUS = [float(Lam*pn[i]) for i in inodes] + [0.0, 1.0, 3.0]
KS = (2, 3, 4, 5, 6, 7)
WAVES = [('same', 'same'), ('opp', 'opp'), ('same', 'opp')]
tag = os.path.basename(fn).replace('.npz', '')
out = os.path.join(HERE, 'cache', 'a1_%s.pkl' % tag)
db = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
db.update(dict(Delta=D, J=J, m=m, eta=eta, beta=float(np.arccosh(eta)), Lam=Lam, mfl=mfl, Dg=Dg, inodes=inodes, NUS=NUS))
print(f"# a1: {tag}: Delta={D:.6f} J={J} m={m:.4f} beta={np.arccosh(eta):.5f}; {len(NUS)} nu values (nu <= {NUMAX})", flush=True)
t0 = time.time()
G = db.setdefault('grav', {})
for nu in NUS:
    for k in KS:
        if (k, nu) in G:
            continue
        A = design(D, J, k, nu)
        row = {}
        for wave in WAVES:
            if wave[0] == 'same' and J % 2 == 1:
                continue
            if 'opp' in wave and J < 4:
                continue
            Tp = Tpred(wave[0], wave[1], D, J, 8.0)
            for kern in KER:
                row[(wave, kern)] = kernel_moment(mo, k, nu, wave, kern, A, Tp)
        G[(k, nu)] = row
    print(f"  grav nu={nu:.4f} done ({time.time() - t0:.0f} s)", flush=True)
    pickle.dump(db, open(out + '.tmp', 'wb')); os.replace(out + '.tmp', out)
# skeleton
nus_ex = np.array([max(nu, 1e-4) for nu in NUS])
Pd = remk.Pders_vec(nus_ex, eta, nmax=0)[0]
db['P'] = Pd
for Df in (3.0, 2.0):
    if ('skel', Df) in db:
        continue
    es = ExactState(D, J, Df); ex = es.values(nus_ex)
    sm = ScalarMoments(D, J, Df, Df, 8000)
    rec = {}
    for k in KS:
        for i, nu in enumerate(NUS):
            sg, lg = sm.moment(k, nu)
            Rp = float(sg*mp.exp(mp.mpf(lg))/sprediction(D, J, Df, Df, k, nu).real)
            rec[(k, nu)] = dict(ex=float(ex[k][i]), R_ex=float(ex[k][i]/(2*Pd[i])), R_pos=Rp)
    db[('skel', Df)] = rec
    print(f"  skeleton Df={Df} done ({time.time() - t0:.0f} s)", flush=True)
    pickle.dump(db, open(out + '.tmp', 'wb')); os.replace(out + '.tmp', out)

def Rg(k, nu, wave, kern):
    meas, pred, shm = G[(k, nu)][(wave, kern)]
    return rel(meas, pred, shm)

print("#  k   nu      x     | ss: FF       E        K*       O      | oo: FF       E        K*       O      | so: |FF|/sqrt  | skel Df3: R_ex   R_pos-R_ex")
for k in KS:
    for nu in sorted(NUS):
        cells = []
        for wave in WAVES[:2]:
            if (wave, 'FF') in G[(k, nu)]:
                cells.append(" ".join(f"{Rg(k, nu, wave, kn).real:+.4f}" for kn in ('FF', 'E', 'KS', 'O')))
            else:
                cells.append(" " * 31 + "-")
        so = ''
        if (WAVES[2], 'FF') in G[(k, nu)]:
            ms, ps, s1 = G[(k, nu)][(WAVES[0], 'FF')]; mo_, po, s2 = G[(k, nu)][(WAVES[1], 'FF')]; mx, px, s3 = G[(k, nu)][(WAVES[2], 'FF')]
            so = f"{np.exp(np.log(abs(mx) + 1e-300) + s3 - 0.5*(np.log(abs(ms)) + s1 + np.log(abs(mo_)) + s2)):.1e}"
        sk = db[('skel', 3.0)][(k, nu)]
        print(f"  {k}  {nu:7.3f}  {nu/m:.4f} | " + " | ".join(cells) + f" | {so:>8s} | {sk['R_ex']:+.5f} {sk['R_pos'] - sk['R_ex']:+.1e}", flush=True)
print(f"# total {time.time() - t0:.0f} s")
