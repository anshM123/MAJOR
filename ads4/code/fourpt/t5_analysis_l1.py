"""
End-to-end analysis for <JJJJ> (l = 1): exact tensor Regge moments of W[v_a, v_b] (cached by t4_compute_l1.py) versus the
factorised leading prediction (predict4.py).  No fitted constant enters the prediction.

Helicity components.  At leading order every mode q = (q1,q3,q2,q4) is
   Pihat^q = Pref * sum_h T_h M^{q1q3}_{h1h3}(nu) M^{q2q4}_{h2h4}(-nu) K_{muL(h),muR(h)}(beta),   T_h = T_{-h} (parity).
The 8 parity classes of h = (h1,h2,h3,h4) are fitted (linear least squares over the 81 modes, each equation weighted by
1/|predicted mode|) with the PREDICTED weights; R_c = T_c(fit)/T_c(pred).  Classes:
   pres-pres (h3 = -h1, h4 = -h2):  'oo' = (+,+,-,-)  [opposite helicity at the vertex, |c^{+-}|^2],  'ss' = (+,-,-,+)  [same, |c^{++}|^2];
   their sum multiplies R_L R_R, their difference g1^2 S_L S_R  (the opposite-sign term):  Sigma_+- := T_oo +- T_ss;
   one pair flipped / both flipped: the remaining six classes.
1/m: R = 1 + c1/m + c2/m^2 + ...  in the Weyl variables (m^2 = alpha alphabar), fits over the masses at fixed beta.
Mutations (must produce c1 != 0): (i) b^(0) -> b^(2) (shadow coefficient of the words AE, BC); (ii) the normalisation A_+^2 of the
'same' component of V_2 replaced by Delta^2 (not Weyl-even).
"""
import os, sys, glob, itertools
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from cases import cases_fixed_beta, info              # noqa: E402
from formulas import log_b                           # noqa: E402
from tmoments4 import TensorMoments4, DATADIR
from predict4 import Mval, prefactor, Kcached, V_l1, weyl
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

HS = list(itertools.product((1, -1), repeat=4))
CLASSES = {}
for h in HS:
    key = h if h[0] == 1 else tuple(-x for x in h)
    CLASSES.setdefault(key, []).append(h)
CNAME = {(1, 1, -1, -1): 'oo', (1, -1, -1, 1): 'ss', (1, 1, 1, 1): 'ff_oo', (1, -1, 1, -1): 'ff_ss',
         (1, 1, -1, 1): 'pf_os', (1, -1, -1, -1): 'pf_so', (1, 1, 1, -1): 'fp_os', (1, -1, 1, 1): 'fp_so'}
QS = list(itertools.product((-1, 0, 1), repeat=4))

def design(Delta, J, k, nu):
    """A[q][class] = Pref * sum_{h in class} M M K  (float complex)."""
    Ap, Am, m, beta = weyl(Delta, J)
    a = 2 + k + 1
    pre = prefactor(1, Delta, J, k, nu)
    A = {}
    for q in QS:
        q1, q3, q2, q4 = q
        row = {}
        for c, hs in CLASSES.items():
            acc = mp.mpc(0)
            for h in hs:
                ML, muL = Mval(1, q1, q3, h[0], h[2], nu, a)
                if ML == 0: continue
                MR, muR = Mval(1, q2, q4, h[1], h[3], -nu, a)
                if MR == 0: continue
                acc += ML*MR*Kcached(nu, muL, muR, beta)
            row[c] = complex(pre*acc)
        A[q] = row
    return A

def hel_vectors(Delta, J):
    """exact helicity structures (t1b_helicity_l1): h_opp = v_1 + J(J+1)(m^2+2) v_2,  h_same = v_1 - (Delta-1)(Delta-2)(m^2-2) v_2"""
    D = mp.mpf(Delta)
    m2 = (D - 1.5)**2 - (J + 0.5)**2
    return (1.0, float(J*(J + 1)*(m2 + 2))), (1.0, float(-(D - 1)*(D - 2)*(m2 - 2)))

def combos(Delta, J):
    o, s_ = hel_vectors(Delta, J)
    return {'v1v1': ((1, 0), (1, 0)), 'v2v2': ((0, 1), (0, 1)), 'v1v2': ((1, 0), (0, 1)),
            'oo': (o, o), 'ss': (s_, s_), 'os': (o, s_)}

def Tpred(cL, cR, Delta, J, mut=None):
    Ap, Am, m, beta = weyl(Delta, J)
    def V(i, h1, h2):
        if mut == 'V2' and i == 2 and h2 == -h1:
            return -mp.mpf(Delta)**2/2
        return V_l1(i, Ap, Am, m, h1, h2)
    def Vc(c, h1, h2):
        return sum(c[i]*V(i+1, h1, h2) for i in range(2))
    T = {cl: Vc(cL, cl[0], cl[1])*Vc(cR, cl[2], cl[3]) for cl in CLASSES}
    if mut == 'b2':
        fac = mp.exp(log_b(mp.mpf(Delta), J, 0) - log_b(mp.mpf(Delta), J, 2))
        T = {c: v*fac for c, v in T.items()}
    return {c: float(v) for c, v in T.items()}

def measured_combo(mom, cL, cR, k, nu, logref):
    res = {}
    for q in QS:
        tot = 0j
        for a in range(2):
            for b in range(2):
                w = cL[a]*cR[b]
                if w == 0: continue
                v, sh = mom.moment(k, nu, (a, b), q)
                tot += w*complex(v)*np.exp(sh - logref)
        res[q] = (-1)**mom.J*tot
    return res

def fit_components(A, meas, Tp, logref):
    """weighted LSQ for the class amplitudes (in units of exp(logref)); returns ratios fit/pred, raw solution, relative residual."""
    cl = list(CLASSES)
    pred = {q: sum(A[q][c]*Tp[c] for c in cl) for q in QS}
    scale = max(abs(v) for v in pred.values())
    rows, rhs = [], []
    for q in QS:
        w = 1.0/max(abs(pred[q]), 1e-9*scale)
        rows.append([A[q][c]*np.exp(-logref)*w for c in cl]); rhs.append(meas[q]*w)
    M = np.array(rows); y = np.array(rhs)
    Mr = np.vstack([M.real, M.imag]); yr = np.concatenate([y.real, y.imag])
    sol, *_ = np.linalg.lstsq(Mr, yr, rcond=None)
    resid = np.linalg.norm(Mr @ sol - yr)/np.linalg.norm(yr)
    Tmax = max(abs(v) for v in Tp.values())
    return {c: (sol[i]/Tp[c] if abs(Tp[c]) > 1e-6*Tmax else None) for i, c in enumerate(cl)}, dict(zip(cl, sol)), resid

def run_case(mom, ks=(4, 6), nus=(0.0, 1.0, 3.0), names=None, mut=None):
    D, J = mom.Delta, mom.J
    cmb = combos(D, J)
    names = names or list(cmb)
    out = {}
    for k in ks:
        for nu in nus:
            A = design(D, J, k, nu)
            for name in names:
                cL, cR = cmb[name]
                Tp = Tpred(cL, cR, D, J, mut)
                big = max(abs(sum(A[q][c]*Tp[c] for c in CLASSES)) for q in QS)
                logref = np.log(big)
                meas = measured_combo(mom, cL, cR, k, nu, logref)
                R, sol, resid = fit_components(A, meas, Tp, logref)
                oo, ss = (1, 1, -1, -1), (1, -1, -1, 1)
                sp_, sm_ = sol[oo] + sol[ss], sol[oo] - sol[ss]
                pp, pm = Tp[oo] + Tp[ss], Tp[oo] - Tp[ss]
                q0 = (1, -1, 1, -1)
                pq = sum(A[q0][c]*Tp[c] for c in CLASSES)*np.exp(-logref)
                out[(k, nu, name)] = dict(R=R, resid=resid, Rp=(sp_/pp if abs(pp) > 1e-6*max(abs(pp), abs(pm)) else None),
                                          Rm=(sm_/pm if abs(pm) > 1e-6*max(abs(pp), abs(pm)) else None), Rmode=(meas[q0]/pq).real,
                                          sol=sol, Tp=Tp)
    return out

def fits(ms, Rs):
    ms = np.array(ms, dtype=float); Rs = np.array(Rs, dtype=float)
    o = np.argsort(ms); ms, Rs = ms[o], Rs[o]
    res = {}
    if len(ms) >= 3:
        m3, R3 = ms[-3:], Rs[-3:]
        A = np.vstack([1/m3, 1/m3**2, 1/m3**4]).T; res['A'] = np.linalg.solve(A, R3 - 1)
        A = np.vstack([1/m3**2, 1/m3**4, 1/m3**6]).T; res['B'] = np.linalg.solve(A, R3 - 1)
    if len(ms) >= 4:
        m4, R4 = ms[-4:], Rs[-4:]
        A = np.vstack([1/m4, 1/m4**2, 1/m4**3, 1/m4**4]).T; res['C'] = np.linalg.solve(A, R4 - 1)
    return res

def fmt(x):
    return "   --  " if x is None else f"{x:.4f}"

LEAD = {'v1v1': [(1, 1, -1, -1), (1, -1, -1, 1), (1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, 1)],
        'v2v2': [(1, 1, -1, -1), (1, -1, -1, 1), (1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, 1)],
        'v1v2': [(1, 1, -1, -1), (1, -1, -1, 1), (1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, 1)],
        'oo': [(1, 1, -1, -1), (1, 1, 1, 1)], 'ss': [(1, -1, -1, 1), (1, -1, 1, -1)], 'os': [(1, 1, -1, 1), (1, 1, 1, -1)]}

if __name__ == "__main__":
    pr("# <JJJJ> end-to-end: exact tensor Regge moments / factorised leading prediction (fourpt/t5_analysis_l1.py)")
    pr("# combos: W[c_L, c_R] with c on (v_1, v_2); 'o' = exact opposite-helicity structure, 's' = exact same-helicity structure (t1b)")
    pr("# classes h = (h1,h2,h3,h4): oo=(+,+,-,-) [pres-pres, |c^{+-}|^2], ss=(+,-,-,+) [pres-pres, |c^{++}|^2], ff_oo=(+,+,+,+), ff_ss=(+,-,+,-) [both flipped],")
    pr("#   pf=(+,+,-,+) [pair 13 preserved, 24 flipped], fp=(+,+,+,-);  Sigma_+- = T_oo +- T_ss: Sigma_+ multiplies R_L R_R, Sigma_- multiplies g1^2 S_L S_R")
    files = sorted(glob.glob(os.path.join(DATADIR, 'modes4_l1_*.npz')))
    moms = [TensorMoments4.load(f) for f in files]
    data = {}
    for mom in moms:
        m, beta, _ = info(mom.Delta, mom.J)
        bkey = round(float(beta), 1)
        res = run_case(mom)
        resb = run_case(mom, ks=(6,), nus=(0.0, 3.0), names=['v1v1', 'v2v2', 'v1v2'], mut='b2')
        resv = run_case(mom, ks=(6,), nus=(0.0, 3.0), names=['v1v1', 'v2v2', 'v1v2'], mut='V2')
        data.setdefault(bkey, []).append((float(m), mom.Delta, mom.J, res, resb, resv, mom))
    for bkey in sorted(data):
        rows = sorted(data[bkey], key=lambda r: r[0])
        pr(f"\n## beta ~ {bkey}")
        for k in (6, 4):
            for nu in (0.0, 1.0, 3.0):
                pr(f"# k={k} nu={nu}: per combo: R of the leading classes | R(Sigma_+) R(Sigma_-) | LSQ residual | ratio of the raw mode (1,-1,1,-1)")
                for (m, D, J, res, resb, resv, mom) in rows:
                    pr(f"  m={m:6.2f} (Delta,J)=({D:.3f},{J})")
                    for name in ('v1v1', 'v2v2', 'v1v2', 'oo', 'ss', 'os'):
                        r = res[(k, nu, name)]
                        pr(f"      {name:5s}: " + " ".join(f"{CNAME[c]}={fmt(r['R'][c])}" for c in LEAD[name])
                           + f" | S+={fmt(r['Rp'])} S-={fmt(r['Rm'])} | {r['resid']:.1e} | {r['Rmode']:.4f}")
                mom = rows[-1][6]
                pr(f"     numerical tail estimate (largest m, W[v1,v1], mode (1,-1,1,-1)): {mom.tail_estimate(k, nu, (0, 0), (1, -1, 1, -1)):.1e}")
        pr(f"\n# Lemma S (two pairs) at beta ~ {bkey}: exact mode symmetries")
        for (m, D, J, res, resb, resv, mom) in rows:
            dev1 = dev2 = 0.0
            v0, s0 = mom.moment(6, 1.0, (0, 1), (1, -1, 1, -1)); ref = abs(v0)
            for q in QS:
                v, s_ = mom.moment(6, 1.0, (0, 1), q)
                v1, s1 = mom.moment(6, 1.0, (1, 0), (q[1], q[0], q[3], q[2]))
                v2, s2 = mom.moment(6, 1.0, (0, 1), tuple(-x for x in q))
                dev1 = max(dev1, abs(v*np.exp(s_ - s0) - v1*np.exp(s1 - s0))/ref); dev2 = max(dev2, abs(v*np.exp(s_ - s0) - v2*np.exp(s2 - s0))/ref)
            pr(f"  m={m:.1f}: max_q |P^q[v1,v2] - P^q'[v2,v1]|/|P^(1,-1,1,-1)| = {dev1:.1e} (q' = (q3,q1,q4,q2));  max_q |P^q - P^-q| = {dev2:.1e}")
        pr(f"\n# 1/m fits at beta ~ {bkey}, k = 6: A = (c1/m + c2/m^2 + c4/m^4; 3 largest m) -> c1;  B = (c2, c4, c6) -> c2;  C = (c1..c4; 4 largest m) -> c1")
        for nu in (0.0, 1.0, 3.0):
            for name in ('v1v1', 'v2v2', 'v1v2', 'oo', 'ss', 'os'):
                for c in LEAD[name][:2] + (['Sp', 'Sm'] if name in ('v2v2', 'v1v2') else []):
                    ms = [r[0] for r in rows]
                    if c == 'Sp': Rs = [r[3][(6, nu, name)]['Rp'] for r in rows]; lab = 'Sigma+'
                    elif c == 'Sm': Rs = [r[3][(6, nu, name)]['Rm'] for r in rows]; lab = 'Sigma-'
                    else: Rs = [r[3][(6, nu, name)]['R'][c] for r in rows]; lab = CNAME[c]
                    if any(x is None for x in Rs): continue
                    f = fits(ms, Rs)
                    line = f"  nu={nu} {name:5s} {lab:7s}: R = " + " ".join(f"{x:.5f}" for x in Rs)
                    if 'A' in f: line += f" | c1(A) = {f['A'][0]:+.4f}  c2(B) = {f['B'][0]:+.1f}"
                    if 'C' in f: line += f"  c1(C) = {f['C'][0]:+.4f}"
                    pr(line)
        pr(f"# mutation tests at beta ~ {bkey}, k = 6: c1(A) of the pres-pres classes with b^(0) -> b^(2), and with A_+^2 -> Delta^2 in V_2")
        for nu in (0.0, 3.0):
            for name in ('v1v1', 'v2v2', 'v1v2'):
                for c in ((1, 1, -1, -1), (1, -1, -1, 1)):
                    ms = [r[0] for r in rows]
                    fb = fits(ms, [r[4][(6, nu, name)]['R'][c] for r in rows]); fv = fits(ms, [r[5][(6, nu, name)]['R'][c] for r in rows])
                    if 'A' in fb:
                        pr(f"  nu={nu} {name} {CNAME[c]}: c1(A) with b^(2): {fb['A'][0]:+.4f};  with Delta^2 in V_2: {fv['A'][0]:+.4f}")
    open(os.path.join(HERE, 'out_t5_analysis_l1.txt'), 'w').write("\n".join(out) + "\n")
