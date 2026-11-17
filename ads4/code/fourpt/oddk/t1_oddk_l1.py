"""
Odd-k (k = 3, 5, 7; with k = 4, 6 for comparison) tensor Regge moments of the exact <JJJJ> partial waves with BOTH pairs spinning,
end to end against the factorised leading formula of SPINNING_4PT.md sect. 4.2 continued in k (predict4: a = k + 3, M = R + g1 S with
the one-pair tables at that a, gamma_a(nu), kernel K; nothing fitted).  Data: the cached mode grids of fourpt (t4_compute_l1), completed
where available by complete_n2 (N2 = 30000 at the dropped tiny-t points); the remaining dropped points are estimated by corner_est
(power law in t; systematic = monotonic bound - estimate).
Per (beta, m, k, nu) and wave: class ratios R_c = T_c(fit)/T_c(pred) of the 8 parity classes (t5_analysis_l1.fit_components),
Sigma_+- = T_oo +- T_ss (W[v2,v2]: both present).  Corner systematic: the spread of R between (corrected) and (corrected + bound).
Usage: python3 t1_oddk_l1.py      Output: out_t1_oddk_l1.txt
"""
import os, sys, glob, itertools
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from cases import info
from tmoments4 import TensorMoments4, DATADIR
import t5_analysis_l1 as T5
from corner_est import moment_and_corner
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
KS = (3, 4, 5, 6, 7); NUS = (0.0, 1.0)
WAVES = ['ss', 'oo', 'v2v2', 'os']

def load_cases():
    moms = []
    for f in sorted(glob.glob(os.path.join(DATADIR, 'modes4_l1_*.npz'))):
        fc = os.path.join(HERE, 'cache', os.path.basename(f).replace('modes4_', 'modes4c_'))
        mom = TensorMoments4.load(fc if os.path.exists(fc) else f)
        mom.completed = os.path.exists(fc)
        moms.append(mom)
    return moms

def measure(mom, cL, cR, k, nu, logref):
    w = {(a, b): cL[a]*cR[b] for a in range(2) for b in range(2) if cL[a]*cR[b] != 0}
    raw, cor, hi, frac = {}, {}, {}, 0.0
    sg = (-1)**mom.J
    for q in T5.QS:
        v, e, bd, sh, sk = moment_and_corner(mom, k, nu, None, q, w)
        f = np.exp(sh - logref)*sg
        raw[q] = complex(v)*f; cor[q] = complex(v + e)*f
        ph = e/abs(e) if abs(e) > 0 else 1.0
        hi[q] = complex(v + ph*max(bd, abs(e)))*f
    return raw, cor, hi

def run_case(mom):
    D, J = mom.Delta, mom.J
    cmb = T5.combos(D, J)
    res = {}
    for k in KS:
        for nu in NUS:
            A = T5.design(D, J, k, nu)
            for name in WAVES:
                cL, cR = cmb[name]
                Tp = T5.Tpred(cL, cR, D, J)
                big = max(abs(sum(A[q][c]*Tp[c] for c in T5.CLASSES)) for q in T5.QS)
                logref = np.log(big)
                raw, cor, hi = measure(mom, cL, cR, k, nu, logref)
                entry = {}
                for lab, meas in (('raw', raw), ('cor', cor), ('hi', hi)):
                    R, sol, resid = T5.fit_components(A, meas, Tp, logref)
                    oo, ss = (1, 1, -1, -1), (1, -1, -1, 1)
                    pp, pm = Tp[oo] + Tp[ss], Tp[oo] - Tp[ss]
                    sc = max(abs(pp), abs(pm))
                    entry[lab] = dict(R=R, resid=resid, Sp=((sol[oo] + sol[ss])/pp if abs(pp) > 1e-6*sc else None),
                                      Sm=((sol[oo] - sol[ss])/pm if abs(pm) > 1e-6*sc else None))
                res[(k, nu, name)] = entry
    return res

def fmt(x, e=None):
    if x is None: return "   --   "
    return f"{x:.4f}" + (f"({abs(e)*1e4:.0f})" if e is not None else "")

if __name__ == "__main__":
    pr("# odd-k two-pair tensor Regge moments of <JJJJ> vs the factorised formula continued in k (fourpt/oddk/t1_oddk_l1.py)")
    pr("# entries: R = fit/pred (corner-corrected), in brackets the corner systematic |R(cor) - R(cor + bound)| in units of 1e-4; [raw] = no corner estimate")
    data = {}
    for mom in load_cases():
        m, beta, _ = info(mom.Delta, mom.J)
        data.setdefault(round(float(beta), 1), []).append((float(m), mom, run_case(mom)))
    summary = {}
    for bkey in sorted(data):
        rows = sorted(data[bkey], key=lambda r: r[0])
        pr(f"\n## beta ~ {bkey}")
        for k in KS:
            for nu in NUS:
                pr(f"# k={k} nu={nu}:   W[s,s] R_ss | W[o,o] R_oo | W[v2,v2] Sigma_+  Sigma_- | W[o,s] R_pf | LSQ resid (v2v2) | completed")
                for (m, mom, res) in rows:
                    def g(name, key):
                        e = res[(k, nu, name)]
                        if key in ('Sp', 'Sm'):
                            a, b, r = e['cor'][key], e['hi'][key], e['raw'][key]
                        else:
                            a, b, r = e['cor']['R'][key], e['hi']['R'][key], e['raw']['R'][key]
                        return a, (b - a if a is not None and b is not None else None), r
                    vals = [g('ss', (1, -1, -1, 1)), g('oo', (1, 1, -1, -1)), g('v2v2', 'Sp'), g('v2v2', 'Sm'), g('os', (1, 1, -1, 1))]
                    for nm, v in zip(('ss', 'oo', 'Sp', 'Sm', 'pf'), vals):
                        summary.setdefault((bkey, k, nu, nm), []).append((m, v[0], v[1]))
                    pr(f"  m={m:6.2f} (Delta,J)=({mom.Delta:.3f},{mom.J}): " + "  ".join(f"{fmt(a, e)} [{fmt(r)}]" for a, e, r in vals)
                       + f" | {res[(k, nu, 'v2v2')]['cor']['resid']:.1e} | {'N2' if mom.completed else '-'}")
    pr("\n# convergence in m: fits of R - 1 over the three largest masses (corrected values): A = (c1/m + c2/m^2 + c4/m^4) -> c1;"
       " B = (c2, c4, c6) -> c2;  and m^2 (R - 1) at the largest two masses")
    for key, lst in sorted(summary.items()):
        lst = sorted(lst)
        ms = [x[0] for x in lst]; Rs = [x[1] for x in lst]
        if any(r is None for r in Rs): continue
        f = T5.fits(ms, Rs)
        bkey, k, nu, nm = key
        pr(f"  beta~{bkey} k={k} nu={nu} {nm:3s}: R = " + " ".join(f"{r:.4f}" for r in Rs) +
           (f" | c1(A) = {f['A'][0]:+.3f}  c2(B) = {f['B'][0]:+.1f}" if 'A' in f else "") +
           f" | m^2(R-1) = {ms[-2]**2*(Rs[-2]-1):+.1f}, {ms[-1]**2*(Rs[-1]-1):+.1f}")
    open(os.path.join(HERE, 'out_t1_oddk_l1.txt'), 'w').write("\n".join(out) + "\n")
