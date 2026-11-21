"""Result Q, Stage 7 (RESULT_Q_STAGE7.md §3): minima of FF/N under Stage 6's non-covariant beyond-LR structures at NATURAL strengths |K| = 0.25,
0.5, 1, 2, 5, per gauge of the C4 derivative piece (H0, H2, H3, H4, SB of Stage 6, plus the Stage-7 extended-grid min-norm gauge XP if its
cache exists), on Stage 6's state set (exact + gap J <= 60 + region A with m <= 30 Lam, J <= 60, bhat < 30) and Dg = 100 / 300 / 1000.
Structures (stage6_stress): (a) Gamma_j = K [nu^2 - (nu - ij)^2], (b) Gamma_j = K (j/3)^2 nu^2  (c2 units; d(dB4) = 2 K u^3 s).
Also: the same with the C5 gauge fixed at G0 (immaterial, Stage 6).  Output: cache/stage7_kscan.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import stage6_model as M
import stage6_stress as ST
import stage4_hj as HJ4

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()
KS = (0.25, 0.5, 1.0, 2.0, 5.0)


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def s4(nu, chi):
    F4 = HJ4.numfuncs('EE3'); d = nu[1] - nu[0]
    R = {j: F4[j](nu.astype(complex)) for j in F4}; Nh = sum(R.values()).real
    mu1 = d*np.sum(chi*Nh*nu**2); M1 = d*np.sum(chi*sum(R[j]*(nu - 1j*j)**2 for j in R)).real
    sb = d*np.sum((chi*sum((j*j/9.)*R[j] for j in R)).real*nu**2)
    return mu1 - M1, sb


def main():
    t0 = time.time()
    c4, c5 = ST.coefficients()
    c4.pop('H1 tol 1e-6', None)                      # not grid-stable (Stage 6)
    for piece, key in (('c4v', 0), ('c4d', 1)):
        fn = os.path.join(CACHE, 'stage7_xprec_%s_d004_N3072.pkl' % piece)
        if os.path.exists(fn):
            d = pickle.load(open(fn, 'rb'))
            chi = np.array([float(mp.mpf(x)) for x in d['weights']]); nu = (np.arange(len(chi)) + 0.5)*0.04
            c4.setdefault('XP (Stage 7)', [None, None])[key] = s4(nu, chi)
    if 'XP (Stage 7)' in c4:
        if None in c4['XP (Stage 7)']:
            c4.pop('XP (Stage 7)')
        else:
            c4['XP (Stage 7)'] = tuple(c4['XP (Stage 7)'])
    pr("# stage7_kscan: min FF/N under non-covariant structures (a), (b) at |K| = " + ", ".join(map(str, KS)) + " (worst sign), per C4 gauge")
    for k, (v, dd) in c4.items():
        pr(f"   C4 {k:16s}: s_a, s_b value ({v[0]:+.4f}, {v[1]:+.4f}); derivative ({dd[0]:+.4f}, {dd[1]:+.4f})")
    from stage6_scan import XS, JCAP, BMAX
    g5 = 'G0 min-norm'
    res = {}
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        flat = M.flat_table(name, XS, JCAP)
        for Dg in (100, 300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = ST.states(name, Dg, flat, Lam)
            FF = np.zeros(len(rows)); NN = np.zeros(len(rows)); U = np.zeros(len(rows)); S = np.zeros(len(rows)); lab = []
            for i, (h, u, J, b2, b3, fw) in enumerate(rows):
                f, n, _, P = M.hyb(F, h, u, float(J), b2, b3, tab, Lam, fw=fw)
                FF[i], NN[i], U[i], S[i] = f, n, u, M.SIG[h]
                lab.append((h, J, 1/np.sqrt(u)))
            pr(f"\n== delta = {delta}, Dg = {Dg}: {len(rows)} cases; K = 0 min FF/N = {np.min(FF/NN):+.4f}")
            pr("   gauge             struct | " + " | ".join(f"|K|={K:<4}" for K in KS) + " | binding at |K| = 1")
            for g4 in c4:
                for st in ('a', 'b'):
                    ix = 0 if st == 'a' else 1
                    sv, sd, s5_ = c4[g4][0][ix], c4[g4][1][ix], c5[g5][ix]
                    dK = -c*2*U**3*sv/Lam**2 - a*(2*U**3*sd - S*U**4*s5_/Lam**2)
                    mins = []
                    for K in KS:
                        r = np.minimum((FF + K*dK)/NN, (FF - K*dK)/NN)
                        mins.append(np.min(r))
                    r1 = np.minimum((FF + dK)/NN, (FF - dK)/NN); i1 = int(np.argmin(r1))
                    res[(delta, Dg, g4, st)] = mins
                    pr(f"   {g4:16s}  ({st})   | " + " | ".join(f"{x:+.4f}" for x in mins) +
                       f" | {lab[i1][0]} J={lab[i1][1]} m={lab[i1][2]:.3f}Lam")
    pickle.dump(dict(c4=c4, res=res, KS=KS), open(os.path.join(CACHE, 'stage7_kscan.pkl'), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, 'stage7_kscan.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
