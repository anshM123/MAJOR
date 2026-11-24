"""
gu1_tol (GAP_UNIFORM.md §3): chart-A tolerance of a near-gap state of Omega^EE on the beyond-hybrid correction Gamma (c2 units),
as a function of the gap.  Read-only import of poserr/large_m_map.state (-> spinfunc/ee1000_common.hyb); writes only to gapuni/cache.

For a state (J, mu = m/Lam, h) the chain's correction is (ee1000_common.dF_state)
    dF = sum_n [ p2_n env(h, x_n, 2) Gamma_C2(x_n) + u p3_n env(h, x_n, 3) Gamma_C3H(x_n) ] / m^2,   x_n = p_n/mu,
with p2_n = psi2(p_n) w_n, p3_n = psi3(p_n) w_n (Lam-independent at fixed (J, mu)), and the state value lambda = (base + dF)/N_s.
Hence  |dF|/N_s <= sup|Gamma| * W,   W = A/(m^2 N_s),  A = sum(|p2| env2 + u |p3| env3)   (A Lam-independent at fixed (J, mu, h)).
Tolerances (sign-free, uniform |Gamma| <= G on every node):
    G_cert = rig/W   with rig = (N_f/5 - LemmaC)/N_s    (certified flat floor minus the Lemma-C curvature bound)
    G_hyb  = lam/W   with lam = base/N_s                (hybrid value, quadrature)
Lemma T (GAP_UNIFORM §3.1): A, N_f are Lam-independent and LemmaC decreases in Lam, so G_cert(Lam) >= G_cert(Lam0) (Lam/Lam0)^2 *
(N_s(Lam)/N_s(Lam0)) ... checked here numerically: the columns W*Lam^2 and G/(Lam/Lam0)^2 at 1000 / 2000 / 5000.
Usage: python3 gu1_tol.py [Jmax]  -> cache/gu1_tol.txt, cache/gu1_tol.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gu_guard  # noqa: F401  (write guard: only gapuni/cache)
sys.path.insert(0, os.path.join(HERE, '..', 'poserr'))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import large_m_map as LM
E = LM.E

C = os.path.join(HERE, 'cache')
OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def main():
    Jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    DGS = (1000.0, 2000.0, 5000.0)
    L0 = E.lam_of(1000.0)
    Js = [J for J in range(0, Jmax + 1) if not (J % 2 == 1 and J < 4)]
    res = []
    t0 = time.time()
    pr("# chart-A tolerances: G_cert = rig/W, G_hyb = lam/W (c2), W = cost of uniform |Gamma| = 1 c2 (N_s units); A = W m^2 N_s")
    pr(f"{'J':>3} {'m/Lam':>8} {'h':>3} | {'lam1000':>8} {'rig1000':>8} | {'W*Lam^2: 1000':>13} {'2000':>9} {'5000':>9} | "
       f"{'Gcert1000':>10} {'Gc2000/r':>9} {'Gc5000/r':>9} | {'Ghyb1000':>9} {'Gh2000/r':>9} {'Gh5000/r':>9}")
    for J in Js:
        mg = E.mgap(1000.0, J)
        mus = sorted(set([round(mg*1.000001, 6)] + [v for v in (1.02, 1.05, 1.1, 1.3) if v > mg]))
        for mu in mus:
            rows = {}
            for Dg in DGS:
                s = LM.state(Dg, J, mu)
                rows[Dg] = s
            for h in (['++'] if J < 4 else (['+-'] if J % 2 else ['++', '+-'])):
                rec = dict(J=J, mu=mu, h=h)
                for Dg in DGS:
                    s = rows[Dg][h]; Lm = E.lam_of(Dg)
                    W = s['W_lo'] + s['W_hi']
                    rec[str(int(Dg))] = dict(Lam=float(Lm), lam=s['lam'], rig=s['rig'], W=W, WL2=W*Lm**2, Ns=s['Ns'], Nf=s['Nf'],
                                             Gcert=s['rig']/W, Ghyb=s['lam']/W, Wlo=s['W_lo'])
                res.append(rec)
                a, b, c = rec['1000'], rec['2000'], rec['5000']
                r2 = (b['Lam']/L0)**2; r5 = (c['Lam']/L0)**2
                pr(f"{J:>3} {mu:8.5f} {h:>3} | {a['lam']:8.4f} {a['rig']:+8.4f} | {a['WL2']:13.5e} {b['WL2']/a['WL2']:9.5f} {c['WL2']/a['WL2']:9.5f} | "
                   f"{a['Gcert']:10.3e} {b['Gcert']/r2/a['Gcert'] if a['Gcert'] > 0 else float('nan'):9.4f} "
                   f"{c['Gcert']/r5/a['Gcert'] if a['Gcert'] > 0 else float('nan'):9.4f} | "
                   f"{a['Ghyb']:9.3e} {b['Ghyb']/r2/a['Ghyb']:9.4f} {c['Ghyb']/r5/a['Ghyb']:9.4f}")
        json.dump(res, open(os.path.join(C, 'gu1_tol.json'), 'w'))
    pr(f"# done {len(res)} states, {time.time() - t0:.0f}s")
    open(os.path.join(C, 'gu1_tol.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
