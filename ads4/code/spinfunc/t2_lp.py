"""
t2_lp (RESULT_D_CONDITIONAL.md §3): trade-off "covered window b-hat <= b_w  vs  bound", re-weighting (psi2, psi3) inside NF02's LP family
with r5's universal-g AdS rows (LARGEJ_OFFDIAG.md §3), for windows not run by r5.  Same LP as `r5_universal.py lp` (r3_redesign.solve),
only the output names differ (cache/t2_lp_<Dg>_b<beta>_dw<dw>_<obj>.pkl/.txt), and the exact 16 L_A bound is printed next to the log form.
Usage: T2_DG=100 T2_BMAX=0.2 T2_DW=0.0,0.1 T2_OBJ=A|bound python3 t2_lp.py
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import r5_universal as R5
import r3_redesign as R3
from t1_band import LA_exact

CACHE = os.path.join(HERE, 'cache')


def main():
    Dg = int(os.environ.get('T2_DG', '100')); bmax = float(os.environ.get('T2_BMAX', '0.2'))
    dws = [float(v) for v in os.environ.get('T2_DW', '0.0,0.1').split(',')]; obj = os.environ.get('T2_OBJ', 'A')
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    G = R5.NL.NFGrid(R5.D, R5.R); v_nf = R3.nf_vector()
    mfls = [float(v) for v in os.environ.get('R5_MFLS', '1.02,1.05,1.1,1.2,1.3,1.5,2.0,3.0,5.0').split(',')]
    bhs = np.unique(np.round(np.geomspace(10 if Dg > 150 else 5, 2*Lm, 60), 3))
    sts = R5.grid_states(Dg, Lm, bmax, mfls, bhs)
    blocks, _ = R5.blocks_for(Dg, sts)
    bh_max = max(b['bh'] for b in blocks)
    pr(f"# t2_lp Delta_g = {Dg} (Lam {Lm:.2f}); r5 universal-g rows, {len(blocks)} states, beta <= {bmax} (b-hat <= {bh_max:.1f}); obj {obj}")
    L0, _ = LA_exact(v_nf, Lm)
    pr(f"  NF02: bound exact {240*np.sqrt(L0)/Lm**4:.4e}")
    R3.summary(blocks, v_nf, pr, 'NF02')
    for dw in dws:
        st, v, info = R3.solve(G, blocks, 0.2, dw, Lm, obj, v0=v_nf)
        if st != 0:
            pr(f"  dw = {dw}: status {st} ({info.get('msg')})"); continue
        L16, _ = LA_exact(v, Lm)
        pr(f"  dw = {dw}: A = {info['A']:.2f}, B_AdS = {info['B_AdS']:.2f}; bound at {Dg}: exact {240*np.sqrt(L16)/Lm**4:.4e}, "
           f"log form {R3.bound(info['A'], info['B_AdS'], Dg):.4e}")
        R3.summary(blocks, v, pr, f'redesign dw={dw}')
        pickle.dump(dict(v=v, info=info, Dg=Dg, bmax=bmax, bh_max=bh_max, obj=obj),
                    open(os.path.join(CACHE, f't2_lp_{Dg}_b{bmax}_dw{dw}_{obj}.pkl'), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(CACHE, f't2_lp_{Dg}_b{bmax}_{obj}.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
