"""
t3_q3check (RESULT_D_CONDITIONAL.md §3): evaluate NF02 and re-weighted (psi2, psi3) on the near-gap states WITH the actually computed
graviton corrections (q3 model M of S1_S4_FIXES: K2 rows, aligned off-diagonal, c_off model, fixed C3; + hybrid extra states), i.e. the
r3_redesign blocks, instead of r5's universal-g rows.  Prints per J the min lambda/N_s and the b-hat of the worst state.
Usage: T3_DG=300 python3 t3_q3check.py <weight pickles in cache/ (keys 'v' or res[dw]['v'])> -> cache/t3_q3check_<Dg>.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('R3_DG', os.environ.get('T3_DG', '300'))
import numpy as np
import r3_redesign as R3

CACHE = os.path.join(HERE, 'cache')


def main():
    Dg = int(os.environ.get('T3_DG', '300'))
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    emods, omods = R3.build_models()
    blocks, Lm = R3.state_rows(Dg, emods, omods)
    pr(f"# t3_q3check Delta_g = {Dg}: {len(blocks)} blocks (q3 model M, corrected even J {sorted(emods)}, odd J {sorted(omods)}; "
       f"hybrid extra m/Lam {os.environ.get('Q3_EXTRA_M')})  [{time.time() - t0:.0f}s]")
    ws = [('NF02', R3.nf_vector())]
    for a in sys.argv[1:]:
        d = pickle.load(open(os.path.join(CACHE, a), 'rb'))
        ws.append((a, d['v'] if 'v' in d else d['res'][0.1]['v']))
    for name, v in ws:
        ls = [(R3.lam(b, v)[0], b) for b in blocks]
        neg = [(l, b) for l, b in ls if l < 0]
        pr(f"\n## {name}: min lambda/N_s {min(l for l, b in ls):+.4f}; negative states {len(neg)} of {len(ls)}"
           + (f"; smallest negative b-hat {min((2*b['J'] + 1)/b['mfl'] for l, b in neg):.1f}" if neg else ""))
        for J in sorted(set(b['J'] for b in blocks)):
            lj = [(l, b) for l, b in ls if b['J'] == J]
            l, b = min(lj, key=lambda t: t[0])
            pr(f"    J={J:2d}: min {l:+.3f} at m/Lam {b['mfl']:.3f} (b-hat {(2*J + 1)/b['mfl']:.1f}, {b['DF']}, corrected {b['cov']})")
    open(os.path.join(CACHE, f't3_q3check_{Dg}.txt'), 'w').write(buf.getvalue())
    pr(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()
