"""
b4_basischk (BATCH4_FOLD.md §1.2; a copy of b5fold/b5_basischk.py, rows 2000 + J): does the structure-basis mapping (b4_basis.py) reproduce the LOCAL points in the chain's own
quantity?  Point set 'fixloc' = 'fix' with the external (mapped) records at the locally recomputed nodes replaced by the local ones
(state and reference, every kind computed locally); Gamma (c2) of the injected own row at those nodes, both sets, every entry
(++, +-, c_off = m Gamma_off; C2 and, where computed, C3H), with the jackknife sigma of the 'fix' value.
Usage: python3 b4_basischk.py J,J,...   -> cache/b4_basischk.txt
"""
import os, sys, glob, pickle, shutil
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_core as C
import numpy as np
import b4_state  # noqa: F401  (only for its K2-basis loader patch; its __main__ is not run)
H = C.H
PTS = C.PTS; LOC = os.path.join(C.B4C, 'loc')
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


if __name__ == '__main__':
    Js = [int(v) for v in sys.argv[1].split(',')]
    os.makedirs(os.path.join(PTS, 'fixloc'), exist_ok=True)
    for f in glob.glob(os.path.join(PTS, 'fix', '*.pkl')):
        shutil.copyfile(f, os.path.join(PTS, 'fixloc', os.path.basename(f)))
    for J in Js:
        D = 2000.0 + J; m = C.mphys(D, J)
        loc = {}
        for kind in ('EE', 'C3H'):
            fl = os.path.join(LOC, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))
            if os.path.exists(fl):
                PL = pickle.load(open(fl, 'rb'))
                P = pickle.load(open(os.path.join(PTS, 'fix', 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'rb'))
                for k in PL:
                    assert k in P, k
                    P[k] = PL[k]
                    loc.setdefault(kind, set()).add(round(k[1]/m, 6))
                pickle.dump(P, open(os.path.join(PTS, 'fixloc', 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'wb'))
        g = {}
        for st in ('fix', 'fixloc'):
            C.SET[0] = st; C._CACHE.clear()
            g[st] = H.grid(J, inject=[D])[D]
        C.SET[0] = 'fix'; C._CACHE.clear()
        reps = [H.grid(J, inject=[D], jk=j)[D] for j in range(8)]
        for x in sorted(loc.get('EE', set())):
            i = int(np.argmin(np.abs(g['fix']['xs'] - x))); i2 = int(np.argmin(np.abs(g['fixloc']['xs'] - x)))
            s = f"J={J} x={x:.4f}:"
            for kd in ('C2', 'C3H'):
                if (kd == 'C2' and x not in loc['EE']) or (kd == 'C3H' and x not in loc.get('C3H', set())):
                    continue
                for h in ('++', '+-', 'off'):
                    k = (kd, h); sc = m if h == 'off' else 1.0
                    a, b = g['fix'][k][i]*sc, g['fixloc'][k][i2]*sc
                    vj = np.array([r[k][int(np.argmin(np.abs(r['xs'] - x)))] for r in reps])*sc
                    e = float(np.sqrt(7/8*((vj - vj.mean())**2).sum()))
                    s += f" | {kd}{h}: mapped {a:+9.1f} local {b:+9.1f} (diff {b - a:+7.2f}; sigma {e:6.1f})"
            pr(s)
    open(os.path.join(C.B4C, 'b4_basischk.txt'), 'a').write("\n".join(lines) + "\n")
