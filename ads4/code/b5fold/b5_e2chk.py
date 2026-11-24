"""
b5_e2chk (BATCH5_FOLD.md §1.2): cross-check of the structure-basis repair with EDGE_R2's joint diagonalisation (edger2/e2_basis.py,
imported READ-ONLY; its e2_guard import resolves to b5_guard): the external records (set 'all') mapped by e2_basis.fix into set
'e2fix', the helicity-label permutation chosen by EDGE_R2's rule (smallest max over (kind, h) of rms Gamma on 0.05 <= x < 0.6), and the
resulting own-row Gamma compared node by node with this note's exact mapping (set 'fix', b5_basis.py).
Usage: python3 b5_e2chk.py J,J,...   -> cache/pts/e2fix/*, cache/b5_e2chk.txt / .json
"""
import os, sys, glob, json, shutil
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_core as C
import numpy as np
import b5_state  # noqa: F401  (K2-basis loader patch only)
sys.path.insert(0, os.path.join(HERE, '..', 'edger2'))
import e2_basis as EB            # read-only
H = C.H; HN = H.HN; NC = H.NC; Q = NC.Q; W = NC.W
PTS = C.PTS
PERMS = [((0, 1), (0, 1)), ((1, 0), (0, 1)), ((0, 1), (1, 0)), ((1, 0), (1, 0))]
lines, res = [], {}


def pr(s=""):
    print(s, flush=True); lines.append(s)


def row(J, D, st):
    C.SET[0] = st; C._CACHE.clear()
    try:
        return H.grid(J, inject=[D])[D]
    finally:
        C.SET[0] = 'fix'; C._CACHE.clear()


if __name__ == '__main__':
    Js = [int(v) for v in sys.argv[1].split(',')]
    dst = os.path.join(PTS, 'e2fix'); os.makedirs(dst, exist_ok=True)
    for f in glob.glob(os.path.join(PTS, 'all', '*.pkl')):
        if not os.path.exists(os.path.join(dst, os.path.basename(f))):
            shutil.copyfile(f, os.path.join(dst, os.path.basename(f)))
    keys = [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-', 'off')]
    for J in Js:
        D = 1000.0 + J; m = C.mphys(D, J)
        Re = Q._strip(W.ee_even(J), J)
        HN.RESID.clear(); HN.set_resid(Re, J, 'EE')
        Tloc = {DD: Q.A.L.basis(J, DD, Re)['T'] for DD in (D, 20000.0)}
        rf = row(J, D, 'fix')
        sc = []
        for pm in PERMS:
            rep = EB.fix(J, D, HN, Tloc, perm=pm, src=os.path.join(PTS, 'all'), dst=dst)
            r = row(J, D, 'e2fix'); s = (r['xs'] >= 0.05) & (r['xs'] < 0.6)
            rms = max(float(np.sqrt(np.mean(np.asarray(r[k])[s]**2))) for k in keys if k[1] != 'off')
            sc.append((rms, pm))
        sc.sort(key=lambda t: t[0])
        pm = sc[0][1]
        rep = EB.fix(J, D, HN, Tloc, perm=pm, src=os.path.join(PTS, 'all'), dst=dst)
        r = row(J, D, 'e2fix')
        com = [(i, int(np.argmin(np.abs(rf['xs'] - x)))) for i, x in enumerate(r['xs']) if np.min(np.abs(rf['xs'] - x)) < 1e-7 and x < 0.94]
        dmax = {}
        for k in keys:
            scl = m if k[1] == 'off' else 1.0
            d = np.array([abs(r[k][i] - rf[k][j])*scl for i, j in com])
            dmax[str(k)] = (float(np.median(d)), float(d.max()))
        o = dict(perm=pm, rms_scores=[(float(a), b) for a, b in sc],
                 resid={f'{DD:g}': dict(off_EE_med=float(np.median(q['off_EE'])), off_EE_max=float(q['off_EE'].max()),
                                         off_C3H_med=float(np.median(q['off_C3H'])), off_C3H_max=float(q['off_C3H'].max())) for DD, q in rep.items()},
                 dGamma=dmax, n=len(com))
        res[J] = o
        pr(f"J={J}: e2_basis perm {pm} (rms scores {[round(a, 1) for a, _ in sc]}); joint-diag residual off/||C||: " +
           "; ".join(f"Delta={DD}: EE med {v['off_EE_med']:.1e} max {v['off_EE_max']:.1e}, C3H med {v['off_C3H_med']:.1e} max {v['off_C3H_max']:.1e}"
                     for DD, v in o['resid'].items()))
        pr(f"      Gamma(e2_basis) - Gamma(exact local mapping) over {len(com)} nodes x < 0.94 (median / max, c2): " +
           ", ".join(f"{k[0]}{k[1]} {a:.1f}/{b:.1f}" for k, (a, b) in ((eval(kk), v) for kk, v in dmax.items())))
    open(os.path.join(C.B5C, 'b5_e2chk.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.B5C, 'b5_e2chk.json'), 'w'), indent=1, default=str)
