"""n3fix/nf_jobs.py (N3C3H_FIX.md §5): the batch-6 external job list -> cache/jobs_batch6.txt (copied to ads4/ext/ by hand).
  G1  the n3 edge + overlap nodes of EDGE_R2's 8 gap-mass states (batch-3 *_n3a / *_n3b lines, same nu, same independent seeds), now with
      the patched evaluator (both kinds, EE and C3H; method n3m<M3>);
  G2  the J = 46, 50 ++ dense-grid recheck (EDGE_R2 §4 (d), ‡): n1 at the midpoints x = 0.3125, 0.3375, ..., < a = 0.8 of the 0.025 bulk
      grid (completing a 0.0125 grid on [0.3, 0.8)), at M + 2 = 13 (method n1m13); kept >= 0.004 from the reference zero.
Usage: python3 nf_jobs.py [M3]"""
import nf_guard  # noqa: F401
import os, sys, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.abspath(os.path.join(HERE, '..', '..', 'ext'))
M3 = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STATES = ('P1_J26_m1.026', 'P1_J30_m1.030', 'P1_J34_m1.033', 'P2_J33_m1.033', 'P2_J38_m1.037', 'P2_J42_m1.041', 'P2_J46_m1.045',
          'P2_J50_m1.049')
lines = open(os.path.join(EXT, 'jobs_batch3.txt')).read().splitlines()
hdr = {}
for ln in lines:
    mm = re.match(r'# (P\d_J\d+_m[\d.]+): (.*)', ln)
    if mm:
        hdr[mm.group(1)] = mm.group(2)
out = ["# ext/jobs_batch6.txt -- N3C3H_FIX.md §5: EDGE-R edge completion with the FIXED n3 C3H evaluator + the J = 46, 50 ++ recheck.",
       "# RUN WITH THE PATCHED RUNNER (NOT run_ext.py / point.py, which import the unpatched code):",
       "#     cd ads4/ext && python3 point_patched.py run 16 jobs_batch6.txt      (then python3 pack.py as usual)",
       "# point_patched.py installs ext/patch_n3c3h.py after hc_points.setup(); even-J jobs also save OUTDIR/basis.pkl (local <TTO>",
       "# structure basis), so no BATCH5_FOLD §1.2 fit is needed for these records (exact map: N3C3H_FIX.md §5).",
       "# format as jobs.txt: TAG METHOD J DELTA SEED nu1,nu2,...; METHOD n3m<K> / n1m<K> = n3 / n1 with HC_M = HC_MREF = K.",
       "# n3 seeds as batch 3 (two independent streams per node: 6000 + 8k and + 4096)."]
nG1 = nG2 = 0
for st in STATES:
    out.append(f"# G1 {st}: {hdr.get(st, '')}")
    for ln in lines:
        f = ln.split()
        if f and f[0] in (st + '_n3a', st + '_n3b'):
            f[0] = 'B6_' + f[0]; f[1] = 'n3m%d' % M3
            out.append('  '.join(f)); nG1 += len(f[5].split(','))
for st, J, D in (('P2_J46_m1.045', 46, 1046.0), ('P2_J50_m1.049', 50, 1050.0)):
    al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
    x0 = float(re.search(r'reference zero x0 = ([\d.]+)', hdr[st]).group(1))
    xs = [0.3125 + 0.025*k for k in range(20)]
    xs = [x if abs(x - x0) >= 0.004 else x0 + 0.004*np.sign(x - x0) for x in xs]
    out.append(f"# G2 {st} ++ recheck: n1 M = 13 at the 0.025-grid midpoints x = {xs[0]:.4f} .. {xs[-1]:.4f} (reference zero {x0}); "
               f"m = {m:.6f}")
    out.append('  '.join(['B6_%s_pp_n1' % st, 'n1m13', str(J), '%.6f' % D, '4242', ','.join('%.6f' % (x*m) for x in xs)]))
    nG2 += len(xs)
out.append(f"# totals: G1 {nG1} n3 jobs, G2 {nG2} n1 jobs")
open(os.path.join(HERE, 'cache', 'jobs_batch6.txt'), 'w').write('\n'.join(out) + '\n')
print(f"G1 {nG1} n3 jobs, G2 {nG2} n1 jobs -> cache/jobs_batch6.txt")
