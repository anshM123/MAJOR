"""
ext/point_patched.py -- ext/point.py with the n3 C3H fix (N3C3H_FIX.md) applied, for batch 6 and later.
Same physics and output format as ext/point.py (hc_points unchanged, imported read-only), plus:
  * ext/patch_n3c3h.py is installed after hc_points.setup(): (a) o_H's second fractional power is replaced by the exact physical-sheet
    identity n4 = n2 X40/X20 (relative sign of the two 2<->4 terms fixed to +1), (b) n3's graviton evaluators take principal logs at
    their saddle (n1's convention) instead of set_branch's W_REF continuation.  n1 values are unchanged by (a) and (b) (checked to
    1e-14 point-wise, N3C3H_FIX.md §2); n3 EE values are unchanged to 1e-9; n3 C3H STATE values change (the original ones are wrong).
  * even J: the local <TTO> structure basis is saved next to the points (OUTDIR/basis.pkl: mons, sol = l5_spin.conserved(D, J)[1] and
    solt = conserved(3 - D, J)[1] for D = Delta and 2e4), so that the BATCH5_FOLD §1.2 repair is exact:
       F_ext[a] = sum_b F_loc[b] Mf[b, a],  Mf = pinv(sol_loc) @ sol_ext   (and Mg likewise from solt),  E_loc = Mf^-T E_ext Mg^-1.
  * NU may be a comma list (one process, several nodes).
Usage: python3 point_patched.py METH J DELTA NU[,NU..] SEED OUTDIR     (METH = n1 | n3 | n1mK | n3mK: K = log2 samples, state and ref)
       python3 point_patched.py run [N] [jobs_batch6.txt]   -- the run_ext.py driver (its source, exec'd with 'point.py' replaced by
                                                               this file; same out/ layout, done markers, restartability)
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
if sys.argv[1] == 'run':
    _H = os.path.dirname(os.path.abspath(__file__))
    _src = open(os.path.join(_H, 'run_ext.py')).read()
    assert _src.count("'point.py'") == 1
    sys.argv = [os.path.join(_H, 'run_ext.py')] + sys.argv[2:]
    _g = sys.modules['__main__'].__dict__             # into __main__ itself: Pool pickles `one` by reference to __main__
    _g['__file__'] = os.path.join(_H, 'run_ext.py')
    exec(compile(_src.replace("'point.py'", "'point_patched.py'"), 'run_ext_patched', 'exec'), _g)
    sys.exit(0)
meth, J, D, nus, seed, out = (sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), [float(v) for v in sys.argv[4].split(',')],
                              int(sys.argv[5]), sys.argv[6])
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'code', 'spinfunc')
import re
_mm = re.fullmatch(r'(n1|n3)(?:m(\d+))?', meth)          # METH may carry the log2 sample count: n1m13 = n1 with HC_M = HC_MREF = 13
assert _mm, meth
meth = _mm.group(1)
if _mm.group(2):
    os.environ['HC_M'] = _mm.group(2); os.environ['HC_MREF'] = _mm.group(2)
os.makedirs(out, exist_ok=True)
os.environ['HC_CACHE'] = os.path.abspath(out)
os.environ['HC_METHOD'] = meth
os.environ['Z2_SEED'] = str(seed)
defaults = ((('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')) if meth == 'n1'
            else (('HC_M', '10'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')))
for k, v in defaults:
    os.environ.setdefault(k, v)
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
sys.path.insert(0, SF)
sys.path.append(HERE)                     # patch_n3c3h (after spinfunc: no name clash)
os.chdir(SF)
import warnings
warnings.filterwarnings('ignore')
import hc_points as HP
import patch_n3c3h as PX
HP.STATE[0] = D
m = HP.mphys(D, J)
HP.setup(odd=(J % 2 == 1))
st = PX.install()
t0 = time.time()
print(f"# point_patched [{meth}] J={J} Delta={D:g} m={m:.3f} nu={nus} x={[round(v/m, 6) for v in nus]} seed={seed} patch={st}", flush=True)
if J % 2 == 0:
    import l5_spin as S
    bas = {'mons': S.conserved(D, J)[0]}
    for DD in (D, 20000.0):
        bas[(DD, 'sol')] = S.conserved(DD, J)[1]; bas[(DD, 'solt')] = S.conserved(3.0 - DD, J)[1]
    with open(os.path.join(out, 'basis.pkl'), 'wb') as f:
        pickle.dump(bas, f)
if meth == 'n1':
    HP.run_n1(J, nus, D, HP.Z2.MQ, HP.Z2.RQ)
else:
    HP.Z2.run(J, nus, [D, 20000.0], HP.Z2.MQ, HP.Z2.RQ)
print(f"# done [{time.time() - t0:.0f}s]", flush=True)
