"""
cl_points (RESULT_D_CLEANUP.md): thin wrapper around hc_points / oo_points (both UNCHANGED, imported read-only) that writes the direct
crystal points to NEW per-state caches with the cl_ prefix, so that no existing cache file is touched:
  mode 'hc' (default; = hc_points.py: n3 or n1 path via HC_METHOD, odd J -> odd 1x1 integrand, even J -> w0 2x2 integrand)
        -> cache/cl_pts_<kind>_J<J>_D<D>.pkl     (+ raw cache/cl_raw_J<J>_D<D>.pkl)
  mode 'oo' (= oo_points.py: even J, n3 path, 3x3 [even(2), odd(1)] p8 integrand, hc_points.setup(odd=True) forced)
        -> cache/cl_oopts_<kind>_J<J>_D<D>.pkl   (+ raw cache/cl_ooraw_J<J>_D<D>.pkl)
Method, seeds, point counts and record format are those of hc_points / oo_points (z2 n3 hybrid contour below/above the fold with
Z2_NUSPLIT = 5, or n1 adaptive re-centred contour; state and reference Delta_ref = 2e4 on common random numbers; one file per state).
CL_FWPATCH=1: also fw_frame.patch_conserved_clean (the structure basis of the FOLD_WINDOW fw rows; used for the (1048, 48) check).
Usage: OMP_NUM_THREADS=1 CL_MODE=oo [OO_M=10 OO_R=8] python3 cl_points.py J D x1,x2,...
       OMP_NUM_THREADS=1 CL_MODE=hc HC_METHOD=n3|n1 [HC_M=.. HC_MREF=.. HC_ADAPT=1] python3 cl_points.py J D x1,x2,...
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cl_guard     # noqa: E402,F401  (no write to a non-cl cache file)
MODE = os.environ.get('CL_MODE', 'hc')
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('Z2_NUSPLIT', '5')
if MODE == 'oo':
    os.environ['HC_M'] = os.environ.get('OO_M', '10')
    os.environ['HC_R'] = os.environ.get('OO_R', '8')
import warnings
warnings.filterwarnings('ignore')
import hc_points as HP      # noqa: E402
Z2 = HP.Z2
C = os.path.join(HERE, 'cache')
PRE = 'cl_oo' if MODE == 'oo' else 'cl_'


def fn(kind, J):
    return os.path.join(C, '%spts_%s_J%d_D%g.pkl' % (PRE, kind, J, HP.STATE[0]))


def rawfn(J):
    return os.path.join(C, '%sraw_J%d_D%g.pkl' % (PRE, J, HP.STATE[0]))


def _save(path, key, val):
    if os.path.basename(path).startswith('z2') and path.endswith('_raw.pkl'):
        path = rawfn(key[2])
    assert os.path.basename(path).startswith('cl_'), ('cl_points: refusing to write a non-cl cache', path)
    return HP._save0(path, key, val)


HP.fn = fn; HP.rawfn = rawfn; HP._save = _save
Z2.fn = fn; Z2._save = _save

if __name__ == '__main__':
    import numpy as np
    J = int(sys.argv[1]); D = float(sys.argv[2])
    HP.STATE[0] = D
    m = HP.mphys(D, J)
    a = sys.argv[3]
    nus = [float(v) for v in a[3:].split(',')] if a.startswith('nu=') else [float(v)*m for v in a.split(',')]
    meth = 'n3' if MODE == 'oo' else os.environ.get('HC_METHOD', 'n3')
    if MODE == 'oo':
        assert J % 2 == 0 and J >= 4
        HP.setup(odd=True)          # forced at even J, as oo_points
    else:
        HP.setup(odd=(J % 2 == 1))
    if os.environ.get('CL_FWPATCH', '0') == '1':     # structure basis of the fw rows (fw_points.setup: patch_conserved_clean)
        import l5_spin as S_
        import n1_pw as PW_
        import fw_frame as FW_
        FW_.patch_conserved_clean(S_)
        FW_.patch_conserved_clean(PW_.S)
    t0 = time.time()
    print(f"# cl_points [{MODE}, {meth}] J={J} Delta={D:g} m={m:.3f}: nu = {[round(v, 4) for v in nus]} "
          f"(x = {[round(v/m, 5) for v in nus]}); M={Z2.MQ} R={Z2.RQ} kinds {Z2.KINDS} NUSPLIT {Z2.NU_SPLIT} -> {fn('EE', J)}", flush=True)
    if meth == 'n1':
        HP.run_n1(J, nus, D, Z2.MQ, Z2.RQ)
    else:
        Z2.run(J, nus, [D, 20000.0], Z2.MQ, Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
