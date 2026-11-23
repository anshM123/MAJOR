"""
sn_points (SMALLNU_WINDOW.md): DIRECT crystal points at small nu (nu < 50) on a near-gap state's row -- spinfunc/hc_points.py,
unchanged and imported read-only, exactly as highspin/hs_points.py and neargap2/ng_points.py use it, with the cache moved to
smallnu/cache (hc_points' own HC_CACHE switch).  State and reference Delta_ref = 2e4 at the same nu on common random numbers, 8 RQMC
replicas, odd J -> parity-odd 1x1 integrand, even J -> w0 2x2 EE integrand + C3H.
Two contours (both hc_points paths):
  n3 : z2's n3 hybrid contour (nu_f = max(nu, 15) above 13, nu_f = 13 below; HYPOTHESES_CLOSE's small-nu protocol: Z2_NUSPLIT = 5,
       HC_M = 10, i.e. 2^10.2.8 points; proposal adapted on the state and re-used for the reference)
  n1 : the n1 real re-centred contour, adaptive Gaussian proposal (HC_M = HC_MREF = 11), as every highspin / neargap2 direct row.
Usage: python3 sn_points.py n3|n1 J Delta nu1,nu2,...  [subdir]
  (explicit nu, not x; 'subdir' = a separate cache under smallnu/cache/ for cross-checks, e.g. the same nu with the other contour.)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import sn_guard  # noqa: F401  (first)
SF = os.path.join(HERE, '..', 'spinfunc')
SNC = os.path.join(HERE, 'cache')
sys.path.insert(0, SF); sys.path.insert(0, HERE)

if __name__ == '__main__':
    meth, J, D = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    nus = [float(v) for v in sys.argv[4].split(',') if v]
    cdir = os.path.join(SNC, sys.argv[5]) if len(sys.argv) > 5 else SNC
    os.makedirs(cdir, exist_ok=True)
    os.chdir(SF)
    os.environ['HC_CACHE'] = cdir
    os.environ['HC_METHOD'] = meth
    if meth == 'n1':
        defaults = (('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5'))
    else:
        defaults = (('HC_M', '10'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5'))
    for k, v in defaults:
        os.environ.setdefault(k, v)
    import warnings
    warnings.filterwarnings('ignore')
    import hc_points as HP
    HP.STATE[0] = D
    m = HP.mphys(D, J)
    HP.setup(odd=(J % 2 == 1))
    t0 = time.time()
    print(f"# sn_points [{meth}] J={J} Delta={D:g} m={m:.3f} nu={nus} (x={[round(v/m, 5) for v in nus]}) -> {HP.fn('EE', J)}", flush=True)
    if meth == 'n1':
        HP.run_n1(J, nus, D, HP.Z2.MQ, HP.Z2.RQ)
    else:
        HP.Z2.run(J, nus, [D, 20000.0], HP.Z2.MQ, HP.Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
