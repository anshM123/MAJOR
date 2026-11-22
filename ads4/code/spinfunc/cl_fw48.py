"""
cl_fw48 (RESULT_D_CLEANUP.md, task 3): the J = 48, Delta_g = 1000 chain of FOLD_WINDOW_LARGEJ §4 (fw_chain48.py, run UNCHANGED in-process:
ee1000 shape chain K = 4, x in [0.05, 0.8], rows 1048, hybrid states, sm_curvfix kappa, fixed-engine fw_k2 rows, generalised helicity
basis) with this agent's DIRECT n1 nodes on (1048, 48) (cl_points.py CL_FWPATCH=1 -> cache/cl_pts_*_J48_D1048.pkl) served through
fw_gamma's direct-node path (fw_gamma.fwpts, replaced right after fw_chain48 imports fw_gamma (import order unchanged); a direct node replaces the K2 node at the same nu; dir estimator).
Env: CL48_SKIP (x of failed nodes left out; default 0.475, whose 2e4 reference has E00/E01 = 82 against 1.3-1.4 at every other node), CL48_XMIN / CL48_XMAX (x-range of the direct nodes used; default all), CL48_JK = i (replica i left out), FW_NODIRECT=1 (pure K2 =
the published chain), CL48_TAG (output tag), CL48_DROPK2 (x of K2 nodes removed from the row, state + reference: 0.475, where the 2e4 reference
matrix is nearly singular (E01 ~ 0) and the chain's off-diagonal becomes non-finite once direct nodes change the helicity basis).   -> cache/cl_fw48_<tag>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cl_guard     # noqa: E402,F401
import numpy as np
C = os.path.join(HERE, 'cache')
TAG = os.environ.get('CL48_TAG', 'x')
XMIN = float(os.environ.get('CL48_XMIN', '0')); XMAX = float(os.environ.get('CL48_XMAX', '9'))
SKIP = [float(v) for v in os.environ.get('CL48_SKIP', '0.475').split(',') if v]   # failed reference point(s), see RESULT_D_CLEANUP §3
JK = None if os.environ.get('CL48_JK', '') == '' else int(os.environ['CL48_JK'])
D48 = 1048.0
os.environ.setdefault('FW_J', '48')


def mphys(D, J):
    a = D - J - 2.0
    return np.sqrt(a*(a + 2*J + 1))


def fwpts(kind, J, D=None):
    if J != 48:
        return {}
    f = os.path.join(C, 'cl_pts_%s_J48_D1048.pkl' % kind)
    P = pickle.load(open(f, 'rb')) if os.path.exists(f) else {}
    m = mphys(D48, 48)
    keep = {nu for (DD, nu), r in P.items() if DD == D48 and XMIN - 1e-9 <= nu/m <= XMAX + 1e-9
            and float(np.max(np.abs(r['err']))/np.max(np.abs(r['E']))) <= 0.1 and not any(abs(nu/m - v) < 1e-6 for v in SKIP)}
    return {k: v for k, v in P.items() if k[1] in keep}


DROPK2 = [float(v) for v in os.environ.get('CL48_DROPK2', '').split(',') if v]


def _cl_patch(FG):
    """FG.fwpts -> this agent's direct nodes; optionally drop the K2 nodes at x in CL48_DROPK2 (state and reference) from the row."""
    FG.fwpts = fwpts
    _gb = FG.gbasis

    def gbasis(J, D, R):              # helicity basis from the pure K2 order-0 records (as hc_chain: basis computed without injected nodes;
        if J != 48:                   # the direct nodes replace K2 records at the same key, so the K2 set is reloaded)
            return _gb(J, D, R)
        return _gb(J, D, FG.Q._strip(dict(FG.fwk2('EE', J)), J))
    FG.gbasis = gbasis
    if DROPK2:
        _k2 = FG.fwk2
        m = mphys(D48, 48)

        def fwk2(kind, J):
            R = _k2(kind, J)
            return {k: v for k, v in R.items() if not any(abs(k[1]/m - x) < 1e-6 for x in DROPK2)} if J == 48 else R
        FG.fwk2 = fwk2


path = os.path.join(HERE, 'fw_chain48.py')
src = open(path).read()
for old, new in (("FG.CUR[J0] = dict(jk=None)", "FG.CUR[J0] = dict(jk=_CL_JK)"),
                 ("import fw_gamma as FG\n", "import fw_gamma as FG\n_CL_PATCH(FG)\n"),
                 ("open(os.path.join(C, 'fw_chain%d_%s.txt' % (J0, VAR)), 'w')", "open(os.path.join(C, 'cl_fw48_%s.txt' % _CL_TAG), 'w')")):
    assert src.count(old) == 1, old
    src = src.replace(old, new)
g = dict(__name__='__main__', __file__=path, __builtins__=__builtins__, _CL_JK=JK, _CL_TAG=TAG, _CL_PATCH=_cl_patch)
nd = sorted(k[1]/mphys(D48, 48) for k in fwpts('EE', 48) if k[0] == D48)
print(f"# cl_fw48 [{TAG}]: direct nodes used (x): {[round(x, 4) for x in nd]} (FW_NODIRECT={os.environ.get('FW_NODIRECT', '0')}, jk {JK})",
      flush=True)
exec(compile(src, path, 'exec'), g)
