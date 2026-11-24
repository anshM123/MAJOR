"""
nf_points (NEARGAP_FILL.md): DIRECT crystal points -- highspin/hs_points.py with the cache moved to ngfill/cache (spinfunc/hc_points.py
unchanged, imported read-only; n1 adaptive re-centred contour, RQMC 2^11*2 points x 8 replicas, state and reference Delta_ref = 2e4 on
common random numbers; odd J -> parity-odd 1x1 integrand, even J -> w0 2x2 EE integrand + C3H).  Defaults HC_ADAPT=1, HC_M=HC_MREF=11,
L1_DEG=6 (= ext/point.py).
x-list: comma list of x = nu/m(Delta, J), 'gk15:a:b[:xcap]', 'dense' (x = 0.025 ... 0.95, the HIGHSPIN 38-node row), 'nu:v1,v2,..'
(explicit nu), or a mix joined by '+'.
Usage: python3 nf_points.py J Delta x-list
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import nf_guard  # noqa: F401  (first)
SF = os.path.join(HERE, '..', 'spinfunc')
HS = os.path.join(HERE, '..', 'highspin')
NFC = os.path.join(HERE, 'cache')
sys.path.insert(0, SF); sys.path.insert(0, HS); sys.path.insert(0, HERE)
import hs_points as HSP                                  # read-only: gk15 / parse_x only
DENSE = [round(0.025*i, 4) for i in range(1, 39)]         # 0.025 ... 0.95


def parse(spec, m):
    xs = []
    for part in spec.split('+'):
        if part == 'dense':
            xs += DENSE
        elif part.startswith('nu:'):
            xs += [float(v)/m for v in part[3:].split(',') if v]
        else:
            xs += HSP.parse_x(part)
    return xs


if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2])
    os.chdir(SF)
    os.environ['HC_CACHE'] = NFC
    for k, v in (('HC_METHOD', 'n1'), ('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')):
        os.environ.setdefault(k, v)
    import warnings
    warnings.filterwarnings('ignore')
    import hc_points as HP
    HP.STATE[0] = D
    m = HP.mphys(D, J)
    xs = parse(sys.argv[3], m)
    HP.setup(odd=(J % 2 == 1))
    t0 = time.time()
    print(f"# nf_points [n1] J={J} Delta={D:g} m={m:.3f} {len(xs)} nodes x={[round(x, 6) for x in xs]} -> {HP.fn('EE', J)}", flush=True)
    for x in xs:                                          # one node at a time (restartable; already-present nodes are skipped by hc_points)
        t1 = time.time()
        HP.run_n1(J, [x*m], D, HP.Z2.MQ, HP.Z2.RQ)
        print(f"  x={x:.6f} nu={x*m:.3f}: {time.time() - t1:.0f}s", flush=True)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
