"""
cl_chain (RESULT_D_CLEANUP.md): runs the UNCHANGED hc_chain.py (modes ee1000 / w2 / z4) or oo_chain.py in-process, with three
redirections only (string-level, each asserted to occur exactly once in the source):
  * direct points: hc_nodes.pts(kind, J, D) and oo_chain.oo_pts(kind, J, D) return the union of the existing per-state cache
    (hc_pts_* / oo_pts_*, read only) and this agent's per-state cache (cl_pts_* / cl_oopts_*, written by cl_points.py); a reference key
    (Delta = 2e4, nu) present in both files is an error (one owner per reference point);
  * failed direct nodes (state max|err|/max|E| > CL_MAXREL = 0.1, e.g. x = 0.9 at Delta = 336, x = 0.95 at Delta = 200/400) are dropped;
  * exact-scalar caches: hc_exact.pkl -> cl_exact.pkl, oo_exact.pkl -> cl_ooexact.pkl (initialised as copies of the originals, so the
    control-variate values of already computed states are identical; new states are appended to the cl_ copies only);
  * outputs: cache/hc_chain_<mode>_<tag>.* -> cache/cl_chain_<mode>_<tag>.*, cache/oo_chain_<tag>.* -> cache/cl_oochain_<tag>.*.
cl_guard refuses any write to a non-cl cache file.
Usage: [env as for hc_chain.py] python3 cl_chain.py hc z4 300      |   [env as for oo_chain.py] python3 cl_chain.py oo
"""
import os, sys, shutil, pickle
import numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cl_guard     # noqa: E402,F401
C = os.path.join(HERE, 'cache')
DREF = 20000.0
for a, b in (('hc_exact.pkl', 'cl_exact.pkl'), ('oo_exact.pkl', 'cl_ooexact.pkl')):
    if not os.path.exists(os.path.join(C, b)):
        shutil.copyfile(os.path.join(C, a), os.path.join(C, b))


MAXREL = float(os.environ.get('CL_MAXREL', '0.1'))     # state points with max|err|/max|E| above this are dropped (failed RQMC nodes)
DROPPED = set()


def merged(files):
    out = {}; owner = {}
    for f in files:
        if not os.path.exists(f):
            continue
        for k, v in pickle.load(open(f, 'rb')).items():
            if k[0] < DREF and float(np.max(np.abs(v['err']))/np.max(np.abs(v['E']))) > MAXREL:
                if (f, k) not in DROPPED:
                    DROPPED.add((f, k))
                    print(f"# cl_chain: dropped failed node {os.path.basename(f)} Delta={k[0]:g} nu={k[1]:.4f} (max rel err "
                          f"{float(np.max(np.abs(v['err']))/np.max(np.abs(v['E']))):.2f} > {MAXREL})", flush=True)
                continue
            if k[0] == DREF:
                assert k not in owner, ('cl_chain: reference key in two files', k, owner.get(k), f)
                owner[k] = f
            assert k not in out, ('cl_chain: point key in two files', k, f)
            out[k] = v
    return out


def hc_pts(kind, J, D):
    return merged([os.path.join(C, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), os.path.join(C, 'cl_pts_%s_J%d_D%g.pkl' % (kind, J, D))])


def oo_pts(kind, J, D):
    return merged([os.path.join(C, 'oo_pts_%s_J%d_D%g.pkl' % (kind, J, D)), os.path.join(C, 'cl_oopts_%s_J%d_D%g.pkl' % (kind, J, D))])


Z3NUMAX = float(os.environ.get('CL_Z3NUMAX', '60'))       # z4 direct small-nu nodes: only nu <= this (edge nodes go via HC_INJ)
Z3OO = {int(a.split(':')[0]): [float(d) for d in a.split(':')[1].split(',')] for a in os.environ.get('CL_Z3OO', '').split(';') if a}


def cap(P):
    return {k: v for k, v in P.items() if k[1] <= Z3NUMAX}


def oo_even(kind, J, D):
    """parity-even 2x2 block of the 3x3 oo points (as oo_chain._oo_even)."""
    out = {}
    for k, v in oo_pts(kind, J, D).items():
        w = dict(v)
        w.update(E=v['E'][:2, :2], err=v['err'][:2, :2], reps=v['reps'][:, :2, :2], method=v.get('method', '') + ' (oo even block)')
        out[k] = w
    return out


def z3pts(kind, J, D):
    """hc_chain z4 mode, HC_Z3 rows: the EE-block small-nu nodes (nu <= CL_Z3NUMAX) from hc_pts+cl_pts, or, for (J, D) in CL_Z3OO,
    from the even 2x2 block of oo_pts+cl_oopts (3x3 points, as oo_chain OO_EESRC=oo)."""
    return cap(oo_even(kind, J, D) if D in Z3OO.get(J, []) else hc_pts(kind, J, D))


# hc_chain injects HC_INJ nodes of kind C3H only through q3._q2_even, which the w2 loader uses for J >= 11; for even J < 11 the C3H
# rows come from W2._q3_load (l9 caches, phase-stripped) -- inject there as well (CL_C3HLOW = 1, default), with the same rule.
_C3HLOW = """_cl_q3l = W2._q3_load


def _cl_q3_load(kind, J, sfx):
    R = _cl_q3l(kind, J, sfx)
    if kind == 'C3H' and J < 11 and J % 2 == 0 and sfx == '_k2':
        R = _inject_dict(R, 'C3H', J, strip=True)
    return R


W2._q3_load = _cl_q3_load
"""


def _sub(src, old, new):
    assert src.count(old) == 1, ('cl_chain: pattern count', old, src.count(old))
    return src.replace(old, new)


if __name__ == '__main__':
    which = sys.argv[1]
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    import hc_nodes as N
    N.pts = hc_pts
    if which == 'hc':
        path = os.path.join(HERE, 'hc_chain.py')
        src = open(path).read()
        src = _sub(src, "Z3.EXPK = os.path.join(C, 'hc_exact.pkl')", "Z3.EXPK = os.path.join(C, 'cl_exact.pkl')")
        src = _sub(src, "            for k, v in N.pts(kind, J, D).items():", "            for k, v in _CL_Z3PTS(kind, J, D).items():")
        if os.environ.get('CL_C3HLOW', '1') == '1':
            src = _sub(src, 'pr(f"# hc_chain {mode}', _C3HLOW + 'pr(f"# hc_chain {mode}')
        src = _sub(src, "'hc_chain_%s_%s.txt'", "'cl_chain_%s_%s.txt'")
        src = _sub(src, "'hc_chain_%s_%s.pkl'", "'cl_chain_%s_%s.pkl'")
    elif which == 'oo':
        path = os.path.join(HERE, 'oo_chain.py')
        src = open(path).read()
        src = _sub(src, "Z3.EXPK = os.path.join(C, 'hc_exact.pkl')", "Z3.EXPK = os.path.join(C, 'cl_exact.pkl')")
        src = _sub(src, "    hc = os.path.join(C, 'hc_exact.pkl')", "    hc = os.path.join(C, 'cl_exact.pkl')")
        src = _sub(src, "OOEXPK = os.path.join(C, 'oo_exact.pkl')", "OOEXPK = os.path.join(C, 'cl_ooexact.pkl')")
        src = _sub(src, "    f = os.path.join(C, 'oo_pts_%s_J%d_D%g.pkl' % (kind, J, D))\n"
                        "    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}",
                   "    return _CL_OOPTS(kind, J, D)")
        src = _sub(src, "(N.pts(kind, J, D) if EESRC == 'hc' else _oo_even(kind, J, D)).items()",
                   "_CL_CAP(N.pts(kind, J, D) if EESRC == 'hc' else _oo_even(kind, J, D)).items()")
        src = _sub(src, "'oo_chain_%s.txt'", "'cl_oochain_%s.txt'")
        src = _sub(src, "'oo_chain_%s.pkl'", "'cl_oochain_%s.pkl'")
    else:
        raise SystemExit('cl_chain.py hc|oo ...')
    g = dict(__name__='__main__', __file__=path, __builtins__=__builtins__, _CL_OOPTS=oo_pts, _CL_Z3PTS=z3pts, _CL_CAP=cap)
    exec(compile(src, path, 'exec'), g)
