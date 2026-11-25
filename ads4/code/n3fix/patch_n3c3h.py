"""
patch_n3c3h (N3C3H_FIX.md): in-process fix of the n3 (hybrid wall-free contour) evaluation of the graviton k = 3 kind C3H, even J
(l9_points.odd_structure via n1_pw.P9) and odd J (p3_odd_action.odd_structure_H via p8_odd_pw.P8).  No project file is edited.

THE BUG.  The parity-odd photon structure o_H (doc 21.6) is the 2 <-> 4 symmetric combination
    o_H = e(Z2,X2,X4,X0,Z0) H40 n2 + H20 e(Z4,X4,X2,X0,Z0) n4,   n2 = (X24 X20 X40^3)^(-1/2),  n4 = (X24 X40 X20^3)^(-1/2),
coded as TWO independent fractional powers.  On the physical sheet n2 X40 = n4 X20 = (X24 X20 X40)^(-1/2), i.e. the relative sign
r = n2 X40 / (n4 X20) = +1.  At the symmetric crystal both power arguments are EXACTLY on the negative real axis (X20 = X40 < 0, X24 < 0).
  * n1 (and the TPS/Laplace engine) take principal logs at the (symmetric) saddle, where the two arguments are equal: r = +1.
  * n3_flow.set_branch continues every log from W_REF, a generic point 1e-4 off the crystal: W_REF's asymmetric displacement can put
    the two arguments on OPPOSITE sides of the cut, and the continuation then carries r = -1 to the whole contour.  Whether it does
    depends on (Delta, nu): at (26, 1026, x = 0.93) the state has r = -1, the Delta_ref = 2e4 reference r = +1.  With r = -1 the
    integrand is the 2 <-> 4 ANTI-symmetric combination, whose integral nearly vanishes: the "state x200-1000 too small, 30-600 %
    errors, references fine" of EDGE_R §2.4.  EE / C2 contain no such pair of powers and are unaffected.
THE FIX (two parts, both exact identities on the physical sheet):
  (a) n4 := n2 X40 / X20 in both point-wise copies (o_H then has r = +1 by construction, whatever branch n2 is on);
  (b) graviton evaluators of n3 take their saddle log branch = principal logs at the full saddle wc -- n1's convention (n1_direct.
      evaluate), so that the overall sign of n2 (and every other graviton log) is the one of the existing n1 / K2 records.  The scalar
      evaluator (contour, flow) keeps set_branch's W_REF continuation unchanged.
Usage:  import patch_n3c3h as PX; PX.install()   AFTER hc_points.setup(odd) (it needs PW.P9 and, for odd J, p8's P8 module).
        PX.enable_or(False) / PX.enable_branch(False) switch the two parts off (diagnostics).
"""
import os, sys
import numpy as np

SF = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'spinfunc'))
if not os.path.isdir(SF):                                    # ads4/ext/ copy: spinfunc is ../code/spinfunc
    SF = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code', 'spinfunc'))
OLD = "    n4 = T.powr(T.mul(T.mul(X24, X40), T.mul(X20, T.mul(X20, X20))), -0.5)\n"
NEW = "    n4 = T.mul(n2, T.div(X40, X20))          # n3fix: (X24 X20 X40)^-1/2 / X20 on n2's branch (relative sign r = +1)\n"
_ORIG = {}
STATE = dict(or_=False, branch=False)


def _func_src(path, name):
    src = open(path).read()
    i = src.index('def %s(' % name)
    j = src.index('\n\n\n', i)
    f = src[i:j] + '\n'
    assert f.count(OLD) == 1, (path, name)
    return f.replace(OLD, NEW)


def _modules():
    import n1_pw as PW
    mods = [(PW.P9, os.path.join(SF, 'l9_points.py'), 'odd_structure')]
    if 'p8_odd_pw' in sys.modules:
        mods.append((sys.modules['p8_odd_pw'].P8, os.path.join(SF, 'p3_odd_action.py'), 'odd_structure_H'))
    return mods


def enable_or(on=True):
    """part (a): swap the module-level o_H function of the point-wise copies (build functions look it up by name at call time)."""
    for M, path, name in _modules():
        key = (id(M), name)
        if key not in _ORIG:
            _ORIG[key] = M.__dict__[name]
            ns = dict(M.__dict__)
            exec(compile(_func_src(path, name), 'n3fix_' + name, 'exec'), M.__dict__)
            _ORIG[key + ('new',)] = M.__dict__[name]
            assert ns.keys() <= M.__dict__.keys()
        M.__dict__[name] = _ORIG[key + ('new',)] if on else _ORIG[key]
    STATE['or_'] = on


def enable_branch(on=True):
    """part (b): wrap n3_flow.set_branch so that graviton evaluators (kind != 'S') get principal logs at their saddle wc."""
    import n3_flow as F
    if '_set_branch' not in _ORIG:
        _ORIG['_set_branch'] = F.set_branch

        def set_branch(ev, w_from=None, nstep=200, fns=('phi', 'pint', 'grav')):
            if STATE['branch'] and ev.kind != 'S' and w_from is None:
                return _ORIG['_set_branch'](ev, w_from=np.array(ev.wc, complex), nstep=1, fns=('grav',))
            return _ORIG['_set_branch'](ev, w_from, nstep, fns)
        F.set_branch = set_branch
    STATE['branch'] = on


def install():
    enable_or(True)
    enable_branch(True)
    return dict(STATE)
