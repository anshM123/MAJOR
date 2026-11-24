"""
gu_guard (GAP_NEARGAP.md): process-local write guard.  Every file opened for writing, every os.replace / os.rename target and every
os.remove inside a gapuni process must lie under ads4/code/gapuni/cache/ (or /dev/null); anything else raises.  It registers itself
as ng_guard / hs_guard / sn_guard / ef_guard, so importing neargap2 / highspin / smallnu / extfold modules installs no second guard
(theirs would refuse writes to gapuni/cache).  Copy of neargap2/ng_guard.py with the root moved.
"""
import builtins, os, io, sys
GUC = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
_open = builtins.open
_replace, _rename, _remove = os.replace, os.rename, os.remove


def _ok(path):
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return True
    return p == '/dev/null' or p.startswith(GUC + os.sep)


def _guarded_open(file, mode='r', *a, **k):
    if any(c in mode for c in 'wax+') and not _ok(file):
        raise PermissionError('gu_guard: refusing to write %s' % file)
    return _open(file, mode, *a, **k)


def _gr(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('gu_guard: refusing to move onto %s' % dst)
    return _replace(src, dst, *a, **k)


def _grn(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('gu_guard: refusing to rename onto %s' % dst)
    return _rename(src, dst, *a, **k)


def _grm(p, *a, **k):
    if not _ok(p):
        raise PermissionError('gu_guard: refusing to remove %s' % p)
    return _remove(p, *a, **k)


builtins.open = _guarded_open
io.open = _guarded_open
os.replace, os.rename, os.remove = _gr, _grn, _grm
for _n in ('ng_guard', 'hs_guard', 'sn_guard', 'ef_guard'):
    sys.modules.setdefault(_n, sys.modules[__name__])
