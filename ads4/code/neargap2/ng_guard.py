"""
ng_guard (NEARGAP2.md): process-local write guard.  Every file opened for writing (open(..., 'w'/'a'/'x'/'+')), every os.replace /
os.rename target and every os.remove inside a neargap2 process must lie under ads4/code/neargap2/cache/ (or /dev/null); anything else
raises.  Project modules are imported read-only; their cache WRITES are redirected by the ng_* wrappers, and this guard makes sure that
no redirection was missed (it never silently writes into spinfunc/cache).
"""
import builtins, os, io
NGC = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
_open = builtins.open
_replace, _rename, _remove = os.replace, os.rename, os.remove


def _ok(path):
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return True                     # file descriptors
    return p == '/dev/null' or p.startswith(NGC + os.sep)


def _guarded_open(file, mode='r', *a, **k):
    if any(c in mode for c in 'wax+') and not _ok(file):
        raise PermissionError('ng_guard: refusing to write %s' % file)
    return _open(file, mode, *a, **k)


def _gr(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('ng_guard: refusing to move onto %s' % dst)
    return _replace(src, dst, *a, **k)


def _grn(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('ng_guard: refusing to rename onto %s' % dst)
    return _rename(src, dst, *a, **k)


def _grm(p, *a, **k):
    if not _ok(p):
        raise PermissionError('ng_guard: refusing to remove %s' % p)
    return _remove(p, *a, **k)


builtins.open = _guarded_open
io.open = _guarded_open
os.replace, os.rename, os.remove = _gr, _grn, _grm
