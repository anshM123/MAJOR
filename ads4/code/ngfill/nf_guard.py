"""
nf_guard (NEARGAP_FILL.md): process-local write guard (a copy of neargap2/ng_guard.py with the root moved).  Every file opened for
writing (open(..., 'w'/'a'/'x'/'+')), every os.replace / os.rename target and every os.remove inside a ngfill process must lie
under ads4/code/ngfill/cache/ (or /dev/null); anything else raises PermissionError.
It also registers itself as 'ng_guard' in sys.modules, so that importing neargap2 modules (which `import ng_guard`) does NOT install
the neargap2 guard (whose root is neargap2/cache) on top of this one: neargap2 and spinfunc are read-only from here.
Import it FIRST in every ngfill script.
"""
import builtins, os, io, sys
HSC = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))
_open = builtins.open
_replace, _rename, _remove = os.replace, os.rename, os.remove


def _ok(path):
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return True                     # file descriptors
    return p == '/dev/null' or p.startswith(HSC + os.sep)


def _guarded_open(file, mode='r', *a, **k):
    if any(c in mode for c in 'wax+') and not _ok(file):
        raise PermissionError('nf_guard: refusing to write %s' % file)
    return _open(file, mode, *a, **k)


def _gr(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('nf_guard: refusing to move onto %s' % dst)
    return _replace(src, dst, *a, **k)


def _grn(src, dst, *a, **k):
    if not _ok(dst):
        raise PermissionError('nf_guard: refusing to rename onto %s' % dst)
    return _rename(src, dst, *a, **k)


def _grm(p, *a, **k):
    if not _ok(p):
        raise PermissionError('nf_guard: refusing to remove %s' % p)
    return _remove(p, *a, **k)


if not getattr(builtins.open, '_nf_guard', False):
    _guarded_open._nf_guard = True
    builtins.open = _guarded_open
    io.open = _guarded_open
    os.replace, os.rename, os.remove = _gr, _grn, _grm
    os.unlink = _grm
sys.modules.setdefault('ng_guard', sys.modules[__name__])
for _alias in ('hs_guard', 'sn_guard', 'ef_guard', 'gn_guard'):
    sys.modules.setdefault(_alias, sys.modules[__name__])
