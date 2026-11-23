"""
dc_guard (DIRCV.md): process-local write guard.  Import FIRST.  Any open() in a writing mode, os.replace/rename destination,
os.remove/unlink, os.mkdir/makedirs whose real path is not inside ads4/code/dircv/cache raises PermissionError.  Also disables
bytecode writing.  /dev/null is allowed.
"""
import os, sys, builtins
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.realpath(os.path.join(HERE, 'cache'))
SPIN = os.path.realpath(os.path.join(HERE, '..', 'spinfunc'))
N2 = os.path.realpath(os.path.join(HERE, '..', 'n2sub'))
ADS4 = os.path.realpath(os.path.join(HERE, '..', '..'))


def _ok(path):
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return True
    return p == '/dev/null' or p == CACHE or p.startswith(CACHE + os.sep)


def _check(path):
    if not _ok(path):
        raise PermissionError('dc_guard: write outside dircv/cache refused: %s' % path)


_open = builtins.open


def _gopen(file, mode='r', *a, **k):
    if isinstance(file, (str, bytes, os.PathLike)) and any(c in mode for c in 'wax+'):
        _check(file)
    return _open(file, mode, *a, **k)


builtins.open = _gopen
import io
io.open = _gopen
_rep, _ren, _rm, _ul, _mk, _mks = os.replace, os.rename, os.remove, os.unlink, os.mkdir, os.makedirs


def _grep(src, dst, *a, **k):
    _check(dst); _check(src); return _rep(src, dst, *a, **k)


def _gren(src, dst, *a, **k):
    _check(dst); _check(src); return _ren(src, dst, *a, **k)


def _grm(path, *a, **k):
    _check(path); return _rm(path, *a, **k)


def _gmk(path, *a, **k):
    _check(path); return _mk(path, *a, **k)


def _gmks(path, *a, **k):
    if k.get('exist_ok') and os.path.isdir(path):      # no-op on an existing directory
        return None
    _check(path); return _mks(path, *a, **k)


os.replace, os.rename, os.remove, os.unlink, os.mkdir, os.makedirs = _grep, _gren, _grm, _grm, _gmk, _gmks
os.environ['MPLCONFIGDIR'] = CACHE
for p in (SPIN, N2):
    if p not in sys.path:
        sys.path.insert(0, p)
sys.path.insert(0, HERE)
