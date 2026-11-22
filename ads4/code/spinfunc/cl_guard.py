"""
cl_guard (RESULT_D_CLEANUP.md): write guard for the cl_* scripts.  Any open() in a writing mode, os.replace / os.rename destination
or os.remove inside ads4/code/spinfunc/cache whose basename does not start with 'cl_' raises -- so no existing cache can be modified by
the unchanged project code that the cl_ wrappers execute in-process.  Paths outside the cache directory are not restricted.
"""
import os, builtins
CACHE = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache'))


def _check(path):
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return
    if os.path.dirname(p) == CACHE and not os.path.basename(p).startswith('cl_'):
        raise PermissionError('cl_guard: write to non-cl cache file refused: %s' % p)


_open = builtins.open


def _gopen(file, mode='r', *a, **k):
    if isinstance(file, (str, bytes, os.PathLike)) and any(c in mode for c in 'wax+'):
        _check(file)
    return _open(file, mode, *a, **k)


builtins.open = _gopen
_rep, _ren, _rm = os.replace, os.rename, os.remove


def _grep(src, dst, *a, **k):
    _check(dst); return _rep(src, dst, *a, **k)


def _gren(src, dst, *a, **k):
    _check(dst); _check(src); return _ren(src, dst, *a, **k)


def _grm(path, *a, **k):
    _check(path); return _rm(path, *a, **k)


os.replace, os.rename, os.remove = _grep, _gren, _grm
