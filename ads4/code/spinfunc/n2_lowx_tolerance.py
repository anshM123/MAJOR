"""
n2_lowx_tolerance (doc 24.8): range of t such that n = 1 at Delta_g = 100 when the graviton-specific correction on x < 0.2 (the region the
direct integration cannot reach, nu <~ 20) is replaced by t x (its K2-chain value), everything on x >= 0.2 direct (n1_analysis rows).
Also t for "correction on x < 0.2 removed" (t = 0) and the per-J minima.  Usage: python3 n2_lowx_tolerance.py -> cache/n2_lowx_tolerance.txt
"""
import os, sys, io, contextlib
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import n1_analysis as N
import n2_neargap_split as NG
tmp = os.path.join(HERE, 'cache', 'n1_tmp'); os.makedirs(os.path.join(tmp, 'cache'), exist_ok=True)
NG.HERE = tmp
NG.A.grid = lambda J, K, kinds=('C2', 'C3H'), **kw: N.grid_direct(J, K, kinds, **kw)
lines = []
def run(t):
    NG.TLOW[0] = t
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        NG.main('C23H_K2', [100])
    out = f.getvalue()
    l = [s for s in out.split('\n') if 'Delta_g = 100' in s][0]
    n = float(l.split('n = ')[1].split(' ')[0]) if 'n = None' not in l else None
    lam = float(l.split('min lam/N_s = ')[1].split(' ')[0])
    return n, lam, out
for t in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
    n, lam, out = run(t)
    jm = [s.strip()[:60] for s in out.split('\n') if s.strip().startswith('J=') and 'min lam' in s][:6]
    s = f"t = {t:+.2f}: n = {n}, min lam/N_s = {lam:+.4f};  " + "; ".join(jm[2:])
    print(s, flush=True); lines.append(s)
def ok(t):
    return run(t)[0] == 1.0
for sgn in (1, -1):
    lo, hi = 1.0, 1.0 + sgn
    while ok(hi) and abs(hi) < 50:
        lo, hi = hi, hi + sgn*abs(hi - 1.0 + sgn)
    for _ in range(12):
        mid = 0.5*(lo + hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    s = f"edge ({'+' if sgn > 0 else '-'}): n = 1 holds up to t = {lo:+.2f}"
    print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 'n2_lowx_tolerance.txt'), 'w').write("\n".join(lines) + "\n")
