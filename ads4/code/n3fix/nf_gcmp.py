"""n3fix/nf_gcmp.py: Gamma (c2, chain jackknife sigma) at J = 26 nodes x = 0.85 ... 0.93 from the nf_gamma variants: local n1 (M = 11),
patched n3 (M = 10, seed 4242), and at x = 0.92 n1 M = 13 / n3 M = 12 seed 7000; at x = 0.9198, 0.9276 batch-3 external n1 vs patched
n3 seed 6000.  Usage: python3 nf_gcmp.py  -> stdout"""
import nf_guard  # noqa: F401
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
V = {v: json.load(open(os.path.join(HERE, 'cache', 'gamma', v, 'gamma_J26.json'))) for v in ('n1val', 'n3val', 'p1', 'p3', 'ovl')
     if os.path.exists(os.path.join(HERE, 'cache', 'gamma', v, 'gamma_J26.json'))}


def g(v, key, x):
    N = V[v]['nodes'][key]
    i = min(range(len(N['x'])), key=lambda i: abs(N['x'][i] - x))
    return (N['G'][i], N['jk'][i]) if abs(N['x'][i] - x) < 1e-5 else (None, None)


for key in V['n1val']['nodes']:
    print(key)
    for x in (0.85, 0.87, 0.90, 0.92, 0.93):
        row = f"  x={x:.2f}: n1 M11 {g('n1val', key, x)[0]:+8.0f}+-{g('n1val', key, x)[1]:4.0f}  n3 M10 {g('n3val', key, x)[0]:+8.0f}+-{g('n3val', key, x)[1]:4.0f}"
        if x == 0.92:
            for v, t in (('p1', 'n1 M13'), ('p3', 'n3 M12 s7000')):
                if v in V:
                    a, e = g(v, key, x); row += f"  {t} {a:+8.0f}+-{e:4.0f}"
        print(row)
    if 'ovl' in V:
        for x in (0.919787, 0.927558):
            a, e = g('n1val', key, x); b, f = g('ovl', key, x)
            print(f"  x={x:.4f}: ext n1 M11 {a:+8.0f}+-{e:4.0f}  n3 M10 s6000 {b:+8.0f}+-{f:4.0f}")
