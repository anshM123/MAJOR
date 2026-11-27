"""md_hhat_total.py -- total stripped light action hhat = hhat^(2) (md_hhat_quad: nu-quadrature of the published weight) + hhat^(3)
(md_k3: exact Q_k), compared with the moment-route values of no_eft/cache/nl1_<Dg>_<K>.txt (read-only)."""
import os, sys, re
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
import sympy as sp
import md_k3 as K3

mp.mp.dps = 40
L = []
for Dg in (1000, 300):
    lam = (Dg - 2)*(Dg - 1)
    Q = K3.Qk(30, lam)
    ref = {}
    fn = os.path.join(HERE, '..', 'no_eft', 'cache', f'nl1_{Dg}_{Dg}.txt')
    if os.path.exists(fn):
        for line in open(fn):
            m = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+sign ([+-])\s+log10\|hhat\| =\s+([\d.]+)', line)
            if m:
                ref[(int(m.group(1)), int(m.group(2)), int(m.group(3)))] = (m.group(4), float(m.group(5)))
    L.append(f'# Delta_g = {Dg}: hhat = hhat2 (quadrature) + hhat3 (exact); no_eft moment-route value in brackets')
    for line in open(os.path.join(HERE, 'out', f'hhat_quad_{Dg}.txt')):
        m = re.match(r'\s*twist\s+(\d+) J\s+(\d+) b (\d): hhat\^\(2\) = (\S+)', line)
        if not m:
            continue
        tw, J, b, h2 = int(m.group(1)), int(m.group(2)), int(m.group(3)), mp.mpf(m.group(4))
        h3 = mp.mpf(sp.Float(K3.hhat3(b, tw + J, J, Q, lam), 30))
        tot = h2 + h3
        r = ref.get((b, tw, J))
        L.append(f'  twist {tw:4d} J {J:3d} b {b}: hhat2 = {mp.nstr(h2, 6):>14s}  hhat3 = {mp.nstr(h3, 6):>14s}  total = {mp.nstr(tot, 8):>16s}'
                 f'  log10|total| = {float(mp.log10(abs(tot))):8.3f}' + (f'   [no_eft: {r[0]} {r[1]:.3f}]' if r else ''))
print('\n'.join(L)); open(os.path.join(HERE, 'out', 'hhat_total.txt'), 'w').write('\n'.join(L) + '\n')
