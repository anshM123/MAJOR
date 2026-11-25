"""
b4_k2cmp (BATCH4_FOLD.md §4): the direct check of "K2 = crystal integral" at Delta_g = 2000 (GAP_NEARGAP §8 caveat 1).  From this
note's outputs only (no computation):
  * state level: lambda/N_s on the own-mass row 2000 + J, K2 row (gapnear's 2e4 rows, ng_bind central value; hs_odd 'K2/<edge>') against
    the direct row (own-row central value; hs_odd 'dir/<edge>'), with the K2 lambda-level step 3|lam_K2 - lam_K1| the K2-only budget
    charged for truncation;
  * node level (odd J, hs_odd node tables): Gamma_dir - Gamma_K2 on the common nodes x <= 0.775, in units of the jackknife sigma and of
    the per-node K2 step |K2 - K1|.
Usage: python3 b4_k2cmp.py   -> cache/b4_k2cmp.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401
import numpy as np
C = os.path.join(HERE, 'cache')
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def ld(p):
    p = os.path.join(C, p)
    return json.load(open(p)) if os.path.exists(p) else None


if __name__ == '__main__':
    pr("# state level, own-row route: lam_dir - lam_K2 (central values) vs the K2 lambda-level truncation charge 3|K2 - K1|")
    for J in (12, 16, 21, 29):
        st = ld('state/b4_state_J%d.json' % J)
        for s in (st or {}).get('states', []):
            k2 = s['k2']
            if not k2 or 2000.0 + J not in k2['rows']:
                continue
            d = s['direct']['lam_central'] - k2['lam_central']
            pr(f"  J={J} m/Lam={s['mfl']:.4f}: K2 {k2['lam_central']:+.4f} direct {s['direct']['lam_central']:+.4f} diff {d:+.4f} "
               f"(jk {s['jk_sigma']:.4f}) | K2 charge 3|K2-K1| = {k2['E_trunc_lam']:.3f} -> |diff| / charge = {abs(d)/k2['E_trunc_lam']:.3f}")
    pr("# state level, hs_odd route: dir/<edge> - K2/<edge> (same edge data; bulk direct vs K2 shape fit)")
    for J in (21, 27, 29):
        h = ld('odd_main/hs_odd_J%d.json' % J)
        for k, o in sorted((h or {}).get('states', {}).items(), key=lambda t: float(t[0])):
            e = o['main_key'].split('/')[1]
            if 'K2/' + e in o:
                pr(f"  J={J} m/Lam={float(k):.4f}: K2/{e} {o['K2/' + e]:+.4f} dir/{e} {o['dir/' + e]:+.4f} diff {o['dir/' + e] - o['K2/' + e]:+.4f} "
                   f"(jk {o['jk']:.4f}); K2 bulk K1 {o['K1']:+.4f} -> K2 charge 3|K2-K1| {3*abs(o['K2/' + e] - o['K1']):.3f}")
    pr("# node level (odd J): Gamma_dir - Gamma_K2 on common nodes 0.05 <= x <= 0.775 (c2)")
    for J in (21, 27, 29):
        h = ld('odd_main/hs_odd_J%d.json' % J)
        if not h:
            continue
        rows = [n for n in h['nodes'] if 0.05 - 1e-9 <= n['x'] <= 0.775 + 1e-9 and n['C2']['K2'] is not None]
        for kd in ('C2', 'C3H'):
            d = np.array([n[kd]['dir'] - n[kd]['K2'] for n in rows]); e = np.array([n[kd]['err'] for n in rows])
            pr(f"  J={J} {kd}: {len(rows)} nodes; rms diff {np.sqrt(np.mean(d**2)):.0f}, max |diff| {np.abs(d).max():.0f} c2; "
               f"chi2/n (jk) {np.mean((d/e)**2):.2f}; rms jk sigma {np.sqrt(np.mean(e**2)):.0f}")
    open(os.path.join(C, 'b4_k2cmp.txt'), 'w').write("\n".join(lines) + "\n")
