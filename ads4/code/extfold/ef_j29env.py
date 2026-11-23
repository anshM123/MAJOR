"""
ef_j29env (EXT_FOLD.md §4.3): the J = 29 sign-free x3 envelope (smallnu/sn_window.envelope / cost_x, read-only, NEARGAP2's gauge
factor; 1e4 below nu = 10) with the pooled n3 nodes, decomposed by node, and what-if scenarios:
  S0  SMALLNU (seed 4242 only)                       S1  pooled (8 super-replica jackknife; = ef_window 'all')
  S2  pooled, sigma from the 31-stream sd/sqrt(31)   S3  nu = 10-17.5 at |Gamma| only (sigma -> 0: infinite n3 statistics there)
  S4  S3 and the nu = 20 n1 node also at |Gamma|      S5  S1 without the nu = 20 node (n3 nodes bridge 17.5 -> 25)
and the uniform bound on [10, 25) that the twist edge tolerates (everything else as in S1).
Usage: python3 ef_j29env.py -> cache/ef_j29env.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C
import numpy as np
import sn_window as SW          # read-only import
J = 29
NEW = json.load(open(os.path.join(C.EFC, 'sn_odd_J29.json')))
OLD = json.load(open(os.path.join(C.EFC, 'old_repro', 'sn_odd_J29.json')))
STR = json.load(open(os.path.join(C.EFC, 'ef_streams.json')))
NU3 = [10.0, 12.5, 15.0, 17.5]


def nodes_of(R, numin=9.99):
    return [r for r in R['nodes'] if r['src'] != 'K2' and r['nu'] >= numin - 1e-6]


def bounds(nodes, sig=None, zero=(), drop=()):
    gx, gb, lab = [], [], []
    for r in nodes:
        nu = round(r['nu'], 1)
        if any(abs(nu - d) < 0.05 for d in drop):
            continue
        b = 0.0
        for k in ('C2', 'C3H'):
            v, e = r[k]
            if sig is not None and any(abs(nu - n) < 0.05 for n in NU3):
                e = sig[(round(nu, 1), k)]
            if any(abs(nu - z) < 0.05 for z in zero):
                e = 0.0
            b = max(b, abs(v) + 3*e)
        gx.append(r['x']); gb.append(b); lab.append(nu)
    return np.array(gx), np.array(gb), lab


if __name__ == '__main__':
    lines, res = [], {}

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    nn = nodes_of(NEW); no = nodes_of(OLD)
    sig31 = {(round(nu, 1), k): STR['nu'][str(nu)][k]['sigma'] for nu in NU3 for k in ('C2', 'C3H')}
    sc = {'S0 SMALLNU (4242 only)': bounds(no), 'S1 pooled (jackknife)': bounds(nn), 'S2 pooled (31-stream sd)': bounds(nn, sig31),
          'S3 nu=10-17.5 at |Gamma|': bounds(nn, zero=NU3), 'S4 S3 + nu=20 at |Gamma|': bounds(nn, zero=NU3 + [20.0]),
          'S5 S1 without the nu=20 n1 node': bounds(nn, drop=[20.0])}
    for key, s in NEW['states'].items():
        mfl = float(key); Eb = s['E_base']; lamc = s['lamc']
        pr(f"# m/Lam = {mfl:.7f}: lam_c {lamc:+.4f}, budget without window {Eb:.3f}; window budget allowed for L = 1: {lamc - Eb:.3f}")
        res[key] = {}
        for name, (gx, gb, lab) in sc.items():
            Ek = SW.cost_x(J, mfl, SW.envelope(gx, gb, 3, 1e4))
            L = lamc/(Eb + Ek)
            res[key][name] = dict(E=Ek, L=L, bounds=dict(zip([str(l) for l in lab], gb.tolist())))
            pr(f"   {name:34s}: env x3 {Ek:.3f} -> L {L:.2f} | node bounds " + ", ".join(f"{l:g}:{b:.0f}" for l, b in zip(lab, gb) if l < 50))
        # per-segment cost of S1
        gx, gb, lab = sc['S1 pooled (jackknife)']
        o = np.argsort(gx); gx, gb = gx[o], gb[o]; lab = [lab[i] for i in o]
        m = mfl*SW.E.lam_of(1000.0)
        seg = []
        for i in range(len(gx) - 1):
            if gx[i]*m >= 50:
                break
            g = 3*max(gb[i], gb[i + 1])
            seg.append((lab[i], lab[i + 1], g, SW.cost(J, mfl, g, gx[i]*m, min(gx[i + 1]*m, 50.0))))
        pr("   S1 per segment (nu_a-nu_b: 3 max(g) -> cost): " + ", ".join(f"{a:g}-{b:g}: {g:.0f} -> {c:.3f}" for a, b, g, c in seg) +
           f"; below 10 (1e4): {SW.cost(J, mfl, 1e4, 0.0, gx[0]*m):.3f}")
        # tolerated uniform bound on [10, 25) with S1 elsewhere
        def cost_with(G):
            f0 = SW.envelope(gx, gb, 3, 1e4)
            return SW.cost_x(J, mfl, lambda x: np.where((x*m >= 10 - 1e-9) & (x*m < 25), G, f0(x)))
        lo, hi = 0.0, 1e4
        for _ in range(40):
            mid = 0.5*(lo + hi)
            if lamc/(Eb + cost_with(mid)) >= 1.0:
                lo = mid
            else:
                hi = mid
        res[key]['G_tol_10_25'] = lo
        pr(f"   tolerated uniform envelope on [10, 25) for L = 1 (x3 envelope elsewhere): {lo:.0f} c2 (i.e. node bounds |Gamma| + 3 sigma <= {lo/3:.0f})")
    open(os.path.join(C.EFC, 'ef_j29env.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.EFC, 'ef_j29env.json'), 'w'), indent=1, default=float)
