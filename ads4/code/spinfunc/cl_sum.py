"""
cl_sum (RESULT_D_CLEANUP.md): summaries of the cl_chain outputs.
  python3 cl_sum.py hc <mode> <Dg> <prefix> <main> [variants...] [--jk=<jkprefix>] [--ref=<variant>] [--J=..] [--mmax=1.3]
      per J: min lambda/N_s over near-gap states (m/Lam <= mmax) and over all states of the main variant; jackknife error (variants
      <jkprefix>0..7) at that state and max over near-gap states; every other variant at the same state and its own near-gap min;
      number of negative eigen-directions (states x helicity eigenvalues < 0) per variant; with --ref: per-state shift main - ref
      (max |shift| near gap, shift at the main minimum).
  python3 cl_sum.py oo <prefix> <main> [variants...] [--jk=jk]      (cl_oochain outputs: OO block and EE block (run0))
"""
import os, sys, pickle
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
MM = float(opt.get('mmax', '1.3'))
JS = [int(v) for v in opt.get('J', '').split(',') if v]


def nneg(r):
    if 'base' not in r:
        return int(r['lam'] < 0)
    Fm = np.diag(r['base']) + r['dF']; n = 1/np.sqrt(r['Ns'])
    return int(np.sum(np.linalg.eigvalsh(n[:, None]*Fm*n[None, :]) < 0))


def jkerr(V):
    V = np.asarray(V, float); n = len(V)
    return float(np.sqrt((n - 1)/n*((V - V.mean())**2).sum())) if n > 1 else np.nan


if args[0] == 'hc':
    mode, Dg, pre, main = args[1], int(args[2]), args[3], args[4]
    others = args[5:]

    def load(v):
        f = os.path.join(C, 'cl_chain_%s_%s_%s.pkl' % (mode, pre, v))
        return pickle.load(open(f, 'rb'))[Dg] if os.path.exists(f) else None

    def key(r):
        return (round(r['mfl'], 6), r['DF'])
    M = load(main)
    jk = [load('%s%d' % (opt['jk'], i)) for i in range(8)] if 'jk' in opt else []
    jk = [j for j in jk if j is not None]
    ref = load(opt['ref']) if 'ref' in opt else None
    for J in sorted(set(r['J'] for r in M)):
        if JS and J not in JS:
            continue
        R = [r for r in M if r['J'] == J]
        near = [r for r in R if r['mfl'] <= MM + 1e-9]
        if not near:
            continue
        w = min(near, key=lambda r: r['lam']); wa = min(R, key=lambda r: r['lam'])
        s = f"J={J:2d} [{pre}_{main}] near-gap min {w['lam']:+.4f} at {w['mfl']:.4f} ({w['DF']})"
        if jk:
            idx = {key(r): r for r in near}
            vals = [{key(r): r['lam'] for r in j if r['J'] == J} for j in jk]
            e0 = jkerr([v[key(w)] for v in vals])
            emax = max(jkerr([v[k] for v in vals]) for k in idx)
            s += f" +- {e0:.4f} (jk {len(jk)}; max over near-gap states {emax:.4f})"
        s += f"; all m {wa['lam']:+.4f} at {wa['mfl']:.3f}; negative directions {sum(nneg(r) for r in R)} of {sum(len(r.get('hs', [0])) for r in R)}"
        print(s)
        if ref is not None:
            rr = {key(r): r for r in ref if r['J'] == J}
            d = [(r['lam'] - rr[key(r)]['lam'], r['mfl']) for r in near if key(r) in rr]
            wr = min((r for r in ref if r['J'] == J and r['mfl'] <= MM + 1e-9), key=lambda r: r['lam'])
            if d:
                dm = max(d, key=lambda t: abs(t[0]))
                print(f"      vs {opt['ref']}: ref near-gap min {wr['lam']:+.4f} at {wr['mfl']:.4f}; shift at main-min state "
                      f"{w['lam'] - rr[key(w)]['lam']:+.4f}; max |shift| near gap {dm[0]:+.4f} at {dm[1]:.4f}; shift of the min "
                      f"{w['lam'] - wr['lam']:+.4f}")
        vs = []
        for v in others:
            o = load(v)
            if o is None:
                continue
            Rv = [r for r in o if r['J'] == J]
            nv = [r for r in Rv if r['mfl'] <= MM + 1e-9]
            if not nv:
                continue
            u = min(nv, key=lambda r: r['lam'])
            same = next((r['lam'] for r in nv if key(r) == key(w)), np.nan)
            vs.append(f"{v} {u['lam']:+.4f}[{u['mfl']:.3f}] (at main state {same:+.4f}; neg {sum(nneg(r) for r in Rv)})")
        if vs:
            print("      variants (near-gap min [m]): " + "; ".join(vs))
        print("      near-gap states: " + " ".join(f"{r['mfl']:.4f}:{r['lam']:+.4f}" for r in sorted(near, key=lambda r: (r['mfl'], r['lam']))
                                             if r is min((q for q in near if abs(q['mfl'] - r['mfl']) < 1e-9), key=lambda q: q['lam'])))
elif args[0] == 'oo':
    pre, main = args[1], args[2]
    others = args[3:]

    def load(v):
        f = os.path.join(C, 'cl_oochain_%s_%s.pkl' % (pre, v))
        return pickle.load(open(f, 'rb'))['states'] if os.path.exists(f) else None
    M = load(main)
    jk = [load('%s%d' % (opt.get('jk', 'jk'), i)) for i in range(8)]
    jk = [j for j in jk if j is not None]
    for (Dg, J), rows in sorted(M.items()):
        near = [t for t in rows if t['mfl'] <= MM + 1e-9]
        for blk, kk in (('OO', 'oo'), ('EE', 'lam_ee')):
            if all(np.isnan(t[kk]) for t in near):
                continue
            w = min(near, key=lambda t: t[kk]); wa = min(rows, key=lambda t: t[kk])
            s = f"J={J:2d} Dg={Dg} {blk} [{pre}_{main}] near-gap min {w[kk]:+.4f} at {w['mfl']:.4f} ({w['DF']})"
            jj = [j for j in jk if not np.isnan(j[(Dg, J)][0][kk])]
            if jj:
                val = lambda s_, m, DF: next(t[kk] for t in s_[(Dg, J)] if abs(t['mfl'] - m) < 1e-9 and t['DF'] == DF)
                e0 = jkerr([val(s_, w['mfl'], w['DF']) for s_ in jj])
                emax = max(jkerr([val(s_, t['mfl'], t['DF']) for s_ in jj]) for t in near)
                s += f" +- {e0:.4f} (jk {len(jj)}; max {emax:.4f})"
            s += f"; all m {wa[kk]:+.4f} at {wa['mfl']:.2f}; negative states {sum(t[kk] < 0 for t in rows)} of {len(rows)}"
            print(s)
            vs = []
            for v in others:
                o = load(v)
                if o is None or (Dg, J) not in o:
                    continue
                nv = [t for t in o[(Dg, J)] if t['mfl'] <= MM + 1e-9]
                if all(np.isnan(t[kk]) for t in nv):
                    continue
                u = min(nv, key=lambda t: t[kk])
                vs.append(f"{v} {u[kk]:+.4f}[{u['mfl']:.3f}]")
            if vs:
                print("      variants: " + "; ".join(vs))
