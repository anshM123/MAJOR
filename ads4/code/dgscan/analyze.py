"""
dgscan/analyze.py (DG_SCAN.md §2-§3): margins of per-state positivity of Omega^EE vs Delta_g from
  (i)   cache/tol_<Dg>.json (tol.py): baseline lambda0, uniform-window thresholds C*, H_LJ window tolerance T_B, bulk tolerance;
  (ii)  the computed chain verdicts at 100 / 300 / 1000 (parsed from spinfunc/cache, read-only):
          y2_chain100_jhi.txt, y2_chain100_jhi_scale.txt (Delta_g = 100, EE nonlocal min-norm, J >= 12 transferred from Delta >= 336 rows),
          w2_chain_final_xc0.txt (300 and 1000, "convention A": K2 rows, spline below the lowest node, local de-shift);
  (iii) the adopted-protocol per-state values at 300 (RESULT_D_STATEMENT §5.3, quoted below with their masses).
Proxies:
  P1  fixed-(J, m/Lam) scaling of the correction: corr(Lam) = corr(300) (298.5/Lam)^p,  p = 2 (bulk, Gamma(x) held), 4 (small-nu window,
      Gamma(nu) held), p_emp from the conv-A t*_J at 300 vs 1000 (same code, own rows at both gaps);
  P2  sign-free H_LJ window: min_J T_B(c = 3) vs a Delta-independent window correction |Gamma_win| = 450 c2 (300, J = 24, 28 fold fits;
      band 200-1000).
Usage: python3 analyze.py -> stdout + cache/analyze.txt, cache/analyze.json
"""
import os, sys, re, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def lam_of(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


# ------------------------------------------------------------------------------------------------ (i) tolerance scan
def load_tol():
    T = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'tol_*.json'))):
        d = json.load(open(fn))
        T[d['Dg']] = d
    return T


def perJ(d, key, agg=min):
    out = {}
    for r in d['states']:
        v = r[key]
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            continue
        out.setdefault(r['J'], []).append(v)
    return {J: agg(v) for J, v in out.items()}


# ------------------------------------------------------------------------------------------------ (ii) chain outputs
def parse_chain(fn, block_pat):
    """-> {J: (min, uncorr, t*_upper)} from the first block whose header matches block_pat."""
    txt = open(os.path.join(SF, 'cache', fn)).read().splitlines()
    res, on = {}, False
    for ln in txt:
        if re.search(block_pat, ln):
            on = True; continue
        if on and re.match(r'\s*(\[|Delta_g =)', ln) and res:
            break
        if on:
            m = re.match(r'\s*J=\s*(\d+): min (?:lam/N_s )?([+-][\d.]+) \(uncorr\. ([+-][\d.]+)\).*t\*_J = \(([^,]+), ([^)]+)\)', ln)
            if m:
                tu = m.group(5).strip()
                res[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float('inf') if 'inf' in tu else float(tu))
    return res


# ------------------------------------------------------------------------------------------------ (iii) adopted 300 values
# RESULT_D_STATEMENT §5.3: (J, min lam/N_s, m/Lam of the minimum; None = gap/twist edge)
ADOPT300 = [(4, 0.2137, 1.013), (5, 0.2038, 1.0166), (6, 0.3075, 1.25), (7, 0.7033, 1.095), (8, 0.4268, 1.055), (9, 0.4506, 1.03),
            (10, 0.6033, 2.0), (11, 0.6687, 1.036), (12, 0.6706, 1.055), (13, 0.7746, 1.10), (14, 0.7487, 1.046), (15, 0.7853, 1.07),
            (16, 0.884, None), (18, 0.872, 4.0), (20, 0.905, 1.065), (22, 0.844, 1.071), (24, 0.820, 1.5), (28, 0.799, 1.5),
            (32, 0.783, 1.15)]


def lam0_at(Dg, J, mfl):
    import ee1000_common as E
    mg = E.mgap(Dg, J)
    mfl = mg if mfl is None else max(mfl, mg)
    hs = ('+-',) if J % 2 else ('++', '+-')
    v = []
    for h in hs:
        st = E.hyb(Dg, J, mfl, h, 256)
        v.append(st['base']/st['Ns'])
    return min(v), mfl


def main():
    T = load_tol()
    Dgs = sorted(T)
    pr("# dgscan/analyze.py")
    # ---------------------------------------------------------------- baseline
    pr("\n## A. baseline (leading-Regge + curvature hybrid states) [M]: min lambda0/N_s over J <= 0.345(Dg-1), m in [m_gap, 1.3 Lam], h")
    pr("| Dg | Lam | #states | min lambda0 (state) | J = 5 gap | min over J >= 6 (state) |")
    pr("|---|---|---|---|---|---|")
    base = {}
    for Dg in Dgs:
        S = T[Dg]['states']
        mn = min(S, key=lambda r: r['lam0'])
        j5 = [r['lam0'] for r in S if r['J'] == 5 and r['gap']]
        mn6 = min((r for r in S if r['J'] >= 6), key=lambda r: r['lam0'])
        base[Dg] = mn['lam0']
        pr(f"| {Dg} | {T[Dg]['Lam']:.1f} | {len(S)} | {mn['lam0']:+.4f} (J={mn['J']}, {mn['h']}, {mn['mfl']:.4f}) | {min(j5):+.4f} |"
           f" {mn6['lam0']:+.4f} (J={mn6['J']}, {mn6['h']}, {mn6['mfl']:.3f}) |")
    # ---------------------------------------------------------------- uniform-window thresholds at fixed J
    pr("\n## B. failure threshold of a uniform correction on nu < 25 (C*_25, c2 units; min |C*| over m <= 1.3 Lam and h), fixed J")
    Js = [5, 8, 12, 16, 20, 24, 28, 32, 40, 48, 64, 100]
    pr("| Dg | " + " | ".join(f"J={J}" for J in Js) + " |")
    pr("|---|" + "---|"*len(Js))
    C25 = {Dg: perJ(T[Dg], 'Cs25', agg=lambda v: max(v)) for Dg in Dgs}       # C* negative: max = smallest |C*|
    for Dg in Dgs:
        pr(f"| {Dg} | " + " | ".join((f"{C25[Dg][J]:.3g}" if J in C25[Dg] else "—") for J in Js) + " |")
    pr("local exponents d log|C*_25| / d log Lam at fixed J (expect -> 4 once 25/Lam << 1):")
    for J in (5, 16, 24, 32):
        s = []
        for a, b in zip(Dgs[:-1], Dgs[1:]):
            if J in C25[a] and J in C25[b]:
                s.append(f"{a}-{b}: {np.log(C25[b][J]/C25[a][J])/np.log(lam_of(b)/lam_of(a)):.2f}")
        pr(f"  J={J}: " + ", ".join(s))
    TA = {Dg: perJ(T[Dg], 'tol_all') for Dg in Dgs}
    pr("bulk (uniform over all nu, psi2 part, 0.9 lambda0) tolerance exponents (expect 2):")
    for J in (5, 16, 24, 32):
        s = []
        for a, b in zip(Dgs[:-1], Dgs[1:]):
            if J in TA[a] and J in TA[b]:
                s.append(f"{a}-{b}: {np.log(TA[b][J]/TA[a][J])/np.log(lam_of(b)/lam_of(a)):.2f}")
        pr(f"  J={J}: " + ", ".join(s))
    # ---------------------------------------------------------------- H_LJ window tolerance
    pr("\n## C. H_LJ sign-free window tolerance T_B(c=3) (C2, worse helicity, min over m <= 1.3 Lam; HLJ_ANALYTIC §2 without W_X)")
    pr("| Dg | J_K2 = sqrt(1.5 Lam) | J_cap | min T_B3 over J >= J_K2 (J) | min T_B3 over all near-gap J (J) | T_B3 at J/Lam = 0.1 / 0.2 / 0.3 | min T_B5 (J >= J_K2) |")
    pr("|---|---|---|---|---|---|---|")
    TB = {}
    for Dg in Dgs:
        Lm = lam_of(Dg)
        tb3 = perJ(T[Dg], 'TB3'); tb5 = perJ(T[Dg], 'TB5')
        TB[Dg] = tb3
        jk2 = int(np.ceil(np.sqrt(1.5*Lm))); jcap = int(0.345*(Dg - 1))
        hi = {J: v for J, v in tb3.items() if J >= jk2 and J >= 4}
        al = {J: v for J, v in tb3.items() if J >= 4}
        hi5 = {J: v for J, v in tb5.items() if J >= jk2 and J >= 4}
        jh = min(hi, key=hi.get) if hi else None; ja = min(al, key=al.get)
        fr = []
        for f in (0.1, 0.2, 0.3):
            Jt = f*Lm
            ks = sorted(tb3)
            Jn = min(ks, key=lambda J: abs(J - Jt))
            fr.append(f"{tb3[Jn]:.0f} (J={Jn})")
        pr(f"| {Dg} | {jk2} | {jcap} | {hi[jh]:.0f} (J={jh}) | {al[ja]:.0f} (J={ja}) | " + " / ".join(fr) +
           f" | {min(hi5.values()):.0f} |")
    # power law of min T_B3 (J >= J_K2) in Lam
    xs, ys = [], []
    for Dg in Dgs:
        Lm = lam_of(Dg); jk2 = int(np.ceil(np.sqrt(1.5*Lm)))
        v = [t for J, t in TB[Dg].items() if J >= jk2 and J >= 4]
        xs.append(np.log(Lm)); ys.append(np.log(min(v)))
    sl = np.polyfit(xs, ys, 1)
    pr(f"fit: min T_B3 (J >= J_K2) ~ {np.exp(sl[1]):.3g} Lam^{sl[0]:.3f} over Dg = {Dgs[0]}..{Dgs[-1]}")
    for G in (200.0, 450.0, 1000.0):
        Lst = np.exp((np.log(G) - sl[1])/sl[0])
        pr(f"  sign-free H_win holds (min T_B3 >= |Gamma_win| = {G:.0f} c2) for Lam >~ {Lst:.0f}, i.e. Dg >~ {Lst + 1.5:.0f}  [M: Gamma_win Delta-independent]")
    # ---------------------------------------------------------------- chains at 100 / 300 / 1000
    pr("\n## D. computed chain verdicts (parsed): per-J (min lam/N_s, uncorrected, t*_J upper)")
    c100 = parse_chain('y2_chain100_jhi.txt', r'\[EE nl min\]')
    c100s = parse_chain('y2_chain100_jhi_scale.txt', r'\[EE nl min\]')
    c300 = parse_chain('w2_chain_final_xc0.txt', r'Delta_g = 300:')
    c1000 = parse_chain('w2_chain_final_xc0.txt', r'Delta_g = 1000:')
    pr("| J | 100 jhi (held) | 100 jhi_scale (Gamma ~ m) | 300 conv A | 1000 conv A | p_emp = ln(t*1000/t*300)/ln(1000/300 Lam) | 300 adopted |")
    pr("|---|---|---|---|---|---|---|")
    pemp = {}
    ad = {J: v for J, v, _ in ADOPT300}
    for J in sorted(set(c100) | set(c300) | set(c1000)):
        f = lambda c: (f"{c[J][0]:+.3f} / t* {c[J][2]:.3g}" if J in c else "—")
        p = ""
        if J in c300 and J in c1000 and np.isfinite(c1000[J][2]) and np.isfinite(c300[J][2]) and c300[J][2] > 0 and J >= 8 and J != 28 and J != 32:
            pe = np.log(c1000[J][2]/c300[J][2])/np.log(lam_of(1000)/lam_of(300)); pemp[J] = pe; p = f"{pe:.2f}"
        pr(f"| {J} | {f(c100)} | {f(c100s)} | {f(c300)} | {f(c1000)} | {p} | {ad.get(J, '—')} |")
    pv = np.array(list(pemp.values()))
    pr(f"p_emp over J = 8..24 (conv A, own rows at both gaps): median {np.median(pv):.2f}, range {pv.min():.2f} .. {pv.max():.2f}"
       " (1000 J = 28, 32 excluded: pre-Jacobi-fix rows, EE_DG1000_LARGEJ §1)")
    # ---------------------------------------------------------------- P1: correction transfer from the adopted 300 values
    pr("\n## E. P1: adopted-300 corrections transferred at fixed (J, m/Lam): margin_J(Dg) = lambda0_J(Dg) - corr_J(300) (298.5/Lam)^p")
    corr = {}
    for J, lam, mf in ADOPT300:
        l0, mfl = lam0_at(300, J, mf)
        corr[J] = (max(l0 - lam, 0.0), mf, l0)
    pr("corr_J(300) = lambda0(300; state) - adopted:  " + ", ".join(f"J={J}: {c[0]:.3f}" for J, c in corr.items()))
    grid = [100, 125, 150, 175, 200, 225, 250, 275, 300, 400, 500, 1000, 2000, 5000]
    res = {}
    pr("| p | " + " | ".join(str(g) for g in grid) + " | threshold Dg* (min_J margin = 0) |")
    pr("|---|" + "---|"*len(grid) + "---|")
    for p in (2.0, 3.0, float(np.median(pv)), 4.0):
        row, mins = [], []
        for Dg in grid:
            Lm = lam_of(Dg)
            m = []
            for J, (c, mf, _) in corr.items():
                l0, _ = lam0_at(Dg, J, mf)
                m.append((l0 - c*(298.5/Lm)**p, J))
            mn = min(m)
            mins.append(mn[0]); row.append(f"{mn[0]:+.3f} (J={mn[1]})")
        # threshold by log-linear interpolation of the sign change
        thr = "none in grid"
        for (a, va), (b, vb) in zip(zip(grid[:-1], mins[:-1]), zip(grid[1:], mins[1:])):
            if va < 0 <= vb:
                thr = f"{a + (b - a)*(-va)/(vb - va):.0f}"
                break
        res[p] = dict(mins=mins, thr=thr)
        pr(f"| {p:.2f} | " + " | ".join(row) + f" | {thr} |")
    pr("(only J <= 32, where 300 anchors exist; J > 32 near the gap is the H_LJ region, section C.)")
    # ---------------------------------------------------------------- consistency: 300-anchored prediction at 100 vs the 100 chain
    pr("\n## F. consistency of P1 with the computed 100 (transfer) and 1000 chains")
    for J in (12, 14, 16, 18, 20, 24, 28):
        if J not in corr:
            continue
        c, mf, _ = corr[J]
        s = []
        for p in (2.0, 3.0, 4.0):
            l0, _ = lam0_at(100, J, mf)
            s.append(f"p={p:.0f}: {l0 - c*(298.5/lam_of(100))**p:+.2f}")
        s1 = [f"p={p:.0f}: {lam0_at(1000, J, mf)[0] - c*(298.5/lam_of(1000))**p:+.3f}" for p in (2.0, 3.0)]
        pr(f"  J={J}: predicted at 100 [{', '.join(s)}] vs 100 chain {c100[J][0] if J in c100 else '—'} (held) / "
           f"{c100s[J][0] if J in c100s else '—'} (Gamma~m);  predicted at 1000 [{', '.join(s1)}] vs 1000 chain {c1000[J][0] if J in c1000 else '—'}")
    # ---------------------------------------------------------------- G. Result Q binding state: s_d(Delta) from the per-order engine data
    # BELOW_FOLD2 §4.2 (H0 gauge, ++ J = 0 at the gap): per-order s_d contributions k = 0..4 at Delta = 300: +1.66/-3.95/+2.52/+0.27/+0.016;
    # Delta = 1000: +1.66/-4.21/+2.96/+0.09/+0.002.  Model per order: k = 0..3: s_k = s_k(inf) + b_k/Delta; k = 4: b_4/Delta^2.
    o300 = np.array([1.66, -3.95, 2.52, 0.27]); o1000 = np.array([1.66, -4.21, 2.96, 0.09])
    bk = (o300 - o1000)/(1/300 - 1/1000); sinf = o1000 - bk/1000
    b4 = 0.016*300**2
    sd = lambda D: sinf.sum() + bk.sum()/D + b4/D**2
    pr(f"\n## G. Result Q binding state (++/OO J = 0 at the gap): s_d(Delta) ~ {sinf.sum():.3f} {bk.sum():+.2f}/Delta + {b4:.0f}/Delta^2 "
       f"[per-order fit to BELOW_FOLD2 §4.2; M]; s_fail = 1.78 (delta = 0.2)")
    pr("  " + ", ".join(f"Dg={D}: s_d = {sd(D):.3f} (s_fail/s_d = {1.78/sd(D):.2f}); de-shift moment ratio 1.9e3/Lam^2 = {1.9e3/lam_of(D)**2:.3f}"
                         for D in (100, 150, 200, 250, 300, 500, 1000, 2000, 5000)))
    # ---------------------------------------------------------------- H. combined grid table
    pr("\n## H. combined table of margins vs Delta_g (most dangerous item per column)")
    pr("| Dg | Lam | baseline min lambda0 [M] | P1: min_{J<=32} margin, p = 2 / 2.83 / 4 [M] | window proxy min_{J<=32} |C*_25|/450 | H_LJ: min T_B3 / 250 (refined) , / 450 (sign-free) | Q: s_fail/s_d |")
    pr("|---|---|---|---|---|---|---|")
    p1 = {}
    for p in (2.0, float(np.median(pv)), 4.0):
        p1[p] = {}
        for Dg in Dgs:
            Lm = lam_of(Dg)
            p1[p][Dg] = min(lam0_at(Dg, J, mf)[0] - c*(298.5/Lm)**p for J, (c, mf, _) in corr.items())
    for Dg in Dgs:
        Lm = lam_of(Dg); jk2 = int(np.ceil(np.sqrt(1.5*Lm)))
        w = min(abs(v) for J, v in C25[Dg].items() if J <= 32 and J >= 4)
        tb = min(t for J, t in TB[Dg].items() if J >= jk2 and J >= 4)
        pr(f"| {Dg} | {Lm:.1f} | {base[Dg]:+.4f} | " + " / ".join(f"{p1[p][Dg]:+.2f}" for p in p1) + f" | {w/450:.2g} | {tb/250:.2g} , {tb/450:.2g} |"
           f" {1.78/sd(Dg):.2f} |")
    json.dump(dict(base=base, pemp=pemp, P1={str(k): v for k, v in res.items()}), open(os.path.join(HERE, 'cache', 'analyze.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'cache', 'analyze.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
