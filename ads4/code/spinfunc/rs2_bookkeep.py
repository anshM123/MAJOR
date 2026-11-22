"""
rs2_bookkeep (HRS_CLOSE.md §2-3): pole / branch / parity bookkeeping of the analytic-radial engine's integrands, order by order.

For the scalar CLS-D k = 2 integrand (l1_crystal / bg_geom.scalar_pieces, f = (v'-u')/(u'v')) and the EE k = 2 / k = 4 integrands
(bg_ee.ee_pieces, f = v'-u', k = 4 = k = 2 x Lemma-B R), assembled exactly as bg_eeq does (bg_geom.assemble x bg_geom.bfac), this
lists for every order a^k (k <= K) and every power (i nu/2)^j the monomials atoms^(E0 + i nu E1 + e) and checks the hypotheses of
the contiguity (integration-by-parts) argument of HRS_CLOSE §2:
  (B1) E1 (the nu-part of the leading exponents) identical for scalar, EE2, EE4, state and reference;
  (B2) E0 differences between kinds lie on the lattice: even in p1, p2, q1, q2 (integer in X = p^2), even in cf, integer in w, u, v;
  (B3) no monomial carries a negative power of s, sf, ca, sa (no pole at interior points / at the crystal angles s = 0, phi = 0,
       alpha = +-pi/2);
  (B4) parity: every monomial has  p1 = p2 = s (mod 2),  q1 = q2 (mod 2),  s + sf + sa = 0 (mod 2) (transverse reflection),
       cf + sa = 0 (mod 2) ((phi, psi) -> (pi - phi, pi - psi)); in particular no monomial with sa odd and s, sf even;
  (B5) every leading exponent has E1 != 0 (non-integer for real nu != 0: the twisted cycles are canonically regularised).
It also reports the smallest real exponent per factor (convergence of the integrals on D without continuation).
Usage: python3 rs2_bookkeep.py Delta J K   -> cache/rs2_bookkeep_D{Delta}_J{J}_K{K}.txt (+ .pkl)
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S, bg_geom as GG, bg_ee as BE

C = os.path.join(HERE, 'cache')
NAMES = S.NAMES


def collect(A, Bf, K):
    """all monomials (exponent rows, |coef|, order k, nu-power j) of corr x B-correction."""
    rows = []
    for j, cj in A['corr'].items():
        Pk = S.lmul(cj, Bf.c, K)
        for k in range(K + 1):
            if len(Pk[k]):
                E = Pk[k].E.astype(np.int64); Cc = np.abs(Pk[k].C)
                rows.append((k, j, E, Cc))
    return rows


def check(name, A, rows, out):
    E0, E1 = A['E0'], A['E1']
    out.append(f'== {name}')
    out.append('   E0 = ' + ' '.join(f'{n}:{v:+.3f}' for n, v in zip(NAMES, E0)))
    out.append('   E1 = ' + ' '.join(f'{n}:{v:+.3f}' for n, v in zip(NAMES, E1)))
    allE = np.vstack([r[2] for r in rows]); allC = np.concatenate([r[3] for r in rows])
    out.append(f'   monomials: {len(allC)} (orders k = 0..{max(r[0] for r in rows)}, nu-powers j = 0..{max(r[1] for r in rows)})')
    out.append('   min exponent shift  ' + ' '.join(f'{n}:{allE[:, i].min():+d}' for i, n in enumerate(NAMES)))
    out.append('   max exponent shift  ' + ' '.join(f'{n}:{allE[:, i].max():+d}' for i, n in enumerate(NAMES)))
    res = {}
    neg = {n: int((allE[:, NAMES.index(n)] < 0).sum()) for n in ('s', 'sf', 'ca', 'sa')}
    res['B3_negative_powers'] = neg
    p1, p2, q1, q2, s_, w, cf, sf, ca, sa, u, v = (allE[:, i] % 2 for i in range(12))
    viol = {
        'p1=p2=s': int(((p1 != p2) | (p1 != s_)).sum()),
        'q1=q2': int((q1 != q2).sum()),
        's+sf+sa even': int(((s_ + sf + sa) % 2 != 0).sum()),
        'cf+sa even': int(((cf + sa) % 2 != 0).sum()),
        'sa odd & s,sf even': int(((sa == 1) & (s_ == 0) & (sf == 0)).sum()),
    }
    res['B4_violations'] = viol
    # weights of the violating classes relative to the largest coefficient of the same order
    for key, m in (('s+sf+sa even', (s_ + sf + sa) % 2 != 0), ('sa odd & s,sf even', (sa == 1) & (s_ == 0) & (sf == 0)),
                   ('cf+sa even', (cf + sa) % 2 != 0)):
        if m.any():
            out.append(f'   ! {key}: {m.sum()} monomials, max |c|/max = {allC[m].max()/allC.max():.2e}')
    # surviving class (what the engine integrates: s, sf, sa even)
    surv = (s_ == 0) & (sf == 0) & (sa == 0)
    Es = allE[surv]
    lat = {'p,q even': int(((Es[:, 0:4] % 2) != 0).any(1).sum()), 'cf even': int((Es[:, 6] % 2 != 0).sum())}
    res['surviving_lattice_violations'] = lat
    res['n_surv'] = int(surv.sum())
    # convergence on D (real parts of the full exponents), surviving monomials
    tot = Es + E0[None, :]
    conv = dict(radial_min=[float(tot[:, i].min()) for i in range(4)], w_min=float(tot[:, 5].min()), cf_min=float(tot[:, 6].min()),
                u_min=float((tot[:, 10]).min()), v_min=float(tot[:, 11].min()))
    res['D_convergence'] = conv
    out.append(f'   B3 negative powers of s, sf, ca, sa: {neg}')
    out.append(f'   B4 parity violations: {viol}')
    out.append(f'   surviving (s, sf, sa even): {surv.sum()} monomials; lattice violations among them: {lat}')
    out.append('   D convergence (surviving, Re exponent; need radial > -2 [p-exponent], w > -1, cf > -1, u, v > -1/2): '
               + ', '.join(f'{k}={np.round(v, 3)}' for k, v in conv.items()))
    # the same with round-off residues removed (|c| < 1e-10 of the largest coefficient of the same (order, nu-power) block)
    keepE = []
    for (k, j, E, Cc) in rows:
        m = (Cc > 1e-10*Cc.max()) & (E[:, 4] % 2 == 0) & (E[:, 7] % 2 == 0) & (E[:, 9] % 2 == 0)
        keepE.append(E[m])
    tf = np.vstack(keepE) + E0[None, :]
    convf = dict(radial_min=[float(tf[:, i].min()) for i in range(4)], w_min=float(tf[:, 5].min()), cf_min=float(tf[:, 6].min()),
                 u_min=float(tf[:, 10].min()), v_min=float(tf[:, 11].min()))
    res['D_convergence_filtered'] = convf
    out.append('   D convergence without round-off residues (|c| > 1e-10 max per block): '
               + ', '.join(f'{k}={np.round(v, 3)}' for k, v in convf.items()))
    zeroE1 = [NAMES[i] for i in (0, 1, 2, 3, 5, 6, 10, 11) if abs(E1[i]) < 1e-12]
    res['B5_E1_zero_atoms'] = zeroE1
    out.append(f'   B5 divisor atoms with E1 = 0: {zeroE1}')
    return res


def run(Dl, J, K):
    t0 = time.time()
    out = [f'rs2_bookkeep Delta={Dl} J={J} K={K}']
    Bf = GG.bfac(K)
    res = {}
    if J == 0:
        logs, polys, G, nm = GG.scalar_pieces(K)
        A = GG.assemble(logs, polys, K)
        res['scalar'] = (A['E0'], A['E1'], check('scalar k=2 (f = (v\'-u\')/(u\'v\'))', A, collect(A, Bf, K), out))
    else:
        import bg_scalarJ  # noqa: F401  (scalar J > 0 uses P_J(t); not re-assembled here)
    P = BE.ee_pieces(K + 1, Dl, J)
    out.append(f'   EE pieces built {time.time() - t0:.0f}s')
    na, nc = len(P['pref']), len(P['pref'][0])
    for a in range(na):
        for c in range(nc):
            for kind in ('EE2', 'EE4'):
                polys = [('vmu', P['vmu']), ('pref', P['pref'][a][c])] + ([('R', P['R'])] if kind == 'EE4' else [])
                A = GG.assemble(P['logs'], polys, K)
                res[(kind, a, c)] = (A['E0'], A['E1'], check(f'{kind} ({a},{c})', A, collect(A, Bf, K), out))
    # R alone: leading form (expect 4 a^4 C^2 (1 + O(a^4)))
    R = P['R']
    out.append(f'== Lemma-B R as an a-series: leading a-power d = {R.d}; nonzero orders: '
               + str([k for k, p in enumerate(R.c) if len(p)]))
    pts = S.atom_points(6, 5)
    X1, X2, Y1, Y2 = pts[:, 0]**2, pts[:, 1]**2, pts[:, 2]**2, pts[:, 3]**2
    Cc = X1*Y2 + X2*Y1 - 2*pts[:, 4]*pts[:, 0]*pts[:, 1]*pts[:, 2]*pts[:, 3]*pts[:, 7]
    ratio = R.c[0].ev(pts)/(4*Cc**2)
    out.append(f'   R.c[0] / (4 C^2) at 6 random points: {np.round(ratio.real, 12)} (imag max {np.abs(ratio.imag).max():.1e})')
    # differences of leading exponents between kinds
    keys = list(res.keys())
    out.append('== E0 / E1 differences (relative to the first kind)')
    E0r, E1r = res[keys[0]][0], res[keys[0]][1]
    for k in keys[1:]:
        dE0 = res[k][0] - E0r; dE1 = res[k][1] - E1r
        out.append(f'   {k}: dE0 = ' + ' '.join(f'{n}:{v:+.3f}' for n, v in zip(NAMES, dE0) if abs(v) > 1e-12)
                   + f' | max|dE1| = {np.abs(dE1).max():.1e}')
    out.append(f'done {time.time() - t0:.0f}s')
    txt = '\n'.join(out)
    print(txt, flush=True)
    tag = 'D%g_J%d_K%d' % (Dl, J, K)
    open(os.path.join(C, f'rs2_bookkeep_{tag}.txt'), 'w').write(txt + '\n')
    pickle.dump({k: (v[0], v[1], v[2]) for k, v in res.items()}, open(os.path.join(C, f'rs2_bookkeep_{tag}.pkl'), 'wb'))


if __name__ == '__main__':
    run(float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
