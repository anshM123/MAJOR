"""
lj_rs_book (HLJ_CLOSE.md §2.4): HRS_CLOSE's bookkeeping hypotheses (B1-B4) checked on the spin-J closed-form engine's integrand, so that
Prop. RS (crystal cycle = c_Gamma(nu) x engine, order by order, c_Gamma common to state, reference and the J = 0 scalar/EE) applies at
J = 28, 32 (and any cached bj pieces).
  B1: E1 (the nu-part of the leading exponents) identical to the J = 0 EE2 one (bg_eeq_D300_J0_K4.pkl) and to the reference's;
  B2: E0 - E0(J = 0 EE2) on the integer lattice of the atoms (even in p, q: integer in X, Y);
  B3: no negative power of s, sf, ca, sa in any monomial of any order (no interior / crystal-angle poles);
  B4: parity of surviving monomials: p1 = p2 = s, q1 = q2 (mod 2), s + sf + sa even and cf + sa even -- violators listed by weight.
Usage: python3 lj_rs_book.py J Delta depth N Kr [p]   -> cache/lj_rs_book_J<J>_D<Delta>_d<depth>.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

C = os.path.join(HERE, 'cache')


def main(J, Dl, depth, N, Kr, patch=False):
    import bj_ee as BJ, bg_ser as S, bg_geom as GG, lj_rs_eeq as LE
    P = pickle.load(open(LE.resolve(Dl, J, depth, patch, pieces=True), 'rb'))
    j0 = pickle.load(open(os.path.join(C, 'bg_eeq_D300_J0_K4.pkl'), 'rb'))['st']['res'][('EE2', 0, 0)]
    E0r, E1r = np.asarray(j0['E0']), np.asarray(j0['E1'])
    pr_ = BJ.pref_sum(P, N, Kr)
    Bf = GG.bfac(Kr)
    lines = [f'# lj_rs_book J={J} Delta={Dl:g} depth={depth} N={N} Kr={Kr}: atoms {S.NAMES}']
    lines.append(f'#   J=0 EE2 reference: E0 {E0r.tolist()}  E1 {E1r.tolist()}')
    for ia in range(len(pr_)):
        for ic in range(len(pr_[0])):
            A = GG.assemble(P['logs'], [('vmu', P['vmu']), ('pref', pr_[ia][ic])], Kr)
            dE0 = A['E0'] - E0r; dE1 = A['E1'] - E1r
            b1 = np.max(np.abs(dE1))
            lat = dE0.copy(); lat[:4] = lat[:4]/2.0
            b2 = np.max(np.abs(lat - np.round(lat)))
            neg = 0; nmon = 0; viol_w = 0.0; tot_w = 0.0; mx = 0.0
            for j, cj in A['corr'].items():
                Pk = S.lmul(cj, Bf.c, Kr)
                for k in range(Kr + 1):
                    E, Cc = Pk[k].E.astype(int), Pk[k].C
                    if not len(Cc):
                        continue
                    nmon += len(Cc)
                    neg += int(np.sum(np.any(E[:, [4, 7, 8, 9]] < 0, axis=1)))
                    par = ((E[:, 0] - E[:, 1]) % 2 != 0) | ((E[:, 0] - E[:, 4]) % 2 != 0) | ((E[:, 2] - E[:, 3]) % 2 != 0) | \
                          ((E[:, 4] + E[:, 7] + E[:, 9]) % 2 != 0) | ((E[:, 6] + E[:, 9]) % 2 != 0)
                    w = np.abs(Cc)
                    mx = max(mx, w.max())
                    viol_w = max(viol_w, w[par].max() if par.any() else 0.0)
            lines.append(f'  entry ({ia},{ic}): B1 max|E1 - E1(J=0)| = {b1:.1e};  E0 - E0(J=0) = {np.round(dE0, 6).tolist()} '
                         f'(B2 off-lattice {b2:.1e});  B3 monomials with negative s/sf/ca/sa power: {neg} of {nmon};  '
                         f'B4 largest parity-violating |c| / max |c| = {viol_w/max(mx, 1e-300):.1e}')
    out = '\n'.join(lines)
    print(out)
    open(os.path.join(C, 'lj_rs_book_J%d_D%g_d%d.txt' % (J, Dl, depth)), 'w').write(out + '\n')


if __name__ == '__main__':
    J, Dl, depth, N, Kr = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    main(J, Dl, depth, N, Kr, len(sys.argv) > 6 and sys.argv[6] == 'p')
