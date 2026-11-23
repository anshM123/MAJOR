"""pair-wise report of rh_t2_heft (1<->3 maps class s <-> t and 1243 <-> 1324, so these pairs are the symmetric units)."""
import sys, pickle
import numpy as np
import rh_common as R
from rh_t2_heft import J0S, NDISP, NFULL, DISP

import os
PFX = '' if os.environ.get('RH_MODE', 'deg') == 'j' else 'd'
PAIRS = {'s+t': ('s', 't'), '1243+1324': ('1243', '1324')}


def load(k, tag):
    return {o: pickle.load(open(f'{R.CACHE}/t2{PFX}_{tag}_k{k}_{o}.pkl', 'rb')) for o in ('ab', 'ba')}


def pairvals(D, k, e, pair):
    sgn = (-1.0)**k
    V = {}
    for j0 in J0S:
        v = 0
        for c in PAIRS[pair]:
            nmax = NDISP if c in DISP else NFULL
            v += np.sum(D['ab']['acc'][(e, j0, c)][:nmax + 1]) + sgn*np.sum(D['ba']['acc'][(e, j0, c)][:nmax + 1])
        V[j0] = v
    Sh = sum(np.sum(D['ab']['sh'][(e, c)]) + sgn*np.sum(D['ba']['sh'][(e, c)]) for c in PAIRS[pair])
    return V, Sh


def report(k, tag):
    D = load(k, tag)
    L = []
    pr = lambda s: (print(s), L.append(s))
    pr(f'pair report, tag {tag}, k = {k}.  V[X^(j)] = symmetrised term-wise value of the remainder after EFT order j-1 (V[X^(0)] = DT[W_h]).')
    pr('  EFT_j := V[X^(j)] - V[X^(j+1)] (term-wise value of the order-j EFT polynomial); TP_R - DT = -V[X^(R+1)].')
    pr('  j* = first order with EFT_j* != 0 (orders below are subcritical: term-wise zero, so V[X^(j)] is constant = DT for j <= j*).')
    for pair in PAIRS:
        for e in sorted(set(x[0] for x in D['ab']['acc'])):
            V, Sh = pairvals(D, k, e, pair)
            j0start = 5 if pair == 's+t' else 3          # below j0start the remainders carry float noise from large subcritical pieces
            js = [j for j in J0S[j0start:-1] if abs(V[j + 1] - V[j]) > 1e-2*abs(V[j])]
            jst = js[0] if js else None
            if jst is None:
                pr(f'  {pair:9s} tau {e[0]:3.0f} J {e[1]}: no non-zero EFT order up to 40'); continue
            DT = V[jst]
            plat = max(abs(V[j]/DT - 1) for j in range(max(j0start, jst - 3), jst + 1))
            lead = V[jst] - V[jst + 1]
            errs = {R1 - 1: abs(V[R1])/abs(DT) for R1 in range(jst + 1, 41)}
            sh_part = {}
            cs = sum(np.cumsum(D['ab']['sh'][(e, c)] + (-1.0)**k*D['ba']['sh'][(e, c)]) for c in PAIRS[pair])
            shconv = abs(cs[-1] - cs[80])/max(abs(cs[-1]), 1e-300)
            pick = [jst, jst + 2, jst + 4, jst + 8, jst + 12, jst + 16, jst + 24, 39]
            pr(f'  {pair:9s} tau {e[0]:3.0f} J {e[1]}: j* = {jst} (plateau spread {plat:.0e}); DT = {DT.real:+.5e}; DT/EFT_j* - 1 = {(DT/lead - 1).real:+.4e};'
               f' |TP_R - DT|/|DT| at R = ' + ', '.join(f'{R}:{errs[R]:.1e}' for R in pick if R in errs) +
               f';  min over R {min(errs.values()):.1e}, max {max(errs.values()):.1e}'
               f';  S_h = {Sh.real:+.5e} [(S_h+DT)/|DT| = {abs(Sh + DT)/abs(DT):.1e}; m-sum change 80->160: {shconv:.1e}]')
    open(f'{R.CACHE}/t2{PFX}_pairs_{tag}_k{k}.txt', 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    report(int(sys.argv[1]), sys.argv[2])
