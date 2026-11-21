"""
bh_k2forms (BELOW_FOLD3.md §2): the 512 c2 K2-vs-engine disagreement (BELOW_FOLD2 §3.3) is a resummation-form question.  Three ways to
turn the Laplace orders (p0, p2, p4) of state (A_k = p_k/p0) and reference (R_k) into D = log(state/ref):
  plain   : log(1 + A2 + A4) - log(1 + R2 + R4)                       (bg_k2check's form)
  cumulant: (A2 - A2^2/2 + A4) - (R2 - R2^2/2 + R4)                   (l8 logser form, truncated at eps^4)
  chain   : log(1 + (A2 - R2) + (A4 - R4 - A2 R2 + R2^2))              (l9_analysis.pred(K=2): the ratio expanded to eps^4 = Result D's rows)
Scalar (Delta = 400, J = 0; l1_crystal.laplace K1, K2 = A2, A4) vs the exact Mellin skeleton; EE k = 2 (w1_pts_EE_J{J}_k2_b rows, M0, M2,
M4; M6 = 0) vs the analytic-radial engine (bg_eeq.images, EE2; state K = 4 depth 5, reference K = 2 depth 3, as bg_k2check).
Usage: python3 bh_k2forms.py scalar | ee J [Delta=400 [numax=125]]   -> appends to cache/bh_k2forms.txt
"""
import os, sys, pickle
os.environ.setdefault('L1_DEG', '6')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np


def forms(A2, A4, R2, R4):
    return dict(plain=np.log(1 + A2 + A4) - np.log(1 + R2 + R4),
                cumulant=(A2 - A2*A2/2 + A4) - (R2 - R2*R2/2 + R4),
                chain=np.log(1 + (A2 - R2) + (A4 - R4 - A2*R2 + R2*R2)),
                K1chain=np.log(1 + (A2 - R2)))


def table(tag, nus, m, base, ref_lab, lap0, F, pr, i0s=(1, 0)):
    """base: the trusted D (exact / engine); lap0: Laplace K0 part (log p0 st/ref incl. logC - logK); F: forms dict."""
    for i0 in i0s:
        f = lambda v: m*m*((v - v[i0]).real)
        pr(f'# {tag}: m^2[D(nu) - D(nu0)], nu0 = {nus[i0]:.2f}; {ref_lab} | Laplace (form) minus {ref_lab}')
        B = f(base)
        for i, nu in enumerate(nus):
            pr(f'   nu = {nu:7.2f}: {ref_lab} {B[i]:+10.2f} | ' + '  '.join(f'{k} {f(lap0 + F[k])[i] - B[i]:+8.2f}' for k in F))


if __name__ == '__main__':
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    mode = sys.argv[1]
    if mode == 'scalar':
        d = pickle.load(open(os.path.join(HERE, 'cache', 'bh_k2scal.pkl'), 'rb'))
        import l1_crystal as Cr
        nus = d['nus'][d['nus'] < 90]; m = d['m']
        A2, A4, R2, R4, l0 = [], [], [], [], []
        for nu in nus:
            a, b = Cr.laplace(nu, 400.0, 0), Cr.laplace(nu, 2e4, 0)
            A2.append(a['K1']); A4.append(a['K2']); R2.append(b['K1']); R4.append(b['K2']); l0.append(a['logY0'] - b['logY0'])
        F = forms(*map(np.array, (A2, A4, R2, R4)))
        table('scalar Delta=400 J=0', nus, m, d['ex'][:len(nus)], 'exact', np.array(l0), F, pr, i0s=(3, 0))
        pr('   (A2, A4, R2, R4 at nu = 19.9 / 39.8: ' + ' '.join(f'{x[0].real:+.5f}/{x[3].real:+.5f}' for x in (A2, A4, R2, R4)) + ')')
    else:
        J = int(sys.argv[2])
        DS = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
        NUMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 125.0
        import bg_eeq as EQ
        d = pickle.load(open(os.path.join(HERE, 'cache', f'w1_pts_EE_J{J}_k2_b.pkl'), 'rb'))
        nus = np.array(sorted(set(k[1] for k in d if k[0] == DS and k[1] < NUMAX)))
        fn = os.path.join(HERE, 'cache', f'bh_k2forms_ee_J{J}.pkl' if DS == 400.0 else f'bh_k2forms_ee_J{J}_D{DS:g}.pkl')
        if os.path.exists(fn):
            st, ref = pickle.load(open(fn, 'rb'))
        else:
            ref = EQ.images(EQ.DREF, J, 2, 3, nus=nus, kinds=('EE2',), verbose=False)
            st = EQ.images(DS, J, 4, 5, nus=nus, kinds=('EE2',), verbose=False)
            pickle.dump((st, ref), open(fn, 'wb'))
        m = d[(DS, nus[0])]['m']
        A2, A4, R2, R4, l0 = [], [], [], [], []
        for nu in nus:
            r, rr = d[(DS, nu)], d[(EQ.DREF, nu)]
            p = [r[k].ravel()[0] for k in ('M0', 'M2', 'M4')]; q = [rr[k].ravel()[0] for k in ('M0', 'M2', 'M4')]
            A2.append(p[1]/p[0]); A4.append(p[2]/p[0]); R2.append(q[1]/q[0]); R4.append(q[2]/q[0])
            l0.append(np.log(p[0]) + r['logC'] - r['logK'] - np.log(q[0]) - rr['logC'] + rr['logK'])
        F = forms(*map(np.array, (A2, A4, R2, R4)))
        for Kt in (3, 4):
            Ist, _ = EQ.value(st, 'EE2', Kt=Kt); Iref, _ = EQ.value(ref, 'EE2')
            table(f'EE k=2 J={J} Delta={DS:g} (engine K={Kt})', nus, m, np.log(Ist/Iref), 'engine', np.array(l0), F, pr, i0s=(1, 0))
        pr('   (A2, A4, R2, R4 at nu = 19.9 / 39.8: ' + ' '.join(f'{x[0].real:+.5f}/{x[1].real:+.5f}' for x in (A2, A4, R2, R4)) + ')')
    open(os.path.join(HERE, 'cache', 'bh_k2forms.txt'), 'a').write('\n'.join(lines) + '\n')
