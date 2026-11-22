"""
bj_conv (BELOW_FOLD_SPINJ.md §2.3): internal convergence of the closed-form engine images: per nu, m^2 x [log(I_a/I_b) - const]
(state/reference ratio, entry (0,0) and max over entries), between truncations (N: P-series terms; Kr: residual geometric order).
Usage: python3 bj_conv.py Delta J depth [patch]
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bj_eeq as BQ, bj_ee as BJ

if __name__ == '__main__':
    Dl, J, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    patch = len(sys.argv) > 4 and sys.argv[4] == 'p'
    st = pickle.load(open(BQ.fn(Dl, J, depth, patch), 'rb'))
    rf = pickle.load(open(BQ.fn(2e4, J, depth, patch), 'rb'))
    al = Dl - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
    keys = sorted(st['res'])
    kr = max(k for k in rf['res'])
    nus = st['res'][keys[0]]['nus']
    def ratio(k):
        I = np.array([[BJ.value(st['res'][k], 'EE2', a, c)[0]/BJ.value(rf['res'][kr], 'EE2', a, c)[0] for c in range(2)] for a in range(2)])
        return np.log(I)
    top = max(keys)
    ref = ratio(top)
    out = [f'# bj_conv ({Dl:g},{J}) depth {depth}: m^2 Re[log ratio(N,Kr) - log ratio(N={top[0]},Kr={top[1]})] minus its nu<1 mean; '
           f'entry (0,0) / max over entries; reference {kr}']
    sel_nu = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11.66, 13.33, 16.66, 20.83, 25]
    idx = [int(np.argmin(np.abs(nus - v))) for v in sel_nu]
    out.append('   (N,Kr)  ' + ' '.join(f'{nus[i]:>9.2f}' for i in idx))
    for k in keys:
        if k == top:
            continue
        d = m*m*(ratio(k) - ref).real
        d = d - d[:, :, nus < 1].mean(-1)[:, :, None]
        out.append(f'   {str(k):8s} ' + ' '.join(f'{d[0, 0, i]:+9.1f}' for i in idx))
        out.append(f'   {"max":8s} ' + ' '.join(f'{np.max(np.abs(d[:, :, i])):9.1f}' for i in idx))
    print('\n'.join(out))
    open(os.path.join(HERE, 'cache', 'bj_conv.txt'), 'a').write('\n'.join(out) + '\n')
