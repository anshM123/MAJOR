import dc_guard, pickle, os, numpy as np
import dc_pipe as DP
d = pickle.load(open(os.path.join(dc_guard.SPIN, 'cache', 'pw_gamma.pkl'), 'rb'))
for est in ('cv', 'dir'):
    for tags in (['_w'], ['_wP'], ['_none']):
        for kf in ((True,) if tags == ['_none'] else (False, True)):
            o = DP.get(100, 1100.0, est, tags, kf, redo=True)
            for k in d:
                if k[0] == 100 and k[1] == 1100.0 and k[4] == est and k[5] == tuple(tags):
                    r = d[k]
                    same = len(r['nus']) == len(o['nus']) and np.allclose(r['nus'], o['nus']) and np.allclose(r[('C2', '++')]['y'], o[('C2', '++')]['y'])
                    if same: print('   reproduces stored', k)
