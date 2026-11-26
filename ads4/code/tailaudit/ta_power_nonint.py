"""tailaudit/ta_power_nonint.py [DV]: Theorem P / Lemma NB spot checks at non-integer twist (a = 998.5, 1234.567) and L = 18 (against
K_osc with the L0 = 18 table, TC_L0TAB=18), using ta_power_check.Shat.  Output: cache/ta_power_nonint.txt"""
import os, sys
os.environ['TC_L0TAB'] = '18'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.dont_write_bytecode = True
import ta_power_check as PC
import tc_power as P
import mpmath as mp
mp.mp.dps = 30
out = []
for a in (998.5, 1234.567):
    for b in (0, 4):
        for L in (18, 20):
            for sig, t in ((0.25, 0.25), (0.25, 45), (1.25, 10), (3.5, 300), (0.25, 510)):
                v = abs(PC.Shat(complex(sig, t), a, b, L))
                bd, _ = P.K_osc(sig, t, t + 0.01, 18)
                assert v <= bd
                out.append(f'a={a} b={b} L={L} beta={sig}+{t}i: |S^| = {v:.5f}  bound(L0=18) {bd:.4g}')
            for sig, t in ((2.5, 0), (1.75, 0.25), (0.25, 0.25)):
                v = abs(PC.Shat(complex(sig, t), a, b, L) - 1)
                bd = P.NB_bound(sig, t)
                assert v <= bd + 2e-11
                out.append(f'a={a} b={b} L={L} NB beta={sig}+{t}i: |S^-1| = {v:.4e}  bound {bd:.4e}')
out.append('all assertions passed')
print('\n'.join(out))
open(os.path.join(HERE, 'cache', 'ta_power_nonint.txt'), 'w').write('\n'.join(out) + '\n')
