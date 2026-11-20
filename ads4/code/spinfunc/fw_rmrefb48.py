"""fw_rmrefb48: remove the 6 J = 48 refb records of cache/fw_k2_EE_J48.pkl that were computed BEFORE the structure-basis patch
(unpatched conserved basis + ee1000 shadow; identified by their shadow_resid = the ee1000 (2e4, 48) value), so that fw_rows recomputes
them in the patched basis.  Only this agent's cache is touched."""
import os, pickle, numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
f = os.path.join(C, 'fw_k2_EE_J48.pkl')
R = pickle.load(open(f, 'rb'))
old = pickle.load(open(os.path.join(C, 'ee1000_shadow_k2_A.pkl'), 'rb'))[(20000.0, 48)][2]
drop = [k for k, r in R.items() if k[0] == 20000.0 and np.allclose(np.ravel(r['shadow_resid']), np.ravel(old), rtol=1e-12, atol=0)
        and k[1] > 5000]
for k in drop:
    R.pop(k)
pickle.dump(R, open(f + '.tmp', 'wb')); os.replace(f + '.tmp', f)
print('removed', len(drop), 'refb records:', [round(k[1], 2) for k in drop])
