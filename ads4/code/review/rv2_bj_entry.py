"""rv2_bj_entry (REVIEW_UNIFIED_2026-11): helicity/entry-resolved version of bj_chain.calibrate's shape check.

bj_chain.calibrate tests only the least-squares scalar projection kappa of the 2x2 engine matrix against the z2 direct records.
Here: per entry (a, c) of the conserved-basis 2x2 matrix, the state/reference double ratio
    R_ac(nu) = [V_st,ac / V_ref,ac] / [I_st,ac / I_ref,ac]
(direct / engine), normalised by the calibrated rho, in c2 units  m^2 Re(R_ac/rho - 1), jackknife over replicas, CV applied.
Read-only: imports bj_chain / z3_nodes from spinfunc; reads their caches; writes nothing to spinfunc.
Usage: python3 rv2_bj_entry.py J depth N Kr [Nref]
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, os.path.abspath(SF))
sys.dont_write_bytecode = True
import numpy as np
import bj_chain as BC
import z3_nodes as Z3

J = int(sys.argv[1]); depth = int(sys.argv[2]); N = int(sys.argv[3]); Kr = int(sys.argv[4])
Nref = int(sys.argv[5]) if len(sys.argv) > 5 else 4
cal = BC.calibrate(J, depth=depth, N=N, Kr=Kr, Nref=Nref, pr=lambda s: print(s))
nus, Ist = cal['nus'], cal['Ist']; nur, Iref = cal['nur'], cal['Iref']
P = BC.direct(J); m = cal['m']; rho = cal['rho']; R = cal['R']
print('# entry-resolved double ratio m^2 Re(R_ac/rho - 1) +- jk  [c2];  entries (0,0) (0,1) (1,0) (1,1);  also |I_ac|/max|I| of engine state')
for v in sorted(set(k[1] for k in P if k[0] == BC.DST and (BC.DREF, k[1]) in P)):
    if v < 8.0:
        continue
    st, rf = P[(BC.DST, v)], P[(BC.DREF, v)]
    try:
        Is, Ir = BC.at(nus, Ist, v), BC.at(nur, Iref, v)
    except AssertionError:
        continue
    Vs = st['reps']*np.exp(st['logN']); Vr = rf['reps']*np.exp(rf['logN'])
    o = Z3.scalar_check(J, v, {BC.DST: st, BC.DREF: rf}, BC.DST)
    sr = (st['repsS']*np.exp(st['logNs']))/(rf['repsS']*np.exp(rf['logNs']))
    cvr = o['cv']*sr.mean()/sr                                     # per replica CV factor
    row = []
    for a in range(2):
        for c in range(2):
            rr = (Vs[:, a, c]/Vr[:, a, c])/(Is[a, c]/Ir[a, c])*cvr       # per replica
            d0 = m*m*(rr.mean()/rho - 1)
            jk = np.array([np.mean(np.delete(rr, i)) for i in range(R)])
            e = m*m*np.sqrt((R - 1)/R*np.sum(np.abs(jk - jk.mean())**2))/abs(rho)
            row.append('%+9.0f%+7.0fi +-%6.0f' % (d0.real, d0.imag, e))
    mag = np.abs(Is)/np.abs(Is).max()
    print('%7.3f: ' % v + ' | '.join(row) + '   |I|: ' + ' '.join('%.2f' % x for x in mag.ravel()))
