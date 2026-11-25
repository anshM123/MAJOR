"""Build point set 'mine' in the sandbox copy: BATCH5's 'fix' set with every N5 record of J replaced by MY independently mapped
record (i3_basis2 constrained fit; E, reps, err from my map; other fields from the raw record)."""
import sys, os, glob, shutil, pickle, numpy as np
import i3_basis2 as B2
SBX = '/home/user/MAJOR/ads4/code/indep3/cache/sbx/ads4/code/b5fold/cache/pts'
J = B2.J; loc3 = sys.argv[2]
os.makedirs(f'{SBX}/mine', exist_ok=True)
for f in glob.glob(f'{SBX}/fix/*.pkl'):
    if not os.path.exists(f'{SBX}/mine/' + os.path.basename(f)):
        shutil.copyfile(f, f'{SBX}/mine/' + os.path.basename(f))
maps = {}
for Dl in (B2.D, B2.REF):
    solG = B2.L5.conserved(Dl, J)[1]; solH = B2.L5.conserved(3 - Dl, J)[1]
    L3 = B2.local(loc3, 'EE'); raw = B2.raw_records('EE')
    key = [k for k in L3 if abs(k[0] - Dl) < 1e-6][0]
    A, Bm, best, res = B2.fit(B2.stack(L3[key]), B2.stack(raw[key]), solG, solH); maps[Dl] = (A, Bm)
for kind in ('EE', 'C3H'):
    fn = f'{SBX}/mine/hc_pts_{kind}_J{J}_D{B2.D:g}.pkl'
    P = pickle.load(open(fn, 'rb')); n = 0
    for f in glob.glob(f'{B2.RAW}/N5_J{J}_*/n1_J{J}_*/hc_pts_{kind}_J{J}_D{B2.D:g}.pkl'):
        for k, rec in pickle.load(open(f, 'rb')).items():
            A, Bm = maps[round(k[0], 6)]; iA, iB = np.linalg.inv(A.T), np.linalg.inv(Bm)
            r = dict(rec); Rm = np.array([iA @ x @ iB for x in rec['reps']])
            r['reps'] = Rm; r['E'] = Rm.mean(0)
            r['err'] = np.sqrt((np.abs(Rm - Rm.mean(0))**2).sum(0)/((len(Rm) - 1)*len(Rm)))
            kk = [q for q in P if abs(q[0] - k[0]) < 1e-6 and abs(q[1] - k[1]) < 1e-6]
            assert len(kk) == 1, k
            P[kk[0]] = r; n += 1
    pickle.dump(P, open(fn, 'wb')); print(kind, 'replaced', n, 'records of', len(P))
