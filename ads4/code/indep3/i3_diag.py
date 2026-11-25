"""diagnose per-key differences between my mapped J records and BATCH5's mapped set (data comparison)."""
import sys, numpy as np
sys.argv = [sys.argv[0]] + sys.argv[1:]
import i3_basis2 as B2
J = B2.J
def run(loc3):
    maps = {}
    for Dl in (B2.D, B2.REF):
        solG = B2.L5.conserved(Dl, J)[1]; solH = B2.L5.conserved(3 - Dl, J)[1]
        raw = B2.raw_records('EE'); L3 = B2.local(loc3, 'EE')
        key = [k for k in L3 if abs(k[0] - Dl) < 1e-6][0]
        A, Bm, best, res = B2.fit(B2.stack(L3[key]), B2.stack(raw[key]), solG, solH)
        maps[Dl] = (A, Bm)
    for kind in ('EE',):
        raw = B2.raw_records(kind); th = B2.keyed(B2.load(f'{B2.THEIRS}/hc_pts_{kind}_J{J}_D{B2.D:g}.pkl'))
        bad = []
        for key, rec in sorted(raw.items()):
            A, Bm = maps[key[0]]
            Em = np.linalg.inv(A.T) @ rec['E'] @ np.linalg.inv(Bm)
            d = B2.rel(Em, th[key]['E']); draw = B2.rel(rec['E'], th[key]['E'])
            if d > 1e-5: bad.append((key, round(rec['nu']/rec['m'], 4), f"{d:.1e}", f"raw-vs-theirs {draw:.1e}", th[key].get('method')))
        print(kind, len(bad), 'keys differ >1e-5'); [print('  ', b) for b in bad[:40]]
        ok = [B2.rel(np.linalg.inv(maps[k[0]][0].T) @ r['E'] @ np.linalg.inv(maps[k[0]][1]), th[k]['E']) for k, r in raw.items()]
        print('  median rel diff', np.median(ok))
run(sys.argv[2])
