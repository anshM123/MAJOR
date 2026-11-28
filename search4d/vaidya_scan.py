"""Physical control scan: concentric annuli in Vaidya-AdS4 (M=1).

For a pool of boundary edges (r_p, t_p) on a spacelike boundary slice
t = t0 + kappa (r - r0), compute every cap branch and every tube branch
(tubes found from BOTH ends as a branch-completeness check), take the minimal
branch, then:
  * K-slacks (SSA/WM instances of Lemma B) on the physical data,
  * SSA/WM/MMI on all unions of intervals for random sub-pools,
  * all 372 facets for random assignments of annuli to 5 parties,
reporting min margin per facet orbit.
Usage: python3 -I vaidya_scan.py dlt t0 kappa npool seed out.json
"""
import sys, os, json, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vaidya4 import Vaidya
from hrt_pieces import Pieces, continuation
from abstract_axis import F, ORB, entropy_vector, entropy, random_labels

INF = 1e300


def main():
    dlt = float(sys.argv[1]); t0 = float(sys.argv[2]); kap = float(sys.argv[3])
    npool = int(sys.argv[4]); seed = int(sys.argv[5]); outp = sys.argv[6]
    rng = np.random.default_rng(seed)
    nec = not (len(sys.argv) > 7 and sys.argv[7] == 'noNEC')
    V = Vaidya(M=1.0, dlt=dlt, NEC=nec)
    P = Pieces(V, cap_grid=(220, 200))
    T = time.time(); npropagated = [0]
    while True:   # well-separated edges (ratio >= 1.07, gap >= 0.03): thinner tubes are near-boundary
        radii = np.sort(np.exp(rng.uniform(np.log(0.25), np.log(2.4), npool)))
        if np.all(np.diff(np.log(radii)) > np.log(1.07)) and np.all(np.diff(radii) > 0.03):
            break
    r0 = radii[0]
    times = t0 + kap * (radii - r0)
    P.cap_scan(vrange=(min(times) - 7.0, max(times) + 0.5), zrange=(0.02, 6.0))
    c = np.full(npool, INF); ncap = []
    for i in range(npool):
        sols = P.caps(radii[i], times[i])
        ncap.append(len(sols))
        if sols:
            c[i] = min(r[2] for x, r in sols)
    L = np.full((npool, npool), INF); nbr = np.zeros((npool, npool), int); mism = 0.0; onesided = []
    allsol = {}; tried = {}
    for i in range(npool):
        for j in range(npool):
            if i != j:
                allsol[i, j] = P.tubes(radii[i], times[i], radii[j], times[j])
        # branch propagation: continue every solution from target j to the other targets
        # (catches branches the grid misses near folds)
        fun = P.tube_fun(radii[i], times[i])
        changed = True; sweeps = 0
        while changed and sweeps < 2:
            changed = False; sweeps += 1
            for j in range(npool):
                if j == i:
                    continue
                for jj in (j - 1, j + 1):
                    if jj < 0 or jj >= npool or jj == i:
                        continue
                    if len(allsol[i, jj]) >= len(allsol[i, j]):
                        continue
                    for x, rr in list(allsol[i, j]):
                        if any(abs(rr[2] - q[2]) < 1e-9 for q in tried.get((i, jj), [])):
                            continue
                        tried.setdefault((i, jj), []).append(rr)
                        res = continuation(fun, x, (radii[j], times[j]), (radii[jj], times[jj]), h0=0.1, hmin=2e-3)
                        if res is None:
                            continue
                        if not any(abs(res[1][2] - q[2]) < 1e-6 for _, q in allsol[i, jj]):
                            allsol[i, jj] = allsol[i, jj] + [res]; changed = True; npropagated[0] += 1
    for i in range(npool):
        for j in range(npool):
            if i == j:
                continue
            sols = allsol[i, j]
            nbr[i, j] = len(sols)
            a = min((r[2] for x, r in sols), default=INF)
            if i < j:
                L[i, j] = a
            else:
                b = L[j, i]
                if (a < INF) != (b < INF):
                    mism = max(mism, 1.0 if min(a, b) < c[i] + c[j] else 0.5)
                    onesided.append([int(j), int(i), float(min(a, b) - c[i] - c[j])])
                elif a < INF:
                    mism = max(mism, abs(a - b))
                L[j, i] = L[i, j] = min(a, b)
    tpieces = time.time() - T
    # K-slacks with L' = min(L, c+c)
    Lp = np.minimum(L, c[:, None] + c[None, :])
    ks = []
    for i, j, k in itertools.combinations(range(npool), 3):
        ks.append(Lp[i, k] + c[j] - Lp[i, j] - c[k]); ks.append(Lp[i, k] + c[j] - c[i] - Lp[j, k])
    for i, j, k, l in itertools.combinations(range(npool), 4):
        ks.append(Lp[i, k] + Lp[j, l] - Lp[i, j] - Lp[k, l]); ks.append(Lp[i, k] + Lp[j, l] - Lp[i, l] - Lp[j, k])
    kmin = float(min(ks)) if ks else None
    # facets on random sub-pools / labellings
    worst = np.full(len(ORB), np.inf); worst_info = [None] * len(ORB)
    nvec = 0
    for trial in range(400):
        m = int(rng.integers(5, min(9, npool) + 1))
        good = np.nonzero(c < INF)[0]   # edges whose cap was not found are excluded
        if len(good) < 5:
            break
        m = min(m, len(good))
        sub = np.sort(rng.choice(good, m, replace=False))
        Ls = L[np.ix_(sub, sub)]; cs = c[sub]
        for _ in range(10):
            lab = random_labels(m, rng)
            v = entropy_vector(lab, Ls, cs)
            if np.any(v >= INF / 2):
                continue
            mg = F @ v; nvec += 1
            for o in range(len(ORB)):
                if o == 2:
                    continue
                k = int(np.argmin(mg[ORB[o]]))
                if mg[ORB[o]][k] < worst[o]:
                    worst[o] = mg[ORB[o]][k]
                    worst_info[o] = dict(sub=sub.tolist(), lab=lab, facet=int(ORB[o][k]))
    res = dict(NEC=nec, dlt=dlt, t0=t0, kappa=kap, radii=radii.tolist(), times=times.tolist(), c=c.tolist(),
               ncap_branches=ncap, L=np.where(L < INF, L, None).tolist(), tube_branches=nbr.tolist(),
               two_end_mismatch=mism, one_sided=onesided, n_propagated=npropagated[0], K_min_slack=kmin, n_vectors=nvec,
               min_margin_per_orbit=[None if not np.isfinite(w) else float(w) for w in worst],
               worst_info=worst_info, seconds_pieces=tpieces, seconds_total=time.time() - T)
    json.dump(res, open(outp, 'w'), indent=1)
    print(json.dumps({k: res[k] for k in ('NEC', 'dlt', 't0', 'kappa', 'ncap_branches', 'two_end_mismatch', 'one_sided', 'n_propagated', 'K_min_slack',
                                          'n_vectors', 'min_margin_per_orbit', 'seconds_total')}))


if __name__ == '__main__':
    main()
