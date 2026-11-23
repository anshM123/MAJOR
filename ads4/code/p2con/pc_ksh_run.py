"""
p2con/pc_ksh_run.py (P2_CONTRACT.md §1.7): exact image and cost of the EE kernel's own O(r) term (pc_ksh) on the leading heavy block
of every sample state (entries as pc_run), de-shifted and costed exactly as in pc_run.  -> cache/ksh_<i0>_<i1>.json/.log
Usage: python3 pc_ksh_run.py i0 i1
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import pc_sym as Sy, pc_contract as Ct, pc_run as Rn, pc_ksh as K
from pc_gensym import qrat


def lead_modes(S, h, hp):
    T = Sy.block_tensor(S, h, hp, 0)
    return {k: sum(v.values()) for k, v in T.items()}


def main():
    R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
    i0, i1 = int(sys.argv[1]), int(sys.argv[2])
    tabN = Rn.get_tables(Rn.NODES, 'nodes')
    logf = open(os.path.join(HERE, 'cache', 'ksh_%d_%d.log' % (i0, i1)), 'w')
    out = []
    for r in R[i0:i1]:
        t0 = time.time()
        q = qrat(r['beta']); beta = float(2*sp.log(q)); m = r['m']; J = r['J']
        S = Sy.numeric_series(Sy.drift_free(Sy.load_q(q))[0], q)
        hs = ('same', 'opp') if J % 2 == 0 else ('odd',)
        shp = Ct.lead_shape(tabN['nus'], beta)
        Nh = {h: float(np.median((Ct.raw_image(Sy.block_tensor(S, h, h, 0), tabN, beta)/shp).real)) for h in hs}
        A = Rn.design(Rn.NODES, beta); AH = Rn.design(Rn.HOLD, beta)
        ent = {}
        for i, h in enumerate(hs):
            for hp in hs[i:]:
                lt = lead_modes(S, h, hp)
                nrm = np.sqrt(Nh[h]*Nh[hp])
                y = K.ksh_raw(Rn.NODES, beta, lt)/(2*np.pi*nrm)
                yH = K.ksh_raw(Rn.HOLD, beta, lt)/(2*np.pi*nrm)
                x, res = Rn.fitF(y, A, Rn.NODES)
                hold = float(np.max(np.abs(AH@x - yH))/max(np.max(np.abs(yH)), 1e-300))
                c = Rn.costs(x, beta, m); c.update(fitres=float(res), hold=hold, maxraw=float(np.max(np.abs(y))), x=[float(v) for v in x])
                ent['%s|%s' % (h, hp)] = c
        row = dict(J=J, mfl=r['mfl'], m=m, beta=beta, ent=ent, time=time.time() - t0)
        out.append(row)
        s = f"J={J:3d} m/L={r['mfl']:.3f} beta={beta:.4f} [{row['time']:.0f}s]  " + ' | '.join(
            f"{k}: max|raw| {v['maxraw']:.2e} fit {v['fitres']:.1e} hold {v['hold']:.1e} total {v['total']:+.3e} (X {v['totalX']:+.3e}) "
            f"kB {v['kB']:+.4f} kC2 {v['kC2']:+.4f} kA {v['kA']:+.3f} kC {v['kC']:+.3f} Q {v['Q']:+.3f} crem {v['crem']:+.2e}" for k, v in ent.items())
        print(s, flush=True); logf.write(s + '\n'); logf.flush()
        json.dump(out, open(os.path.join(HERE, 'cache', 'ksh_%d_%d.json' % (i0, i1)), 'w'))


if __name__ == '__main__':
    main()
