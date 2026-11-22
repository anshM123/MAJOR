"""
qf_eeq (RESULT_Q_FWD_CLOSE.md): production driver for the analytic-radial EE k = 4 images needed to close H_OO^Q and H_fwd^Q.
Read-only imports: po_ee_odd.images (parity-odd symmetric <TTO> structure, 1 x 1) and bg_eeq.images (parity-even structures, all
(a, c) entries).  Nothing of po_* / bg_* / bh_* is modified; results go to cache/qf_eeq_{odd|even}_D{Delta}_J{J}_K{K}.pkl.
  mode odd : parity-odd structure on both sides (OO at even J; +- at odd J)
  mode even: parity-even conserved structures (J >= 4: 2 x 2 prefactor matrix -> ++ and +- after the helicity basis)
Usage: python3 qf_eeq.py odd|even Delta J K depth
"""
import os, sys, pickle, time, resource
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
C = os.path.join(HERE, 'cache')


def fn(mode, Dl, J, K):
    return os.path.join(C, 'qf_eeq_%s_D%g_J%d_K%d.pkl' % (mode, Dl, J, K))


def load(mode, Dl, J, K):
    """qf_ cache, else the equivalent existing po_/bh_/bg_ cache (read-only)."""
    p = fn(mode, Dl, J, K)
    if os.path.exists(p):
        return pickle.load(open(p, 'rb'))
    alt = {'odd': ['po_eeq_D%g_J%d_K%d.pkl'], 'even': ['bh_eeq_D%g_J%d_K%d.pkl']}[mode]
    for a in alt:
        q = os.path.join(C, a % (Dl, J, K))
        if os.path.exists(q):
            return pickle.load(open(q, 'rb'))
    if mode == 'even' and K == 2 and abs(Dl - 2e4) < 1e-6:          # bg_eeq stores the reference inside the state pickle
        import glob
        for q in sorted(glob.glob(os.path.join(C, 'bg_eeq_D*_J%d_K4.pkl' % J))):
            return pickle.load(open(q, 'rb'))['ref']
    raise FileNotFoundError(p)


if __name__ == '__main__':
    mode, Dl, J, K, depth = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    t0 = time.time()
    if mode == 'odd':
        import po_ee_odd as PO
        rec = PO.images(Dl, J, K, depth)
    else:
        import bg_eeq as EQ
        rec = EQ.images(Dl, J, K, depth, kinds=('EE4',))
    rec['qf_mode'] = mode
    p = fn(mode, Dl, J, K)
    pickle.dump(rec, open(p + '.tmp', 'wb')); os.replace(p + '.tmp', p)
    print(f'qf_eeq {mode} D={Dl} J={J} K={K} depth={depth}: done {time.time() - t0:.0f}s, maxrss '
          f'{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f} GB', flush=True)
