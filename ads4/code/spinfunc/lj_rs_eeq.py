"""
lj_rs_eeq (HLJ_CLOSE.md): bj_eeq.run (closed-form spin-J EE engine images, unchanged code) with its two cache files redirected to
cache/lj_rs_bj_eeq_*.pkl / cache/lj_rs_bj_pieces_*.pkl (no existing cache is read for writing or overwritten).  Existing bj_ pieces of the
same (Delta, J, depth, basis) are NOT reused (they would be copied, not recomputed, if present: set LJ_REUSE=1).
Usage: [BJ_PATCH=1] BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 python3 lj_rs_eeq.py Delta J depth N[,N2] [nus-spec]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import bj_eeq as BQ

_fn0 = BQ.fn


def fn(Dl, J, depth, patch=None):
    p = _fn0(Dl, J, depth, patch)
    return os.path.join(os.path.dirname(p), 'lj_rs_' + os.path.basename(p))


def resolve(Dl, J, depth, patch=None, pieces=False):
    """existing bj_ file if present, else the lj_rs_ one (for readers)."""
    for p in (_fn0(Dl, J, depth, patch), fn(Dl, J, depth, patch)):
        q = p.replace('bj_eeq_', 'bj_pieces_') if pieces else p
        if os.path.exists(q):
            return q
    return None


if __name__ == '__main__':
    BQ.fn = fn
    Dl, J, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Ns = [int(v) for v in sys.argv[4].split(',')]
    nus = BQ.grid(sys.argv[5] if len(sys.argv) > 5 else 'win')
    BQ.run(Dl, J, depth, Ns, nus, ('EE2',))
