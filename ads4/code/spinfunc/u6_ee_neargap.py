"""
u6 (EE_FUNCTIONAL.md §5): the S1_S4 near-gap chain (q3_analysis: K2 linear chain, rows Delta = 336/400/800/1600, aligned off-diagonals,
c_off = m Gamma_off model, combined near-gap states) with the C^_2 crystal rows of EE (u5) in place of q^2's, for the J where EE rows
exist (env U6_JEE, default "16,18,20"); C^_3 rows (C3H, o_H q) unchanged.  Runs the SAME chain twice (q^2 rows / EE rows) and prints the
per-J minima side by side.  The helicity basis of each row is the rank-1 basis of the order-0 crystal matrix, identical for EE and q^2
(u1: E = q at the crystal point; u5 checks M0_EE = M0_q2 per point).
The local k = 2 de-shift of the chain is the spinning (1 + u)/2 with u = u_T (q^2) or u_EE (exact rational, u2) respectively; all nodes
used at Delta_g = 300, 1000 have nu >= 15 (far from the real zeros 1.17, 6.64 of the EE forward normalisation).
Usage: OMP_NUM_THREADS=1 U6_JEE=16,18,20 Q3_JEVEN=16,18,20 Q3_JODD=5 python3 u6_ee_neargap.py 300 1000 -> cache/u6_ee_neargap.txt
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import q3_analysis as Q
import u5_ee_points as U5
import u2_hclosed_struct as U2

JEE = [int(v) for v in os.environ.get('U6_JEE', '16,18,20').split(',')]
_q3_load = Q.A._load
MODE = {'ee': False}
_F = U2.numfuncs()
_uq = Q.u_spin
_uE = lambda v: (lambda a, b: (a - b)/(a + b))(_F['EE'][0](np.asarray(v, complex)), _F['EE'][1](np.asarray(v, complex)))


def _load(kind, J, sfx):
    if MODE['ee'] and kind == 'C2' and J in JEE:
        if sfx != '_k2':
            return {}
        return Q._strip(U5.load_all(J), J)
    return _q3_load(kind, J, sfx)


Q.A._load = _load
_orig_basis = Q.A.L.basis


def _basis(J, D, R, *a, **k):
    # helicity basis of each row (rank-1 basis of the order-0 crystal matrix, x >= 0.15): taken from the q^2 records in BOTH runs, so that
    # the EE and q^2 chains are evaluated in the same (physical) helicity basis; at m -> infinity M0_EE = M0_q2 (u5 diag), and the EE
    # reference rows only exist at the state nu's (x_ref <= 0.064), too few for the basis fit.
    Rq = dict(Q.A.L.records(J)); Rq.update(_q3_load('C2', J, '_k2')); Rq.update(_q3_load('C2', J, ''))
    return _orig_basis(J, D, Rq, *a, **k)


Q.A.L.basis = _basis
_shift_q = Q.A.shift


def _set_mode(ee):
    MODE['ee'] = ee
    Q.u_spin = _uE if ee else _uq                        # q3._shift_generic reads the module-level name u_spin at call time
    Q.ALIGN_LOG.clear()


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000]
    t0 = time.time()
    res = {}
    for ee in (False, True):
        _set_mode(ee)
        pr(f"\n######## C^_2 rows: {'EE (u5)' if ee else 'q^2 (q2/l8 caches)'};  J with EE rows: {JEE}")
        grids, out = Q.run(Dgs, pr)
        res[ee] = out
    pr("\n# side by side (per J: min lam/N_s with corrections | max |dF_+-|/sqrt(NN) | max |dF_hh|/N_h):")
    for Dg in Dgs:
        for J in JEE:
            a = res[False][Dg]['perJ'].get(J); b = res[True][Dg]['perJ'].get(J)
            if a and b:
                pr(f"   Delta_g = {Dg}, J = {J}: q^2 {a['min']:+.4f} | {a['off']:.4f} | {a['diag']:.4f}   EE {b['min']:+.4f} | {b['off']:.4f} | "
                   f"{b['diag']:.4f}   (uncorrected {a['unc']:+.4f})")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'u6_ee_neargap.txt'), 'w').write(buf.getvalue())
