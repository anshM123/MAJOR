"""
fw_points (FOLD_WINDOW_LARGEJ.md): direct all-orders crystal points (kinds EE = C^'_EE w0 integrand, C3H = C^_3 o_H q; scalar on the
same samples) at large J near the gap, state and heavy reference Delta_ref = 2e4 on COMMON RANDOM NUMBERS -- z2_points (EE_DG300_SMALLNU)
unchanged in method (n3 hybrid contour above the fold, nu_f = 13 below 13; proposal adapted on the state, re-used for the reference; Q_J
overflow rescaling QSC), with two in-process precision fixes needed at J >~ 40 (fw_frame):
  (1) spin-J contraction: Fourier coefficients below the round-off floor zeroed (edge-mode noise ~ eps^2 4^J);
  (2) scalar Q_J by the Legendre recurrence instead of the power basis (cancellation ~ 2^J);
  (3) conserved <TTO> null space: over/underflowing random sample points skipped (needed at J >= 150).
Caches: cache/fw_pts_<kind>_J<J>.pkl keyed (Delta, round(nu, 9)) in z2/n1_run format; raw: cache/fw_raw.pkl.  Resumable (finished nodes
skipped).  Usage: OMP_NUM_THREADS=1 [FW_M=11 FW_R=8] python3 fw_points.py J nu1,nu2,... Dstate
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ['Z2_M'] = os.environ.get('FW_M', '11')
os.environ['Z2_R'] = os.environ.get('FW_R', '8')
import z2_points as Z2
import fw_frame as FW

C = Z2.C
TAG = os.environ.get('FW_TAG', '')


STATE = [None]          # state Delta of this run: its reference points go to a per-state file (the z2 key (Delta, nu) does not
                        # identify the state, and two rows of the same J sharing a nu overwrote each other's reference -- §2.6)


def fn(kind, J):
    if STATE[0] is not None:
        return os.path.join(C, 'fw%s_pts_%s_J%d_D%g.pkl' % (TAG, kind, J, STATE[0]))
    return os.path.join(C, 'fw%s_pts_%s_J%d.pkl' % (TAG, kind, J))


Z2.fn = fn
Z2.ZTAG = 'fw%s' % TAG          # raw cache -> cache/z2fw_raw.pkl is avoided below
_save = Z2._save


def _save_fw(path, key, val):
    if path.endswith('z2fw%s_raw.pkl' % TAG):
        path = os.path.join(C, 'fw%s_raw.pkl' % TAG)
    return _save(path, key, val)


Z2._save = _save_fw


def setup():
    import n1_pw as PW
    FW.patch_qrec_pw(PW, HERE)          # QSC + recurrence (z2.patch_qscale then returns early)
    Z2.patch()
    FW.install_pw(PW)
    import l1_crystal as Cr
    FW.patch_qrec_tps(Cr, HERE)         # TPS copy (saddle/flow set-up uses it for J > 0 scalar pieces)
    import l5_spin as S
    FW.patch_S(S)
    FW.patch_conserved_clean(S)
    FW.patch_conserved_clean(PW.S)       # J >~ 150: over/underflowing sample points skipped (null space unchanged)


if __name__ == '__main__':
    setup()
    J = int(sys.argv[1])
    nus = [float(v) for v in sys.argv[2].split(',')]
    Ds = [float(v) for v in sys.argv[3].split(',')]
    STATE[0] = Ds[0]
    if len(Ds) == 1:
        Ds.append(20000.0)
    Z2.run(J, nus, Ds, Z2.MQ, Z2.RQ)
