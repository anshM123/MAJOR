"""
fw_chain48 (FOLD_WINDOW_LARGEJ.md §J = 48): the ee1000 J = 48 chain (EE_DG1000_LARGEJ §3a; shape model K = 4, x in [0.05, 0.8], rows 1048,
hybrid states, sm_curvfix kappa-image subtraction, stage post) re-run on the J = 48 K2 rows RECOMPUTED with the fixed engine
(fw_rows: spin-J Fourier floor + Q_J recurrence; state + reference X35 and refb), instead of the ee1000_pts records.
Checks printed: rank-1 basis quality (sv3, imag) at the state and the reference; the smooth-offset fit of SMOOTH_BOUND §3 (sm_offset.fitr,
x in [0.1, 0.45] and [0.1, 0.7]); old-vs-new record differences.
Records are loaded through fw_gamma (fw_fix basis consistency; generalised helicity basis gbasis, state x in [FW_XB0, FW_XB1] =
[0.15, 0.6] here, reference x_ref in [0.3, 0.6]).  Variant FW48_VAR = 'k2' (all K2, new records); the old direct sub-row (ee1000_n1, old
engine and old structure basis) is not used.
Usage: [FW_J=48|64] FW48_VAR=k2 python3 fw_chain48.py  -> cache/fw_chain<J>_<var>.txt
"""
import os, sys, io, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
VAR = os.environ.get('FW48_VAR', 'k2')
os.environ['EE1000_JEVEN'] = os.environ.get('FW_J', '48')
for k, v in dict(EE1000_XMAX='0.8', EE1000_K='4', EE1000_MODEL='shape', EE1000_XMIN='0.05', W2_DMIN='1000', EE1000_JEVEN='48',
                 EE1000_JODD='', SM_CURV='row', OMP_NUM_THREADS='1', L1_DEG='6', W2_BASIS='ee', W2_MODES='ee').items():
    os.environ.setdefault(k, v)
os.environ.setdefault('FW_XB1', '0.6')
import numpy as np
import sm_curvfix as CF
import ee1000_chain as CH
import sm_offset as SO
import l7_dc2_analysis as An
import l8_analysis as L8

import fw_gamma as FG
os.environ['W2_DESPIKE'] = '1'                 # fw_gamma switches the despike off; the ee1000 chain keeps it
W, Q, C = CH.W, CH.Q, CH.C
J0 = int(os.environ.get('FW_J', '48'))          # FW_J = 64: the (1064, 64) row (same pipeline)
DROW = 1000.0 + J0
os.environ['EE1000_JEVEN'] = str(J0)
FG.CUR[J0] = dict(jk=None)                     # J = 48 records from fw_gamma's loaders (fw_k2 + fw_fix), generalised basis
DREF = 2e4
RBMAX = float(os.environ.get('FW48_REFBMAX', '0.75'))
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


def basis_report():
    R = Q._strip(FG.ee_even(J0), J0)
    for D in (DROW, DREF):
        FG.gbasis(J0, D, R)
    pr("# generalised basis: " + "; ".join(f"Delta={k[1]:g}: " + ", ".join(f"{kk}={vv:.1e}" if isinstance(vv, float) else f"{kk}={vv}"
                                                                         for kk, vv in v.items()) for k, v in FG.GB_LOG.items()))
    pr(f"# fw_fix conversions (records whose stored b differs from the current fw shadow): {len(FG.FIXLOG)}")


def offset_report():
    g = Q.even_grid(J0)
    r = g[DROW]; xs, m = r['xs'], r['m']
    pr(f"# smooth-offset fit (SMOOTH_BOUND §3 / sm_offset.fitr), row {DROW:g}, {len(xs)} nodes, after the kappa image (rows already corrected "
       f"by sm_curvfix: SM_CURV={CF.MODE})")
    for kd in ('C2', 'C3H'):
        for h in ('++', '+-'):
            y = np.asarray(r[(kd, h)], float)
            if CF.MODE == 'none':
                y = y - CF.curv_image(kd, h, J0, m, xs, Q.A.env, Q.A.SIG)
            s = []
            for xh in (0.45, 0.7):
                sel = (xs >= 0.1 - 1e-9) & (xs <= xh + 1e-9)
                rr, er, res, rms = SO.fitr(xs, y, h, J0, sel)
                s.append(f"x in [0.1,{xh}] (n={sel.sum()}): r = {rr:+6.0f} +- {er:4.0f} (resid {res:5.0f} of {rms:5.0f})")
            pr(f"   {kd:3s} {h}: " + "   ".join(s))
    return g


if __name__ == '__main__':
    t0 = time.time()
    pr(f"######## fw_chain48 variant {VAR}: J = 48 rows recomputed with the fixed engine (fw_k2_*_J48); SM_CURV={CF.MODE}")
    basis_report()
    offset_report()
    recs, info = CH.run(pr)
    w = min(recs, key=lambda r: r['lam'])
    pr(f"# J = {J0} (Delta_g = 1000): min lam/N_s = {w['lam']:+.4f} at m/Lam = {w['mfl']:.4f}; states {len(recs)}")
    pr("# " + CF.summary()[:400])
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'fw_chain%d_%s.txt' % (J0, VAR)), 'w').write(buf.getvalue())
