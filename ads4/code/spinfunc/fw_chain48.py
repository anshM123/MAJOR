"""
fw_chain48 (FOLD_WINDOW_LARGEJ.md §J = 48): the ee1000 J = 48 chain (EE_DG1000_LARGEJ §3a; shape model K = 4, x in [0.05, 0.8], rows 1048,
hybrid states, sm_curvfix kappa-image subtraction, stage post) re-run on the J = 48 K2 rows RECOMPUTED with the fixed engine
(fw_rows: spin-J Fourier floor + Q_J recurrence; state + reference X35 and refb), instead of the ee1000_pts records.
Checks printed: rank-1 basis quality (sv3, imag) at the state and the reference; the smooth-offset fit of SMOOTH_BOUND §3 (sm_offset.fitr,
x in [0.1, 0.45] and [0.1, 0.7]); old-vs-new record differences.
Variants (FW48_VAR): 'k2' (all K2, new records), 'd48' (EE1000_DSTATE=1: the old ee1000_n1 direct state sub-row with the NEW K2 references
and bases).  FW48_REFBMAX (default 0.75): refb records with x_ref above this are excluded from the reference basis (x_ref = 0.8: order-0
action asymmetric at 3e-2 with the fixed engine).
Usage: FW48_VAR=k2 python3 fw_chain48.py  -> cache/fw_chain48_<var>.txt
"""
import os, sys, io, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
VAR = os.environ.get('FW48_VAR', 'k2')
for k, v in dict(EE1000_XMAX='0.8', EE1000_K='4', EE1000_MODEL='shape', EE1000_XMIN='0.05', W2_DMIN='1000', EE1000_JEVEN='48',
                 EE1000_JODD='', SM_CURV='row', OMP_NUM_THREADS='1', L1_DEG='6', W2_BASIS='ee', W2_MODES='ee').items():
    os.environ.setdefault(k, v)
if VAR == 'd48':
    os.environ['EE1000_DSTATE'] = '1'
import numpy as np
import sm_curvfix as CF
import ee1000_chain as CH
import sm_offset as SO
import l7_dc2_analysis as An
import l8_analysis as L8

W, Q, C = CH.W, CH.Q, CH.C
J0 = 48
DREF = 2e4
RBMAX = float(os.environ.get('FW48_REFBMAX', '0.75'))
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


def fwk2(kind, J):
    R = pickle.load(open(os.path.join(C, 'fw_k2_%s_J%d.pkl' % (kind, J)), 'rb'))
    mr = CH.E.mphys(DREF, J)
    drop = [k for k in R if k[0] == DREF and k[1]/mr > RBMAX + 1e-9]
    return {k: v for k, v in R.items() if k not in drop}, drop


_ee = W.ee_even
_q2 = Q._q2_even


def ee_even(J):
    return fwk2('EE', J)[0] if J == J0 else _ee(J)


def q2_even(kind, J):
    if J == J0 and kind == 'C3H':
        return fwk2('C3H', J)[0]
    return _q2(kind, J)


W.ee_even = ee_even
Q._q2_even = q2_even
if CH.DSTATE:
    W._q2_even_orig = lambda kind, J: q2_even(kind, J)


def basis_report():
    R = Q._strip(ee_even(J0), J0)
    for D in (1048.0, DREF):
        m = CH.E.mphys(D, J0)
        pts = sorted((k[1]/m, R[k]) for k in R if k[0] == D and k[1]/m >= 0.15)
        sv, vs = An.rank1([p[1]['M0'] for p in pts])
        asym = max(abs(p[1]['M0'][0, 1] - p[1]['M0'][1, 0])/np.abs(p[1]['M0']).max() for p in pts)
        pr(f"# basis Delta={D:g}: {len(pts)} records x in [{pts[0][0]:.3f}, {pts[-1][0]:.3f}]: sv3 = {sv[2]:.1e} (sv2 {sv[1]:.1e}), "
           f"imag = {max(L8.realvec(v)[1] for v in vs):.1e}, max asym {asym:.1e}")
    old = pickle.load(open(os.path.join(C, 'ee1000_pts_EE_J48_k2.pkl'), 'rb'))
    new, drop = fwk2('EE', J0)
    d = [np.abs(new[k]['M0'] - old[k]['M0']).max()/np.abs(old[k]['M0']).max() for k in new if k in old and k[0] == 1048.0]
    pr(f"# old (ee1000) vs new (fixed engine) EE records at 1048: {len(d)} common, rel. M0 difference median {np.median(d):.1e}, max {max(d):.1e}; "
       f"refb dropped (x_ref > {RBMAX}): {[round(k[1]/CH.E.mphys(DREF, J0), 3) for k in drop]}")


def offset_report():
    g = Q.even_grid(J0)
    r = g[1048.0]; xs, m = r['xs'], r['m']
    pr(f"# smooth-offset fit (SMOOTH_BOUND §3 / sm_offset.fitr), row 1048, {len(xs)} nodes, after the kappa image (rows already corrected "
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
    pr(f"# J = 48 (Delta_g = 1000): min lam/N_s = {w['lam']:+.4f} at m/Lam = {w['mfl']:.4f}; states {len(recs)}")
    pr("# " + CF.summary()[:400])
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'fw_chain48_%s.txt' % VAR), 'w').write(buf.getvalue())
