"""mk_table: markdown version of pcyc_compare's combined table (for PCYC_CLOSE.md section 3.4)."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import io, contextlib
import numpy as np
import pcyc_compare as PC
with contextlib.redirect_stdout(io.StringIO()):
    rows = PC.main()
groups = {}
for r in rows:
    groups.setdefault((r['Delta'], r['nu'], r['cont']), {}).setdefault(r['f'], {})[r['est']] = r
print('| Δ | ν | contour | runs | `⟨Q⟩` rel. error per run | `D′` error | Re T | Im T | σ | `|T|/σ` | thr. | `(|T| + 2σ)/thr` |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for k in sorted(groups):
    runs = [(d.get('CV', d.get('plain')) if k[2] == 'n3' else d['plain']) for d in groups[k].values()]   # CVs help on n3, not on n4
    Ts = np.array([r['T'] for r in runs]); sj = np.array([r['sT'] for r in runs]); n = len(runs)
    Tm = Ts.mean()
    sc = float(np.sqrt(np.sum(np.abs(Ts - Tm)**2)/(n - 1))) if n > 1 else 0.0          # sample std of the runs
    sig = max(float(np.sqrt(np.sum(sj**2))/n), sc/np.sqrt(n) if n > 1 else 0.0, runs[0]['qerr'])
    thr = 2.1e-4*(300.0/k[0])**2
    est = '/'.join('%.0e' % float(np.sqrt(max(r['sT']**2 - r['qerr']**2, 0))) for r in runs)
    print(f'| {k[0]:g} | {k[1]:g} | {k[2]} | {n} | {est} | {runs[0]["qerr"]:.0e} | {Tm.real:+.1e} | {Tm.imag:+.1e} | {sig:.1e} | '
          f'{abs(Tm)/sig:.2f} | {thr:.1e} | {(abs(Tm) + 2*sig)/thr:.2g} |')
