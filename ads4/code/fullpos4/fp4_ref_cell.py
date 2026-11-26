"""
fullpos4/fp4_ref_cell.py (FULLPOS4.md §5.2) [ball arithmetic]: naive interval enclosure of the §78 reference over a twist cell
a in [a0, a0 + delta] at fixed integer J (all of c, d, h, w as balls; W_j and the jet propagate the cell).  Prints the output ball.
Usage: fp4_ref_cell.py J a0 b nu M prec delta1,delta2,...
"""
import sys, os, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import fp4_ref as R

J = int(sys.argv[1]); a0 = arb(sys.argv[2]); b = int(sys.argv[3]); nu = float(sys.argv[4]); M = int(sys.argv[5]); ctx.prec = int(sys.argv[6])
for dl in sys.argv[7].split(','):
    t0 = time.time()
    delta = arb(dl)
    acell = a0 + delta/2 + arb(0, 1)*0 if False else arb(a0 + delta/2, delta/2)
    P = R.params(J, acell, b)
    W, big = R.Wcoef(P['c'], P['h'], P['d'], P['s'], M)
    r = R.kref(3, nu, P, W)
    line = dict(J=J, a0=sys.argv[2], delta=dl, nu=nu, M=M, prec=ctx.prec, kref=r.str(10, radius=True), rad=float(r.rad()), secs=time.time() - t0)
    print(json.dumps(line), flush=True)
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'fp4_ref_cell.jsonl'), 'a').write(json.dumps(line) + '\n')
