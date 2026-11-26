"""rt2_contour.py La Lb: the four Regge-contour pieces (psi2 dressed, psi3) of Route T' (rt_cert.piece_*, unchanged), as ratios to F."""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import rt_cert as RC
La, Lb = float(sys.argv[1]), float(sys.argv[2]); k0 = int(sys.argv[3]) if len(sys.argv) > 3 else 2
d = math.exp(-La); t0 = time.time()
res = {}
for (name, p, Q, Qp, dress) in (('psi2', 3, RC.Q2, RC.Q2p, True), ('psi3', 4, RC.Q3, RC.Q3p, False)):
    bot = RC.piece_bottom(La, p, Q, Qp, dress)
    top, Y = RC.piece_top(La, p, Q, dress, k0)
    left = RC.piece_left(La, p, Q, dress, Y)
    right = RC.piece_right(La, p, Q, dress, Y)
    fac = d**3*(2/RC.LAMf)*2 if name == 'psi2' else d**3*(d*d/RC.TAU_COEF)*RC.LAMup*2
    res[name] = dict(bottom=fac*bot, left=fac*left, top=fac*top, right=fac*right, Y=Y)
    print(name, res[name], f'[{time.time()-t0:.0f}s]', flush=True)
F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
res['F'] = F
json.dump(res, open(os.path.join(common2.HERE, 'cache', f'contour_{La}_{Lb}.json'), 'w'))
print('F', F, 'E2/F', sum(v for k, v in res['psi2'].items() if k != 'Y')/F, 'E3/F', sum(v for k, v in res['psi3'].items() if k != 'Y')/F)
