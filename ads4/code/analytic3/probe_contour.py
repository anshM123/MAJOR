import sys, math, time
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import rt_cert as RC
La, Lb, k0 = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
d = math.exp(-La); F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
for (name, p, Q, Qp, dress) in (('psi2', 3, RC.Q2, RC.Q2p, True), ('psi3', 4, RC.Q3, RC.Q3p, False)):
    fac = d**3*(2/RC.LAMf)*2 if name == 'psi2' else d**3*(d*d/RC.TAU_COEF)*RC.LAMup*2
    for nm, fn in (('bottom', lambda: RC.piece_bottom(La, p, Q, Qp, dress)), ('left', lambda: RC.piece_left(La, p, Q, dress, (p-0.5)-0.25+k0)),
                   ('top', lambda: RC.piece_top(La, p, Q, dress, k0)[0]), ('right', lambda: RC.piece_right(La, p, Q, dress, (p-0.5)-0.25+k0))):
        t0 = time.time()
        try:
            v = fac*fn()/F; print(name, nm, f'{v:.3e}', f'{time.time()-t0:.0f}s', flush=True)
        except AssertionError as e:
            print(name, nm, 'ASSERT', e, flush=True)
