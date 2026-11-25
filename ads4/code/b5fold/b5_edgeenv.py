"""
b5_edgeenv (BATCH5_FOLD.md §2.3): DIAGNOSTIC (not used for any grade).  The even-J twist edges J >= 22 are blocked by ng_bind's edge
term: |Gamma| <= G_edge = max(1e4, 10 max|Gamma| of the last 4 own-row nodes), flat, on x > x_last, for every entry of the 2x2
(++, +-, and the off-diagonal c_off, which is measured at <~ 10 c2).  The odd-J protocol (hs_odd E_miss, accepted for J = 21-31) instead
covers the unmeasured edge with the light-cone envelope 3 A_e ((1 - x_l)/(1 - x))^2, A_e = max |Gamma_e| of the last three valid nodes.
Here: the same sign-free errmat structure as ng_bind (|weights| of the 256-node hybrid, Frobenius norm of the normalised 2x2), with
the per-entry odd-J envelope on x > x_last, to show what the even edge would cost under the odd-J rule.
Usage: python3 b5_edgeenv.py J,J,...   -> cache/b5_edgeenv.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_core as C
import numpy as np
import ng_bind as NB
H, V = C.H, C.V
NC = H.NC
C.SET[0] = 'fix'
lines, res = [], {}


def pr(s=""):
    print(s, flush=True); lines.append(s)


if __name__ == '__main__':
    Js = [int(v) for v in sys.argv[1].split(',')]
    for J in Js:
        D = 1000.0 + J
        st = json.load(open(os.path.join(C.B5C, 'state', 'b5_state_J%d.json' % J)))
        bad = set(st['states'][0]['failed']) | set(st['states'][0]['rz'])
        _rec0 = H.HN.records
        H.HN.records = lambda kind, J_, D_, jk=None, **kw: {k: v for k, v in _rec0(kind, J_, D_, jk, **kw).items()
                                                             if round(k[1]/H.HN.mphys(D_, J_), 6) not in bad}
        V.INJ[J] = [D]; NB.NU_LOW = 50.0
        g2, _ = NB.grid_for(J, 2, False)
        own = g2[D]; xs = own['xs']; xl = xs[-1]
        A = {}
        for h in ('++', '+-', 'off'):
            sc = own['m'] if h == 'off' else 1.0
            A[h] = max(float(np.max(np.abs(np.asarray(own[(kd, h)])[-3:])))*sc for kd in ('C2', 'C3H'))
        for s in st['states']:
            mfl = s['mfl']
            d2 = NC.dens(1000.0, J, mfl, NB.model_of(J, g2), 256)
            d2['abs'] = NB.absdens(1000.0, J, mfl, 256)
            x = d2['x']; sel = x > xl
            env = lambda h: 3*A[h]*np.maximum(1.0, ((1 - xl)/(1 - x))**2)
            n = len(d2['hs']); M = np.zeros((n, n)); Ns = d2['Ns']
            for i, h in enumerate(d2['hs']):
                a2, a3 = d2['abs'][h]
                M[i, i] = ((a2 + a3)*env(h))[sel].sum()/Ns[i]
            a2, a3 = d2['abs']['off']
            M[0, 1] = M[1, 0] = ((a2 + a3)*env('off'))[sel].sum()/np.sqrt(Ns[0]*Ns[1])
            Ee = float(np.sqrt((M*M).sum()))
            Eflat = s['direct']['E_edge']
            res[f'{J}:{mfl:.7f}'] = dict(J=J, mfl=mfl, xl=float(xl), A=A, E_env=Ee, E_flat=Eflat)
            pr(f"J={J} m/Lam={mfl:.7f}: x_last {xl:.4f}, A(++, +-, off) = {A['++']:.0f}, {A['+-']:.0f}, {A['off']:.1f} c2 | edge: ng_bind flat "
               f"{Eflat:.3f} -> odd-J envelope {Ee:.3f}")
        V.INJ.pop(J, None); H.HN.records = _rec0
    open(os.path.join(C.B5C, 'b5_edgeenv.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.B5C, 'b5_edgeenv.json'), 'w'), indent=1)
