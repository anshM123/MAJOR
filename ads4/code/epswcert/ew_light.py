"""
epswcert/ew_light.py (EPSW_CERT.md §5): light-side constants (K = K0 + dK0, Einstein value E = 16 Lam V3, moments M_j, tau_0) of
Omega^EE's k = 2 piece with the PUBLISHED and the REPAIRED continuum weight, in ball arithmetic.

Uses ads4/code/kcert/kc_dk0.py and kc_common.py READ-ONLY (imported, no bytecode, nothing of theirs written: kc_dk0.main() only writes
when run as a script).  The repaired weight is substituted by rebinding kc_dk0.Chi to a subclass whose constructor reads
cache/k2_continuum_D_<Dg>_rep.json instead of prl/data/k2/k2_continuum_D_<Dg>.json (same format; only x_re[0] differs).
tau_0(K) = (240/Lam^4) sqrt(E/K) (kcert/kc_bounds.py); H17 is not recomputed (its sensitivity is through M_2, M_1, which move by
< 1e-12 relative, printed below).
Usage: python3 ew_light.py   -> cache/ew_light.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kcert'))
import kc_common, kc_dk0
from flint import arb, ctx

LINES = []


def pr(s=''):
    print(s, flush=True); LINES.append(s)


class ChiRep(kc_common.Chi):
    def __init__(self, Dg, part='re'):
        import numpy as np
        d = json.load(open(os.path.join(HERE, 'cache', 'k2_continuum_D_%d_rep.json' % Dg)))
        assert d['Delta_g'] == Dg
        cs = np.array(d['centres_sigma'], float)
        self.c = cs[:, 0]; self.s = cs[:, 1]
        self.x = np.array([float(v) for v in d['x_' + part]])
        self.y = [float(v) for v in d['y_' + part]]
        self.xmax = float(np.abs(self.x).max())
        self.d = d
        self._prec = None


def main():
    pr(__doc__.split('Usage')[0])
    PUB = kc_dk0.Chi
    out = {}
    for Dg in (300, 1000):
        rows = {}
        for wt in ('pub', 'rep'):
            kc_dk0.Chi = PUB if wt == 'pub' else ChiRep
            _, r = kc_dk0.main(Dg)
            ctx.prec = 128
            L = r['Lam']; E = 16*L*r['V3']; Kc = r['K0'] + r['dK0']
            tau = lambda K: 240/L**4*(E/K).sqrt()
            rows[wt] = dict(dK0=r['dK0'], V3=r['V3'], E=E, Kc=Kc, tau1=tau(arb(1)), tauc=tau(Kc), M=[r['M%d' % j] for j in range(4)])
        kc_dk0.Chi = PUB
        p, q = rows['pub'], rows['rep']
        pr(f'\n=== Dg = {Dg} ===')
        for k in ('dK0', 'V3', 'E', 'Kc', 'tau1', 'tauc'):
            rel = (q[k] - p[k])/p[k]
            pr(f'  {k:5s}: published {p[k].str(12)}   repaired {q[k].str(12)}   relative change {rel.str(3)}')
        for j in range(4):
            pr(f'  M_{j}  : published {p["M"][j].str(12)}   repaired {q["M"][j].str(12)}   rel. change {((q["M"][j] - p["M"][j])/p["M"][j]).str(3)}')
        out[Dg] = {wt: {k: (v.str(15) if not isinstance(v, list) else [x.str(15) for x in v]) for k, v in rows[wt].items()} for wt in rows}
    json.dump(out, open(os.path.join(HERE, 'cache', 'ew_light.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'cache', 'ew_light.txt'), 'w').write('\n'.join(LINES) + '\n')


if __name__ == '__main__':
    main()
