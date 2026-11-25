"""Symbolic-s detector amplitudes (zeroth and first moments, source helicity +2 and -2, all components m) of a canonical
<TTO> structure at fixed integer J, from the ANGULAR_PROOFS symbolic-J vectors and the read-only pipeline.
Result: dict  (moment, source) -> (K, {m: P_m(s) sympy}),  M_m = G_K(s) P_m(s),
   G_K = -4 pi i e^{i pi e14} pi^{5/2} 2^{4-2 e34} / (Gamma(e14+Ka) Gamma(e13+Kb) Gamma(e34+KB) Gamma(e34-1/2+Kg)),
   e13 = 5 - s (+1/2 odd), e14 = e34 = s (+1/2 odd).   Cached in oddz/cache/amps_<name>_J<J>.pkl."""
import os, pickle, time
from common import *

def low_spin(J, parity):
    """J = 0, 2: conserved structure solved by the read-only kernels code (symbolic s)."""
    monos, basis = kl.conserved(J, parity)
    assert len(basis) == 1
    return monos, basis[0]

def amps(name, J, moments=('zero', 'mid', 'orig'), sources=('plus', 'minus')):
    fn = guard(os.path.join(CACHE, f'amps_{name}_J{J}.pkl'))
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    parity = 'even' if name in ('pp', 'pm', 'even') else 'odd'
    if J in (0, 2):
        monos, coeffs = low_spin(J, parity)
    else:
        monos, coeffs = structure_at_J(name, J)
    out = {'monos': monos, 'coeffs': coeffs}
    for src in sources:
        kin = pl.Kin(src)
        e0 = [pl.scri_mono(m, kin) for m in monos]
        ex = {'zero': e0}
        if 'mid' in moments:
            ex['mid'] = [moment_expr(e, parity, 'mid') for e in e0]
        if 'orig' in moments:
            ex['orig'] = [moment_expr(e, parity, 'orig') for e in e0]
        for mo in moments:
            K = pl.Kbounds(ex[mo], parity)
            polys = [pl.amplitude(e, kin, parity, K=K)[0] for e in ex[mo]]
            tot = kl.combine(polys, coeffs)
            P = {m: sp.expand(kl.poly_to_sym(v)) for m, v in pl.tcoeffs(tot, J).items()}
            out[(mo, src)] = (K, P)
    pickle.dump(out, open(fn, 'wb'))
    return out

if __name__ == '__main__':
    import sys
    name, J = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    r = amps(name, J)
    print(name, J, 'done', f'{time.time()-t0:.0f}s', {k: (v[0], sorted(v[1])) for k, v in r.items() if isinstance(k, tuple)})
