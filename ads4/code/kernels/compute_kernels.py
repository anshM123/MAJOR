"""
Symbolic kernels: for each conserved, Bose-symmetric <TTO> structure a (spin J, generic Delta) compute
    M^{(a)}_m(theta=0) = G(s) * P^{(a)}_m(s),   s = (Delta+J)/2,   source T(e_+) (helicity +2),
with G as in pipeline.amplitude.  Prints factorised P_m and the null tests (m = 1, 2, 3 must vanish:
int E = H and int e^{-+ i theta} E = P^x -+ i P^y).
usage: compute_kernels.py J [parity]
"""
import sys, time, pickle, os
import sympy as sp
import kernels_lib as kl, pipeline as pl
S = kl.S

def kernels(J, parity, source='plus'):
    monos, basis = kl.conserved(J, parity)
    if not basis:
        return monos, basis, None, []
    polys, K = kl.mono_amplitudes(J, parity, monos, source=source)
    res = []
    for b in basis:
        tot = kl.combine(polys, b)
        res.append({m: sp.factor(kl.poly_to_sym(v)) for m, v in pl.tcoeffs(tot, J).items()})
    return monos, basis, K, res

if __name__ == '__main__':
    J = int(sys.argv[1])
    pars = [sys.argv[2]] if len(sys.argv) > 2 else ['even', 'odd']
    for par in pars:
        t0 = time.time()
        monos, basis, K, res = kernels(J, par)
        print(f"=== J = {J}, parity {par}: {len(basis)} conserved structure(s); K = {K}  [{time.time()-t0:.0f}s]")
        for a, (b, r) in enumerate(zip(basis, res)):
            print(f"  structure {a}: coefficients")
            for m_, c in zip(monos, b):
                if c != 0:
                    print(f"     {m_}: {c}")
            print("  P_m(s):")
            for m in range(J, -J - 1, -1):
                print(f"     m={m:+d}: {r.get(m, 0)}")
            assert all(r.get(m, 0) == 0 for m in (1, 2, 3)), "Ward identity null test failed"
        fn = os.path.join(kl.CACHE, f"kern_J{J}_{par}.pkl")
        pickle.dump((monos, basis, K, res), open(fn, 'wb'))
