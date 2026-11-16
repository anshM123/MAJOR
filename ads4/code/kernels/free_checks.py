"""
Validation of the pipeline against free fields (d = 3):
  position-space Wick contraction -> expansion in embedding structures -> light transform + Fourier transform
  -> normalised amplitudes ahat^O_{m,h}, to be compared with the Fock-space results of fock_free.py:
     ahat^T_{2,2} = 1 (Ward identity), ahat^T_{-2,2} = c = +1 (scalar), -1 (fermion),
     |ahat^{J_s}_{+-s,2}| = 1, ahat^{J_s}_{m,2} = 0 for |m| < s, |ahat^{phi^2}_{0,2}| = |ahat^{psibar psi}_{0,2}| = 1.
Run: python3 free_checks.py [maxspin]
"""
import sys, time
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
import structures as st, wick_free as wf, pipeline as pl, kernels_lib as kl
mp.mp.dps = 30

def fit(theory, J, Delta, c4, G4=None):
    sig = Fr(Delta + J, 2)
    out = {}
    for par in ('even', 'odd'):
        integer_shift = (sig - Fr(1, 2)).denominator == 1 if par == 'even' else sig.denominator == 1
        if not integer_shift:
            continue
        monos = st.independent_subset(st.monomials(J, par))[0]
        if not monos:
            continue
        if theory == 'scalar':
            corr = lambda x4, z4, x1, z1, x3, z3: wf.scalar_3pt(c4, wf.c_scalar_T(), wf.c_scalar_T(), x4, z4, x1, z1, x3, z3)
        else:
            corr = lambda x4, z4, x1, z1, x3, z3: wf.fermion_3pt(c4, G4, wf.c_fermion_T(), 'z', wf.c_fermion_T(), 'z',
                                                               x4, z4, x1, z1, x3, z3)
        if par == 'even':
            expo = (int(sig - Fr(9, 2)), int(Fr(1, 2) - sig), int(Fr(1, 2) - sig))
        else:
            expo = (int(sig) - 5, -int(sig), -int(sig))
        try:
            lam = wf.project_general(corr, monos, expo)
        except ValueError:
            continue                      # no solution: correlator has the other parity
        out[par] = (monos, lam)
    return sig, out

def amplitudes(J, sig, par, monos, lam, source, n=(1, 0)):
    polys, K = kl.mono_amplitudes(J, par, monos, source=source, sig=sig, n=n)
    tot = pl.ZERO
    for p, c in zip(polys, lam):
        c = sp.Rational(c)
        tot = tot + pl.C(int(c.p), int(c.q))*p
    tot = pl.redI(tot)
    G = kl.G_prefactor(par, K, mp.mpf(sig.numerator)/sig.denominator)
    res = {}
    for m, v in pl.tcoeffs(tot, J).items():
        val = complex(0)
        z = mp.mpc(0)
        for (et, ed, es, eI), c in v.to_dict().items():
            assert int(ed) == 0 and int(es) == 0
            z += mp.mpf(int(c.p))/int(c.q)*(1j)**int(eI)
        res[m] = G*z
    return res

def run(theory, name, J, Delta, c4, G4, nuO, calib):
    """calib = (prefactor of the 3-pt function) / sqrt(norm of O * C_T) in the units of wick_free"""
    sig, fits = fit(theory, J, Delta, c4, G4)
    for par, (monos, lam) in fits.items():
        # conserved subspace membership
        cm, cb = kl.conserved(J, par, sval=sig)
        assert cm == monos
        M = sp.Matrix(cb + [lam])
        in_span = M.rank() == len(cb)
        line = f"{theory:7s} {name:10s} parity {par:4s} s={sig}: #conserved structures at this s = {len(cb)}, free correlator in span: {in_span}"
        print(line)
        for source, h in (('plus', 2), ('minus', -2)):
            amp = amplitudes(J, sig, par, monos, lam, source)
            NT = kl.norm_T_plus()
            outs = []
            for m in range(-J, J + 1):
                a = amp.get(m, mp.mpc(0))
                NO = kl.norm_formula(Delta, J, m)
                if NO == 0:
                    outs.append(f"m={m}: M={mp.nstr(abs(a), 3)} (norm 0)")
                    continue
                ahat = 2*mp.pi*a*calib/mp.sqrt(NO*nuO*NT)
                outs.append(f"m={m}: ahat={mp.nstr(mp.chop(ahat), 12)}")
            print(f"    h={h:+d}: " + "; ".join(outs))
    return fits

if __name__ == '__main__':
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t0 = time.time()
    fourpi = 4*mp.pi
    # scalar: W = (4 pi)^{-3} sum l_a S_a, N_O = nu_O/(4 pi)^2, C_T = (3/2)/(4 pi)^2
    cal_s = lambda: fourpi**(-3)/mp.sqrt(fourpi**(-2)*mp.mpf(3)/2*fourpi**(-2))
    nuS = {0: 2, 2: mp.mpf(3)/2, 4: 645120, 6: 30656102400}
    # fermion (Dirac = 2 Majorana): kappa c_F = 1/(8 pi) up to sign (C_T^Dirac = 3/(16 pi^2), nu_T = 12)
    cal_f = lambda: -1/(8*mp.pi)/mp.sqrt(12)      # kappa c_F / sqrt(nu_T); sign fixed by ahat^T_{2,2} = +1
    nuF = {0: 2, 2: 12, 4: 1290240, 6: 61312204800}
    print("== free real scalar ==")
    run('scalar', 'T', 2, 3, wf.c_scalar_T(), None, nuS[2], cal_s())
    run('scalar', 'phi^2', 0, 1, wf.c_scalar_J(0), None, nuS[0], cal_s())
    for sp_ in range(4, smax + 1, 2):
        run('scalar', f'J_{sp_}', sp_, sp_ + 1, wf.c_scalar_J(sp_), None, nuS[sp_], cal_s())
        print(f"   [{time.time()-t0:.0f}s]")
    print("== free fermion ==")
    run('fermion', 'T', 2, 3, wf.c_fermion_T(), 'z', nuF[2], cal_f())
    run('fermion', 'psibar psi', 0, 2, {(0, 0): Fr(1)}, '1', nuF[0], cal_f())
    for sp_ in range(4, smax + 1, 2):
        run('fermion', f'J_{sp_}', sp_, sp_ + 1, wf.c_fermion_J(sp_), 'z', nuF[sp_], cal_f())
        print(f"   [{time.time()-t0:.0f}s]")
