r"""
ee5_sd_run (EESTEP_CLOSE.md §5): F_NORM's "deciding test" for (P) -- the engine's ++ J = 0 derivative-piece correction s_d for heavy
states of growing mass at FIXED functional.

The chain (repro_engine/sdchain.py, imported read-only) is: rho(nu) = I_state(D)/I_ref(D_ref = 2e4) from the direct 7-d quadrature of the
exact EE k = 4 integrand (repro_engine/direct_rho.py, imported read-only), Gamma_K = m^2 (c T_ref rho - T_st) - ..., and
s_d = 0.04 sum chi_d Re Gamma_K with the single (Delta_g-independent) Stage 6 H0 weight chi_d.  Nothing in it depends on Delta_g: s_d(D) is the
c2-unit correction of a J = 0 state of dimension D under the fixed functional.  Cached runs exist for D = 300, 1000 (and the reference);
this script adds D = 3000 (m = 10.0 Lambda_300) [and optionally 6000] on the same grid, writing ONLY to eestep/cache (direct_rho's save
directory is redirected by setting its module attribute HERE).
Usage:  python3 ee5_sd_run.py 3000 8,8,18,18,6,20,20 [numax]
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
RE = os.path.join(HERE, '..', 'repro_engine')
sys.path.insert(0, RE)
import direct_rho as DR

if __name__ == '__main__':
    D = float(sys.argv[1]); ns = [int(v) for v in sys.argv[2].split(',')]
    numax = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
    DR.HERE = HERE                     # redirect the .npz into eestep/cache (read-only use of direct_rho)
    t0 = time.time()
    fn = DR.run(D, ns, 'f64', numax)
    print(f'ee5_sd_run D={D:g} n={ns}: {time.time() - t0:.0f}s -> {fn}', flush=True)
