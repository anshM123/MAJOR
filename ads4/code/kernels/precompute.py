"""precompute conserved structures and per-monomial amplitudes (cached).  usage: precompute.py J parity [source]"""
import sys, time
import kernels_lib as kl
J = int(sys.argv[1]); par = sys.argv[2]; src = sys.argv[3] if len(sys.argv) > 3 else 'plus'
t0 = time.time()
monos, basis = kl.conserved(J, par)
print(f"J={J} {par}: {len(monos)} independent monomials, {len(basis)} conserved structures  [{time.time()-t0:.0f}s]", flush=True)
polys, K = kl.mono_amplitudes(J, par, monos, source=src)
print(f"amplitudes done, K={K}  [{time.time()-t0:.0f}s]", flush=True)
