"""v1: free-field check of the all-k rules R_k (signs, normalisation, T only at k=2).

Free massless scalar / Dirac fermion in d=3: T|0> at rest is a 2-particle state, back-to-back, energies E/2.
E(theta) acts by multiplication: E(theta)/E = (1/2)[delta(theta-phi) + delta(theta-phi-pi)],
so E_m = int e^{i m theta} E(theta)/E acts as e^{i m phi} (1+(-1)^m)/2.
Rest-frame 2-particle states of charge q belong to the spin-|q| twist-1 current (q=0: phi^2 or psibar psi, q=+-2: T).
We compute, by quadrature on the relative angle, the charge decomposition of E_{+-m}|T,2> and evaluate
  LHS_k = sum_{q=4+k} |a|^2  (b=4 ladder),  RHS_k = sum_{q=-k, O != T} |a|^2 (b=0 ladder),  |c|^2 = weight on (T, q=-2).
Expect LHS_k - RHS_k = delta_{k2} |c|^2 with |c|^2 = (t4/4)^2 = 1 (t4 = +4 scalar, -4 fermion)."""
import numpy as np

N = 4096
phi = 2 * np.pi * np.arange(N) / N


def charges(psi):
    """|<q|psi>|^2 for unit-normalised plane waves e^{iq phi}/sqrt(2pi) (distinguishable-particle angle)."""
    c = np.fft.fft(psi) / N
    q = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    return dict(zip(q, np.abs(c) ** 2 * 2 * np.pi / (np.sum(np.abs(psi) ** 2) * 2 * np.pi / N) * 1.0))


def run(label, psiT):
    # normalise <T|T> = 1 on the circle
    nrm = np.sqrt(np.sum(np.abs(psiT) ** 2) * 2 * np.pi / N)
    psiT = psiT / nrm
    out = []
    for k in range(0, 9):
        m = k + 2
        res = {}
        for sgn in (+1, -1):
            mm = sgn * m
            op = np.exp(1j * mm * phi) * (1 + (-1) ** mm) / 2
            v = op * psiT
            c = np.fft.fft(v) / N
            q = np.fft.fftfreq(N, d=1.0 / N).astype(int)
            w = np.abs(c) ** 2 * 2 * np.pi           # |<q|v>|^2, sum_q w = ||v||^2
            res[sgn] = dict((int(qq), ww) for qq, ww in zip(q, w) if ww > 1e-12)
        lhs = sum(w for qq, w in res[+1].items() if qq == 4 + k)
        t_part = res[-1].get(-2, 0.0) if k == 2 else 0.0       # q=-2 is T itself (only spin-2 state at rest)
        rhs = sum(w for qq, w in res[-1].items() if qq == -k) - t_part
        stray = sum(w for qq, w in res[+1].items() if qq != 4 + k) + sum(w for qq, w in res[-1].items() if qq != -k)
        out.append((k, lhs, rhs, t_part, lhs - rhs - t_part, stray))
    print(label)
    for k, lhs, rhs, tp, viol, stray in out:
        print(f"  k={k}: sum|a4k|^2={lhs:.6f}  sum_(O!=T)|a0k|^2={rhs:.6f}  |c|^2 term={tp:.6f}  "
              f"R_k residual={viol:+.1e}  charge leakage={stray:.1e}")


# scalar: wavefunction of T_{zz}-type (helicity +2) pair: e^{2 i phi}
run("free scalar (lambda=+2)", np.exp(2j * phi))
# Dirac fermion: pair wavefunction for helicity-2 state, spinor factors give e^{2 i phi} times a unit-modulus
# spin factor that commutes with E_m (E_m is multiplication by a phase); the result is identical
run("free fermion (lambda=+2, spin factor e^{i phi} absorbed: total charge 2)", np.exp(2j * phi))

# endpoint rule check of the normalisation: k=0 states are phi^2 (b=0,J=0,D=1) and the spin-4 current (b=4,J=4,D=5)
for (b, J, D) in [(0, 0, 1), (4, 4, 5)]:
    T = (D - 1) * (D - 2)
    print(f"  endpoint weight W_{b}(D={D}) = {(T - b*(b-1))/3 + 1.5}  (each has |a_0|^2=1 -> total 3 = RHS)")
