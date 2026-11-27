"""P4: size and sign of the back-to-back kernel.  In m1_remainder, Psi_nu(z) = z P_{-1/2+i nu}(X)/(1-z)^2 with
X = (1+z)/(1-z); the back-to-back point z = -1 is X = 0 (theta = pi/2), and z in (-1, 0) is X = cos(theta) in (0, 1).

Checks (60 digits; P_nu(0) from the hypergeometric form DLMF 14.3.1 with mu = 0, x = 0):
  * DLMF 14.5.1 at mu = 0:  P_{-1/2+i tau}(0) = sqrt(pi) / |Gamma(3/4 + i tau/2)|^2   (exact)
  * DLMF 14.20.15:          P_{-1/2+i tau}(cos th) ~ (th/sin th)^{1/2} I_0(tau th)  ->  e^{pi tau/2}/sqrt(2 pi tau) at th = pi/2
  * DLMF 14.20.9 (Mehler-Dirichlet): P_{-1/2+i tau}(cos th) = (2/pi) int_0^th cosh(tau phi)/sqrt(2(cos phi - cos th)) dphi > 0,
    i.e. the kernel is a POSITIVE superposition of cosh(tau phi), phi in [0, th]: no cancellation can come from the kernel.
Output: out_p4_conical.txt
"""
import mpmath as mp

mp.mp.dps = 60

if __name__ == "__main__":
    out = []
    for tau in (1, 10, 100, 300):
        tau = mp.mpf(tau)
        nu = -mp.mpf(1) / 2 + 1j * tau
        lhs = mp.hyp2f1(-nu, nu + 1, 1, mp.mpf(1) / 2, maxterms=10**6)   # Ferrers P_nu(0)
        exact = mp.sqrt(mp.pi) / abs(mp.gamma(mp.mpf(3) / 4 + 1j * tau / 2)) ** 2
        asym = mp.exp(mp.pi * tau / 2) / mp.sqrt(2 * mp.pi * tau)
        md = (2 / mp.pi) * mp.quad(lambda p: mp.cosh(tau * p) / mp.sqrt(2 * mp.cos(p)), [0, mp.pi / 4, mp.pi / 2]) \
            if tau <= 100 else mp.nan
        out.append(f"tau={int(tau):5d}  P(0)={mp.nstr(mp.re(lhs), 12):>22}  sqrt(pi)/|G|^2={mp.nstr(exact, 12):>22}  "
                   f"e^(pi tau/2)/sqrt(2 pi tau)={mp.nstr(asym, 12):>22}  Mehler-Dirichlet={mp.nstr(md, 12)}")
    txt = "\n".join(out)
    print(txt)
    open(__file__.replace("p4_conical_check.py", "out_p4_conical.txt"), "w").write(txt + "\n")
