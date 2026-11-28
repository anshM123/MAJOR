#!/usr/bin/env python3
"""Moduli-space criterion  C_T(origin) >= (7/8) C_T(free IR theory on moduli space)  on known SCFTs,
and the F-theorem number for the 28-chiral candidate.  Units: C_T / C_T(real scalar) unless stated."""
import numpy as np
from scipy.integrate import quad
from scipy.special import polygamma

# ---- N=2 chiral of R-charge D on squashed sphere: F_b = -log s_b(iQ(1-D)/2)
#      log s_b(iy) = int_0^inf dw/w [ sinh(2yw)/(2 sinh(bw) sinh(w/b)) - y/w ]
def F_chiral(D, b):
    Q = b + 1 / b
    y = Q * (1 - D) / 2
    def integrand(w):
        if w < 1e-4:   # series to avoid cancellation
            # expand to O(w): sinh(2yw)/(2 sinh(bw)sinh(w/b)) - y/w
            return y * (4 * y * y - b * b - 1 / b ** 2) / 6
        r = np.exp((2 * y - Q) * w) * (-np.expm1(-4 * y * w)) / (np.expm1(-2 * b * w) * np.expm1(-2 * w / b))
        return (r - y / w) / w
    v, _ = quad(integrand, 0, np.inf, limit=500, epsabs=1e-13)
    return -v


def CT_chiral(D, h=1e-3):
    """C_T in units of a real scalar: 1910.03600 eq. C_T=(48/pi^2) d^2F/db^2 with C_T(scalar)=3/2."""
    d2 = (F_chiral(D, 1 + h) - 2 * F_chiral(D, 1) + F_chiral(D, 1 - h)) / h ** 2
    return (48 / np.pi ** 2) * d2 / 1.5


if __name__ == "__main__":
    print("F(free chiral) =", F_chiral(0.5, 1.0), " expect log2/2 =", np.log(2) / 2)
    print("C_T(free chiral) =", CT_chiral(0.5), " expect 4")
    cxyz = 3 * CT_chiral(2 / 3)
    print("C_T(XYZ) =", cxyz, " (= 3 chirals at R=2/3);  ratio to free chiral per field", CT_chiral(2 / 3) / 4)

    rows = []
    # (name, C_T(origin) in scalar units, n_free_chirals at generic moduli point)
    rows.append(("XYZ / N=2 SQED_1 (branch: 1 chiral)", cxyz, 1))
    for n in [2, 3, 4, 5, 10, 50]:
        ct = 12 * (2 * n * n * polygamma(1, n / 2 + 1) + (np.pi ** 2 * n + 4) * n + 4) / (np.pi ** 2 * (n + 1)) / 1.5
        rows.append((f"N=4 SQED_{n} (Higgs: {n-1} hypers)", ct, 2 * (n - 1)))
    rows.append(("T[SU(3)] (Higgs: 3 hypers)", 75.5329 / 1.5, 6))
    rows.append(("min SU(3) SCFT (54-192/pi^2)", (54 - 192 / np.pi ** 2) / 1.5, 4))
    rows.append(("min SO(5) SCFT (h=3, 2 hypers)", (54 - 64 / np.pi ** 2) / 1.5, 4))
    rows.append(("min SO(8) SCFT (h=6, 5 hypers)", (736 / 7 - 200 / np.pi ** 2) / 1.5, 10))
    rows.append(("T_3 / E6 min (11 hypers)", 160.246 / 1.5, 22))
    # N=8 (c_T already in scalar units); moduli: Sym^N(C^4/Z_k) -> N free N=8 multiplets = 4N chirals
    for name, ct, N in [("U(1)_k x U(1)_-k ABJM", 16.0, 1), ("U(2)_2 x U(1)_-2 ABJ", 64 / 3, 1),
                        ("U(2)_1 x U(2)_-1 ABJM", 37.3333, 2), ("U(2)_2 x U(2)_-2 ABJM", 42.6667, 2),
                        ("SU(2)_3 x SU(2)_-3 BLG", 46.9998, 2), ("SU(2)_4 x SU(2)_-4 BLG", 50.3575, 2),
                        ("SU(2)_5 x SU(2)_-5 BLG", 52.9354, 2)]:
        rows.append((name, ct, 4 * N))
    rows.append(("CANDIDATE (C_T in (92,93))", 92.5, 28))

    print(f"\n{'theory':42s} {'C_T(0)':>9s} {'C_T(IR)':>8s} {'ratio':>7s} {'R_lower=3.5n/C_T':>17s}")
    for name, ct, nch in rows:
        ctir = 4 * nch
        print(f"{name:42s} {ct:9.3f} {ctir:8.1f} {ct/ctir:7.3f} {3.5*nch/ct:17.4f}")
    print("\nconjecture requires ratio >= 7/8 = 0.875, i.e. R_lower <= 1")
    print("candidate: R >= 98/93 =", 98 / 93, " to 98/92 =", 98 / 92)
    Fs = np.log(2) / 8 - 3 * 1.2020569031595942 / (16 * np.pi ** 2)
    Fd = np.log(2) / 4 + 3 * 1.2020569031595942 / (8 * np.pi ** 2)
    print("F_scalar, F_Dirac, F_chiral =", Fs, Fd, 2 * Fs + Fd, " 28 chirals:", 28 * (2 * Fs + Fd))
