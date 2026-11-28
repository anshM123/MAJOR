"""Shared kinematics for the independent check of the all-k angular rules.

Conventions (re-derived, see README.md):
  state |T,lambda=+2> at rest, E_m = int dphi e^{i m phi} E(phi) / E_total.
  O at rest with SO(2) charge q.  <O,q|E_m|T,2> != 0 only for q = 2+m.
  m = +(k+2) -> q = 4+k  : ladder b=4, index k, s=|q|=4+k  (needs J >= 4+k)
  m = -(k+2) -> q = -k   : ladder b=0, index k, s=|q|=k    (needs J >= k)
  |a_{k+1}|^2/|a_k|^2 = ((k+4)/(k+1))^2 (C - s(s+1))/(T_D - s(s+1)),  s=b+k,
  C = J(J+1), T_D = (D-1)(D-2)   [paper, supplement.tex Sec. sec:angular item (b)].
  endpoint weight W_b = (T_D - b(b-1))/3 + 3/2, budget sum W |a_0|^2 = 3 (D=0).
  rule R_k ([E_{k+2},E_{-k-2}]=0):  sum_O |a_{O,4,k}|^2 - sum_{O!=T} |a_{O,0,k}|^2 = delta_{k2} |c|^2,
  |c|^2 = |<T,-2|E_{-4}|T,2>|^2 = tau^2/16.
"""
from fractions import Fraction as Fr
import numpy as np


def ladder_exact(b, J, D, kmax):
    """|a_k|^2/|a_0|^2 for k=0..kmax as Fractions (D rational)."""
    D = Fr(D)
    C = J * (J + 1)
    T = (D - 1) * (D - 2)
    r = [Fr(1)]
    for k in range(kmax):
        s = b + k
        if k + 1 > J - b:        # component |q| = s+1 does not exist
            r.append(Fr(0))
            continue
        r.append(r[-1] * Fr((k + 4) ** 2, (k + 1) ** 2) * (C - s * (s + 1)) / (T - s * (s + 1)))
    return r


def ladder_float(b, J, D, kmax):
    C = J * (J + 1.0)
    T = (D - 1.0) * (D - 2.0)
    r = np.zeros(kmax + 1)
    r[0] = 1.0
    for k in range(kmax):
        s = b + k
        if k + 1 > J - b:
            break
        r[k + 1] = r[k] * ((k + 4.0) / (k + 1.0)) ** 2 * (C - s * (s + 1)) / (T - s * (s + 1))
    return r


def W(b, D):
    T = (D - 1) * (D - 2)
    return (T - b * (b - 1)) / 3 + Fr(3, 2) if isinstance(D, Fr) else (T - b * (b - 1)) / 3.0 + 1.5


def allowed(b, J, odd_b0=False):
    if b == 4:
        return J >= 4
    return (J % 2 == 0) or odd_b0
