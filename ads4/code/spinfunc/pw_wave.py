"""
pw_wave (HWIN_THIMBLE.md): exact travelling-wave split of the Legendre polynomial, P_J = W_+ + W_-,
    W_+-(t) = (1/2) [P_J(t) -+ (2i/pi) Q_J(t)]      (Q_J: Ferrers function of the second kind, analytic on C minus (-inf,-1], [1, inf))
so that on the cut t = cos(theta), J theta >> 1:  W_+- ~ e^{+-i((J+1/2) theta - pi/4)} / sqrt(2 pi J sin theta)   [checked in _tests].
W_- = (i/pi) Q^{(3)}_J(t + i0) continued, W_+ = -(i/pi) Q^{(3)}_J(t - i0) continued: each wave has only the LOGARITHMIC singularity of
Q at t = +-1 (theta = 0, pi); the sum is entire.
Vectorised upward recurrence for P_n, Q_n (stable on/near the cut: both solutions are of comparable size there); derivatives in t
up to order K from (1 - t^2) y' = J (y_{J-1} - t y_J) and the differentiated Legendre equation
    (1 - t^2) y^{(k+2)} = 2(k+1) t y^{(k+1)} - (J(J+1) - k(k+1)) y^{(k)}.
Functions: pq(J, t) -> (P_J, Q_J, P_{J-1}, Q_{J-1});  waves(J, t) -> (W+, W-);  wave_derivs(J, t, sgn, K) -> array (K+1, ...).
sgn = +1 / -1 for W_+-, 0 for P_J itself.
"""
import numpy as np


def pq(J, t):
    t = np.asarray(t, complex)
    P0 = np.ones_like(t); Q0 = 0.5*np.log((1 + t)/(1 - t))
    if J == 0:
        return P0, Q0, None, None
    P1 = t.copy(); Q1 = t*Q0 - 1.0
    Pm, Qm, Pc, Qc = P0, Q0, P1, Q1
    for n in range(1, J):
        Pn = ((2*n + 1)*t*Pc - n*Pm)/(n + 1)
        Qn = ((2*n + 1)*t*Qc - n*Qm)/(n + 1)
        Pm, Qm, Pc, Qc = Pc, Qc, Pn, Qn
    return Pc, Qc, Pm, Qm


def _comb(P, Q, sgn):
    if sgn == 0:
        return P
    return 0.5*(P - sgn*(2j/np.pi)*Q)


def waves(J, t):
    P, Q, _, _ = pq(J, t)
    return _comb(P, Q, 1), _comb(P, Q, -1)


def wave_derivs(J, t, sgn, K):
    """y^{(k)}(t), k = 0..K, for y = W_sgn (sgn = +-1) or P_J (sgn = 0)."""
    t = np.asarray(t, complex)
    P, Q, Pm, Qm = pq(J, t)
    y = _comb(P, Q, sgn); ym = _comb(Pm, Qm, sgn)
    one = 1.0 - t*t
    out = [y]
    if K >= 1:
        out.append(J*(ym - t*y)/one)
    for k in range(0, K - 1):
        out.append((2*(k + 1)*t*out[k + 1] - (J*(J + 1) - k*(k + 1))*out[k])/one)
    return np.array(out)


def _tests():
    import mpmath as mp
    mp.mp.dps = 40
    rng = np.random.default_rng(3)
    worst = 0.0
    for J in (4, 40, 100, 300):
        for th in (0.003, 0.02, 0.3, 1.2):
            t = np.cos(th + 1e-3j*rng.normal()) + 1e-6j*rng.normal()
            P, Q, _, _ = pq(J, np.array([t]))
            Pe = complex(mp.legenp(J, 0, mp.mpc(t), type=2)); Qe = complex(mp.legenq(J, 0, mp.mpc(t), type=2))
            e = max(abs(P[0] - Pe)/abs(Pe), abs(Q[0] - Qe)/abs(Qe))
            worst = max(worst, e)
            Wp, Wm = waves(J, np.array([t]))
            ths = complex(mp.acos(mp.mpc(t)))
            asy = np.exp(1j*((J + .5)*ths - np.pi/4))/np.sqrt(2*np.pi*(J + .5)*np.sin(ths))
            d = wave_derivs(J, np.array([t]), 1, 4)[:, 0]
            h = 1e-6*(1 - t.real)
            fd = (waves(J, np.array([t + h]))[0][0] - waves(J, np.array([t - h]))[0][0])/(2*h)
            print(f"J={J:3d} th={th}: rel err P,Q vs mpmath {e:.1e};  W+/asym {Wp[0]/asy:.4f}  W-/conj-asym "
                  f"{Wm[0]/np.conj(asy) if abs(t.imag) < 1e-12 else Wm[0]*asy*2*np.pi*(J+.5)*np.sin(ths):.4f};  W+' rel err {abs(d[1]/fd - 1):.1e}")
    print('worst P/Q rel err', worst)


if __name__ == '__main__':
    _tests()
