"""
fullpos/fp_kernel.py (FULLPOS_CERT.md §4.3): the finite-mass angular action as a frequency KERNEL, [DV] (mpmath, high precision).
  H_ang,b = int g(nu) Kfin(nu) dnu + (k=3 part),
  Kfin(nu; J, T, b) = (1-r)^3 sum_{k<=J-b} D_k^2 U_{b,k} [ z_{k+2}(nu)(1/(k+3) - r/(k+1)) + (2 y_{k+2}(nu) + V_{b,k} z_{k+2}(nu))/T ],
  KLR(nu; r)        = -(4x+9)(4x+49)/(4608(4x+1)) [R_-(nu) P_{-1/2+i(nu+i)}(X) + R_+(nu) P_{-1/2+i(nu-i)}(X)],  X = (1+r)/(1-r)
(identity I2 of fp_norm.py: KLR is the T -> oo limit of Kfin at fixed r).  The individual terms are ~ e^{(pi/2) nu}; the kernel is O(1).
Also: the §31 identity 2 Y_k = ((k+3)^2/2) L_k + ((k+4)^2/2) L_{k+1}, checked exactly at the polynomial level (fmpq) for k <= 60.
"""
import sys, os, json
sys.dont_write_bytecode = True
import mpmath as mp
import fp_common as FC


def kfin(nu, J, T, b, dps):
    mp.mp.dps = dps
    nu = mp.mpf(nu); s = nu*nu + mp.mpf(1)/4
    n = J - b; C = J*(J + 1); T = mp.mpf(T); r = C/T
    f = [mp.mpf(0), mp.mpf(1)]
    for m in range(1, n + 6):
        f.append((-(s - 2*m*m)*f[m] - (m - 1)**2*f[m - 1])/(m + 1)**2)
    U = mp.mpf(1); tot = mp.mpf(0); big = mp.mpf(0)
    for k in range(n + 1):
        m = k + 2
        z = 12*f[m] - (m - 1)*f[m - 1] + (m + 1)*f[m + 1]
        y = (mp.mpf(12 - m*m)/2*f[m] + mp.mpf(13 - 7*m)/2*f[m - 1] + mp.mpf(13 + 7*m)/2*f[m + 1]
             + mp.mpf((m - 2)**2)/4*f[m - 2] + mp.mpf((m + 2)**2)/4*f[m + 2])
        V = mp.mpf(4*k)/(k + 3) if b == 0 else 4 + mp.mpf(12)/(k + 1)
        t = mp.binomial(k + 3, 3)**2*U*(z*(mp.mpf(1)/(k + 3) - r/(k + 1)) + (2*y + V*z)/T)
        tot += t; big = max(big, abs(t))
        u = (b + k)*(b + k + 1)
        U = U*(C - u)/(T - u)
    return (1 - r)**3*tot, (1 - r)**3*big, r


def klr(nu, r, dps=40):
    mp.mp.dps = dps
    kj = json.load(open(os.path.join(FC.DATA, 'kernels.json')))['D_k2_EE']
    nu = mp.mpf(nu); x = nu*nu; r = mp.mpf(r); X = (1 + r)/(1 - r)
    ev = lambda e: mp.mpmathify(eval(e, {'nu': nu, 'I': mp.mpc(0, 1)}))
    Pm = mp.legenp(-mp.mpf(1)/2 + 1j*(nu + 1j), 0, X, type=3)
    Pp = mp.legenp(-mp.mpf(1)/2 + 1j*(nu - 1j), 0, X, type=3)
    return mp.re(-(4*x + 9)*(4*x + 49)/(4608*(4*x + 1))*(ev(kj['R_minus'])*Pm + ev(kj['R_plus'])*Pp))


def check_31(K=60):
    from flint import fmpq
    import fp_moments as FM
    Z, Y, f = FM.zy_fmpq(K + 1)
    L = [Z[0]*fmpq(1, 3)] + [(Z[k] - fmpq(k, k + 3)*Z[k - 1])*fmpq(1, k + 3) for k in range(1, K + 2)]
    for k in range(K + 1):
        assert 2*Y[k] == fmpq((k + 3)**2, 2)*L[k] + fmpq((k + 4)**2, 2)*L[k + 1], k
    return True


def main():
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr(f'§31 identity 2Y_k = ((k+3)^2/2)L_k + ((k+4)^2/2)L_(k+1), polynomial level, k <= 60: {check_31()}')
    for (J, tw) in ((40, 1000), (200, 1000), (600, 1000), (1000, 1000), (2000, 1000)):
        Delta = J + tw; T = (Delta - 1)*(Delta - 2)
        for nu in (0.5, 5, 50, 300, 1000):
            dps = int(0.75*nu + 60 + 0.1*J)
            kf, big, r = kfin(nu, J, T, 0, dps)
            kl = klr(nu, r)
            pr(f'J={J:5d} Delta=J+{tw} r={float(r):.4f}  nu={nu:7.1f}  Kfin={mp.nstr(kf, 10):>16s}  KLR={mp.nstr(kl, 10):>16s}  '
               f'(Kfin-KLR)*T/(1+nu^2)={mp.nstr((kf-kl)*T/(1+mp.mpf(nu)**2), 4):>10s}  max|term|={mp.nstr(big, 3)}')
    open(os.path.join(FC.CACHE, 'fp_kernel.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
