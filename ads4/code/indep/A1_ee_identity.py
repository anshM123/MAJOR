"""Target A (i): the EE per-nu kernel as a finite sum of b-shifted scalar continuous-spin partial waves.
Independent implementation (statement only: EESTEP_CLOSE §1, LRC_CLOSE §1, H7 §18.1).

Frame (3d metric (-,+,+)): yhat = (1,0,0), xhat = (cosh beta, sinh beta, 0); |x| = 1, |y| = r;
  x1 = -r yhat, x3 = r yhat, x0 = r n_theta, n_theta = (0, cos th, sin th); x2 = -xhat, x4 = xhat.
  celestial z0 = (1, cos ph, sin ph) (rest frame of yhat).
Scalar per-nu integrand (CLS-D j'=0, as written in LRC_CLOSE §1.1, with m_ij = |X_ij|, Nn = V_{0,42} X24, c10 = |z0.x10|,
c30 = |z0.x30|, c13 = |x10.x30| = 2 r^2):
  S_nu = |V_{0,42}|^l (|X24|/|X20 X40|)^((1+l)/2) c13^(5/4 + i nu/2) (c10 c30)^(-delta/2),  l = -1/2 + i nu, delta = 1/2 + i nu.
EE dressing:  P = tau13 * D,  tau13 = (w1.x10)^2 (w3.x30)^2,  D = [(H24 - 2 V_{2,04} V_{4,02})/(X20 X40)]^2 (CPPR conventions).
Claims tested:
  (a) P_even(theta) = sum_{b=-4..4} c_b t^b exactly, t = X20/X40  (so K^EE = sum_b c_b Psi^(b), Psi^(b) = int S_nu t^b)
  (b) Psi^(b)(x) = N(nu) r^p v^q Ghat_{b'}(x)  with p = 2, q = -b/2, b' = -b  (EESTEP §1.4), Ghat from cblocks.py (own recursion)
  (c) end to end: K^EE = N r^2 sum_b c_b v^{-b/2} Ghat_{-b}, N from the scalar kernel at one point.
"""
import os, sys, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cblocks

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A1_ee_identity.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

G3 = np.diag([-1.0, 1.0, 1.0])
def d3(u, v):  # works on (...,3) arrays
    return -u[..., 0] * v[..., 0] + u[..., 1] * v[..., 1] + u[..., 2] * v[..., 2]

# ---------------------------------------------------------------- quadrature nodes
def de_arc(h, tmax):
    """double-exponential nodes on (0, pi) for u in (0, pi): u = pi/2 (1 + tanh(pi/2 sinh t))."""
    t = np.arange(-tmax, tmax + h / 2, h)
    s = np.pi / 2 * np.sinh(t)
    u = np.pi / 2 * (1 + np.tanh(s))
    w = h * np.pi / 2 * np.pi / 2 * np.cosh(t) / np.cosh(s)**2
    # distance to the nearest endpoint, computed without cancellation
    dist = np.pi / 2 * (1 - np.tanh(np.abs(s)))
    dist = np.pi / (1 + np.exp(2 * np.abs(s)))
    keep = dist > 1e-300
    return u[keep], w[keep], dist[keep]

class Frame:
    def __init__(self, r, beta, ntheta=256, h=1/48., tmax=3.6):
        self.r, self.beta = r, beta
        self.yh = np.array([1.0, 0, 0]); self.xh = np.array([math.cosh(beta), math.sinh(beta), 0])
        th = 2 * np.pi * (np.arange(ntheta) + 0.5) / ntheta
        self.th, self.wth = th, 2 * np.pi / ntheta
        u, wu, dist = de_arc(h, tmax)
        # phi = th + u on arc 1, th + pi + u on arc 2; |sin(ph - th)| = sin(u) computed via dist for accuracy
        self.u = np.concatenate([u, u + np.pi]); self.wu = np.concatenate([wu, wu])
        sinu = np.sin(u); sinu = np.where(dist < 1e-3, np.sin(dist), sinu)
        self.abs_sin = np.concatenate([np.abs(sinu), np.abs(sinu)])
        self.geom()

    def geom(self):
        r = self.r
        th = self.th[:, None]; ph = th + self.u[None, :]
        n = np.stack([np.zeros_like(self.th), np.cos(self.th), np.sin(self.th)], -1)       # (T,3)
        x0 = r * n
        x1 = -r * self.yh; x3 = r * self.yh; x2 = -self.xh; x4 = self.xh
        self.x0, self.x1, self.x2, self.x3, self.x4 = x0, x1, x2, x3, x4
        z0 = np.stack([np.ones_like(ph), np.cos(ph), np.sin(ph)], -1)                         # (T,U,3)
        x20 = x2[None, :] - x0; x40 = x4[None, :] - x0                                       # x_i - x_0
        self.X20 = d3(x20, x20); self.X40 = d3(x40, x40); self.X24 = d3(x2 - x4, x2 - x4)
        # V_{0,42} = [(z0.(x4-x0)) X02 - (z0.(x2-x0)) X04]/X24
        V = (d3(z0, x40[:, None, :]) * self.X20[:, None] - d3(z0, x20[:, None, :]) * self.X40[:, None]) / self.X24
        self.V = V
        x10 = x1[None, :] - x0; x30 = x3[None, :] - x0
        self.x10, self.x30 = x10, x30
        self.c13 = np.abs(d3(x10, x10[:, :] * 0 + x30))                                         # = 2 r^2
        # c10 c30 = r^2 sin^2(ph - th) exactly; use the accurate |sin|
        self.c1030 = r * r * self.abs_sin[None, :]**2
        self.t = self.X20 / self.X40
        self.s = 2 * r * math.sinh(self.beta) / (1 - r * r)

    def phi_integral(self, nu):
        """Phi(theta) = int dphi S_nu  (T,)"""
        l = -0.5 + 1j * nu; dl = 0.5 + 1j * nu
        S = (np.abs(self.V)**l * (np.abs(self.X24) / np.abs(self.X20 * self.X40))[:, None]**((1 + l) / 2)
             * self.c13[:, None]**(1.25 + 0.5j * nu) * self.c1030**(-dl / 2))
        return (S * self.wu[None, :]).sum(1)

def null(alpha):
    return np.array([1.0, math.cos(alpha), math.sin(alpha)])
def boost(beta):
    ch, sh = math.cosh(beta), math.sinh(beta)
    return np.array([[ch, sh, 0], [sh, ch, 0], [0, 0, 1.0]])

def embed(x, z=None):
    # P = (P+, P-, x) = (1, x^2, x);  Z = (0, 2 x.z, z);  P.Q = -(P+Q- + P-Q+)/2 + x.y
    if z is None:
        return ('P', x)
    return ('Z', x, z)

def dressing(F, pols):
    """P(theta) = tau13 * D on the theta grid (complex), pols = (a1, a3, a2, a4) angles."""
    a1, a3, a2, a4 = pols
    w1, w3 = null(a1), null(a3)                          # rest frame of yhat
    Bx = boost(F.beta)
    z2, z4 = Bx @ null(a2), Bx @ null(a4)                # boosted with xhat
    tau = d3(w1[None, :], F.x10)**2 * d3(w3[None, :], F.x30)**2
    # embedding invariants: P_i.P_j = -X_ij/2, Z_i.P_j = z_i.(x_j - x_i), Z_i.Z_j = z_i.z_j
    x0 = F.x0; x2, x4 = F.x2, F.x4
    X24 = F.X24; X20 = F.X20; X40 = F.X40
    PP24 = -X24 / 2; PP20 = -X20 / 2; PP40 = -X40 / 2
    Z2P4 = d3(z2, x4 - x2); Z4P2 = d3(z4, x2 - x4); Z2Z4 = d3(z2, z4)
    Z2P0 = d3(z2[None, :], x0 - x2[None, :]); Z4P0 = d3(z4[None, :], x0 - x4[None, :])
    H24 = -2 * (Z2Z4 * PP24 - Z2P4 * Z4P2)
    # V_{i,jk} = [(Z_i.P_j)(P_i.P_k) - (Z_i.P_k)(P_i.P_j)]/(P_j.P_k)
    V2_04 = (Z2P0 * PP24 - Z2P4 * PP20) / PP40
    V4_02 = (Z4P0 * PP24 - Z4P2 * PP40) / PP20
    D = ((H24 - 2 * V2_04 * V4_02) / (X20 * X40))**2
    return tau * D

def laurent_coeffs(F, Pth, deg=8):
    """P_even(theta) (X20 X40)^4 = sum_k q_k cos^k theta (Chebyshev fit) -> Laurent coefficients c_b, b = -4..4."""
    c = np.cos(F.th)
    Peven = 0.5 * (Pth + Pth[::-1])          # theta -> -theta maps node k to N-1-k (midpoint grid symmetric)
    f = Peven * (F.X20 * F.X40)**4
    cheb = np.polynomial.chebyshev.chebfit(c, f, 12)
    tail = np.max(np.abs(cheb[deg + 1:])) / np.max(np.abs(cheb))
    q = np.polynomial.chebyshev.cheb2poly(cheb[:deg + 1])
    s, r = F.s, F.r
    # (1-t)^k (1+t)^(8-k) t^-4 / (256 (1-r^2)^8) s^-k
    P = np.polynomial.polynomial
    cb = np.zeros(9, complex)
    for k in range(deg + 1):
        poly = P.polymul(P.polypow([1, -1], k), P.polypow([1, 1], 8 - k))
        cb += q[k] * s**(-k) * np.asarray(poly, complex) / (256 * (1 - r * r)**8)
    # cb[i] multiplies t^(i-4)
    # direct residual of the Laurent representation on the grid
    rep = sum(cb[i] * F.t**(i - 4) for i in range(9))
    res = np.max(np.abs(rep - Peven)) / np.max(np.abs(Peven))
    odd_over_even = np.max(np.abs(0.5 * (Pth - Pth[::-1]))) / np.max(np.abs(Peven))
    return cb, tail, res, odd_over_even

def zzb(rho, rhb):
    return 4 * rho / (1 + rho)**2, 4 * rhb / (1 + rhb)**2

def main():
    t0 = time.time()
    pts = [(rho, beta) for rho in (0.25, 0.55, 0.85) for beta in (0.6, 1.4, 2.2)]
    nus = [0.4, 1.3, 2.7]
    polsets = [(0.3, 2.1, -0.7, 1.9), (1.2, -2.5, 2.8, 0.4), (-1.0, 0.8, -2.2, -2.9)]
    say('Target A (i): EE kernel = sum_b c_b x shifted scalar partial waves. Points (rho, beta):', pts)
    # quadrature convergence
    rho, beta = 0.85, 0.6; r = rho * math.exp(-beta)
    for nu in (1.3, 2.7):
        Fa = Frame(r, beta, 256, 1 / 48., 3.6); Fb = Frame(r, beta, 512, 1 / 96., 3.9)
        A = (Fa.phi_integral(nu) * Fa.wth).sum(); Bv = (Fb.phi_integral(nu) * Fb.wth).sum()
        say(f'  quadrature self-convergence (rho=0.85, beta=0.6, nu={nu}): rel diff {abs(A/Bv-1):.2e}')
    NMAX = 230
    blocks = {}
    for nu in nus:
        for b in range(0, 5):
            for sgn in (1, -1):
                if b == 0 and sgn == -1: continue
                blocks[(nu, sgn * b)] = cblocks.recursion(nu, sgn * b, NMAX)
    say(f'  blocks: own recursion, N_max = {NMAX}; built in {time.time()-t0:.0f} s')
    # per point data
    data = {}
    for (rho, beta) in pts:
        r = rho * math.exp(-beta); eta = math.cosh(beta)
        F = Frame(r, beta, 384, 1 / 64., 3.8)
        rhb = r * math.exp(-beta)
        z, zb = zzb(rho, rhb)
        v = (1 - z) * (1 - zb)
        for nu in nus:
            Phi = F.phi_integral(nu)
            Psi = {b: (Phi * F.t**b * F.wth).sum() for b in range(-4, 5)}
            P = cblocks.legendre_table(nu, eta, NMAX + 2)
            Gh = {bp: cblocks.Ghat(nu, bp, r, eta, B=blocks[(nu, bp)], P=P) for bp in range(-4, 5)}
            KEE = {}
            for ip, pol in enumerate(polsets):
                Pth = dressing(F, pol)
                cb, tail, res, oe = laurent_coeffs(F, Pth)
                K = (Phi * Pth * F.wth).sum()
                KEE[ip] = (K, cb, tail, res, oe)
            data[(rho, beta, nu)] = dict(r=r, v=v, Psi=Psi, Gh=Gh, KEE=KEE)
    say(f'  quadratures done ({time.time()-t0:.0f} s)')

    # (a) exact reduction
    say('\n(a) exact reduction P_even = sum_b c_b t^b and K^EE = sum_b c_b Psi^(b)')
    mt = max(d['KEE'][ip][2] for d in data.values() for ip in d['KEE'])
    mr = max(d['KEE'][ip][3] for d in data.values() for ip in d['KEE'])
    moe = max(d['KEE'][ip][4] for d in data.values() for ip in d['KEE'])
    say(f'    Chebyshev tail beyond degree 8 of P_even (X20 X40)^4: max {mt:.1e};  Laurent representation residual: max {mr:.1e}')
    say(f'    sin-odd part / even part (pointwise, max): {moe:.2f}  (it integrates to zero against the even Phi)')
    worst = 0; worst_canc = 0
    for key, d in data.items():
        for ip, (K, cb, *_r) in d['KEE'].items():
            pred = sum(cb[i] * d['Psi'][i - 4] for i in range(9))
            canc = sum(abs(cb[i] * d['Psi'][i - 4]) for i in range(9)) / abs(K)
            worst = max(worst, abs(pred / K - 1)); worst_canc = max(worst_canc, canc)
    say(f'    max |sum_b c_b Psi^(b) / K^EE - 1| over 9 points x 3 nu x 3 pol sets = {worst:.2e}; max cancellation in the b-sum {worst_canc:.1e}')

    # (b) identification: ratio test Psi^(b)(x)/Psi^(b)(x)' vs blocks; spread over x of R(x) = Psi^(b)/(r^2 v^{-b/2} Ghat_{b'})
    say('\n(b) identification of Psi^(b): spread over the 9 points of R_b(x) = Psi^(b)(x) / (r^2 v^(-b/2) Ghat_{b\'}(x)), divided by its mean')
    say('    columns: [1] b\' = -b, prefactor v^(-b/2) (claim, direct series);  [2] same object via Euler identity')
    say('             Ghat_{-b} = v^b Ghat_b, i.e. prefactor v^(+b/2) Ghat_{+b} (only for b > 0);  [3] b\' = +b with v^(-b/2) (wrong sign, control);')
    say('             [4] scalar block b\' = 0 with v^(-b/2) (control).  Also |Psi^(b)/Psi^(-b) - 1| (max over points).')
    for nu in nus:
        for b in range(-4, 5):
            row = []
            for (bp, q) in ((-b, -b / 2), (abs(b), abs(b) / 2), (b, -b / 2), (0, -b / 2)):
                R = np.array([data[(rho, beta, nu)]['Psi'][b] /
                              (data[(rho, beta, nu)]['r']**2 * data[(rho, beta, nu)]['v']**q * data[(rho, beta, nu)]['Gh'][bp])
                              for (rho, beta) in pts])
                row.append(np.max(np.abs(R / R.mean() - 1)))
            sym = max(abs(data[(rho, beta, nu)]['Psi'][b] / data[(rho, beta, nu)]['Psi'][-b] - 1) for (rho, beta) in pts)
            say(f'    nu={nu:<4} b={b:+d}:  {row[0]:.2e}   {row[1]:.2e}   {row[2]:.2e}   {row[3]:.2e}    Psi^(b)/Psi^(-b)-1: {sym:.1e}')
    # fitted exponents p, q: log|Psi^(b)/Ghat_{-b}| = p log r + q log v + const
    say('\n    fitted prefactor exponents: log|Psi^(b)/Ghat_{-b}| = p log r + q log v + c (least squares over 9 points)')
    for nu in (1.3,):
        for b in (-4, -3, -1, 0):
            Y = np.array([math.log(abs(data[(rho, beta, nu)]['Psi'][b] / data[(rho, beta, nu)]['Gh'][-b])) for (rho, beta) in pts])
            A = np.array([[math.log(data[(rho, beta, nu)]['r']), math.log(data[(rho, beta, nu)]['v']), 1] for (rho, beta) in pts])
            sol, *_ = np.linalg.lstsq(A, Y, rcond=None)
            say(f'    nu={nu} b={b:+d}: p = {sol[0]:.8f}, q = {sol[1]:.8f} (claim q = {-b/2}), max residual {np.max(np.abs(A@sol-Y)):.1e}')
    # phase constancy of Psi^(b)/(r^2 v^-b/2 Ghat_-b)
    # (c) end to end
    say('\n(c) end to end: K^EE(x) vs N(nu) r^2 sum_b c_b v^(-b/2) Ghat_{-b}(x), N(nu) from the scalar kernel Psi^(0) at (rho,beta) = (0.25,0.6)')
    for ip in range(len(polsets)):
        worst = 0; ctrl = 0; canc = 0
        for nu in nus:
            d0 = data[(0.25, 0.6, nu)]
            N = d0['Psi'][0] / (d0['r']**2 * d0['Gh'][0])
            for (rho, beta) in pts:
                d = data[(rho, beta, nu)]
                K, cb = d['KEE'][ip][0], d['KEE'][ip][1]
                # b = i-4; Psi^(b) = N r^2 v^(-b/2) Ghat_{-b}; for -b < 0 use the Euler identity Ghat_{-b} = v^b Ghat_b
                terms = [cb[i] * (d['v']**(-(i - 4) / 2) * d['Gh'][-(i - 4)] if i - 4 <= 0 else d['v']**((i - 4) / 2) * d['Gh'][i - 4])
                         for i in range(9)]
                pred = N * d['r']**2 * sum(terms)
                worst = max(worst, abs(K / pred - 1))
                canc = max(canc, sum(abs(x) for x in terms) / abs(sum(terms)))
                # control: scalar nu-profile only (all b replaced by the b = 0 block)
                pc = N * d['r']**2 * sum(cb[i] * d['v']**(-(i - 4) / 2) * d['Gh'][0] for i in range(9))
                ctrl = max(ctrl, abs(K / pc - 1))
        say(f'    pol set {ip}: max |K/pred - 1| = {worst:.2e};  cancellation in the b-sum up to {canc:.1e};  control (b=0 block for all b): {ctrl:.2f}')
    say(f'\n(elapsed {time.time()-t0:.0f} s)')
    open(OUT, 'w').write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
