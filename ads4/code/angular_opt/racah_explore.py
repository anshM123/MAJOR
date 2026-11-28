"""Explore the Racah-family dual: lam_k D_k^2 = A (-n)_k (n+s)_k / (k! (s+1)_k).
Then Lambda^(0) = A 4F3(-n,n+s,-J,J+1; s+1, 2-Delta, Delta-1; 1) and Lambda^(4) = A 4F3(-n,n+s,4-J,J+5; s+1,6-Delta,Delta+3;1)."""
import mpmath as mp, sys
mp.mp.dps = 60

def Lam(b, J, D, n, s):
    # sum_k c_k * prod_{i<k} (C - y_{b+i})/(T - y_{b+i})
    C = mp.mpf(J*(J+1)); T = (D-1)*(D-2)
    tot = mp.mpf(1); term = mp.mpf(1)
    for k in range(n):
        if k+1 > J-b: break
        y = (b+k)*(b+k+1)
        term = term*(-n+k)*(n+s+k)/((k+1)*(s+1+k))*(C-y)/(T-y)
        tot += term
    return tot

def Wf(b, D):
    T = (D-1)*(D-2); return (T - b*(b-1))/mp.mpf(3) + mp.mpf(3)/2

def scan(g, n, s, Jmax, ts):
    g = mp.mpf(g)
    DB = g+4; A = -Wf(4, DB)
    worst = mp.inf; where=None
    for J in range(0, Jmax+1):
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            for t in ts:
                D = J + g + t
                v = A*Lam(b, J, D, n, s)
                m = (Wf(b, D) + v) if b == 4 else (Wf(b, D) - v)
                m = m/Wf(b, D)
                if m < worst: worst = m; where = (b, J, float(t))
    c2 = mp.rf(-n, 2)*mp.rf(n+s, 2)/(2*mp.rf(s+1, 2))
    L = -A*c2/100
    return worst, where, L

if __name__ == '__main__':
    g = int(sys.argv[1]); s = int(sys.argv[2])
    ts = [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 1000, 10000]
    for n in range(2, int(1.2*g)):
        w, wh, L = scan(g, n, s, Jmax=8*g, ts=ts)
        tau = 4*mp.sqrt(3/L)
        print(g, s, n, 'minmargin=%.4f'%float(w), wh, 'L=%.4g'%float(L), 'tau=%.4g'%float(tau), 'tau*g^3=%.1f'%float(tau*g**3), flush=True)
