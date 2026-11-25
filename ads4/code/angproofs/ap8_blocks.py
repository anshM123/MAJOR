"""Block lower bound (Lemma B2):  for |delta| <= 1, S >= 4, h in [S, 2S]:
    k^{(delta)}_{2h}(1 - S^{-2}) := z^h 2F1(h+delta, h+delta; 2h; z) >= c(delta) 4^h S^{1/2 + 2 delta},
    c(delta) = (1 - 2/S0) (1 - 1/(S0-1)) (1 + 1/(2 S0))^{-1/2} (2 sqrt(pi))^{-1} (1 - 3/S0)^{2 S0} 2^{-delta-1} 3^{-max(delta,0)}
at S0 = 4 (all factors increase with S0).  Numerical check of the ratio for delta in {-1,-1/2,0,1/2,1}."""
import mpmath as mp
mp.mp.dps = 30

def c_delta(delta, S0=4):
    S0 = mp.mpf(S0)
    return ((1 - 2/S0)*(1 - 1/(S0 - 1))*(1 + 1/(2*S0))**-0.5/(2*mp.sqrt(mp.pi))*(1 - 3/S0)**(2*S0)
            * mp.mpf(2)**(-delta - 1)*mp.mpf(3)**(-max(delta, 0)))

def k(h, delta, z):
    return z**h*mp.hyp2f1(h + delta, h + delta, 2*h, z)

if __name__ == '__main__':
    for delta in (-1, -0.5, 0, 0.5, 1):
        cmin = mp.inf
        for S in (4, 5, 8, 16, 64, 256, 1024, 4096):
            z = 1 - mp.mpf(S)**-2
            for f in (1, 1.25, 1.5, 1.75, 2):
                h = mp.mpf(S)*f
                r = k(h, delta, z)/(4**h*mp.mpf(S)**(0.5 + 2*delta))
                cmin = min(cmin, r)
        print(f"delta = {delta:5}: min ratio over grid = {mp.nstr(cmin, 5)}   proven constant c(delta) = {mp.nstr(c_delta(delta), 5)}"
              f"   bound holds: {cmin >= c_delta(delta)}")
    for S0 in (8, 16):
        print(f"   (constant at S0 = {S0}, delta = 0: {mp.nstr(c_delta(0, S0), 5)})")
