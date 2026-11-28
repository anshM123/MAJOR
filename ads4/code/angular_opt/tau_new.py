"""New angular bound (A2' limit, D and the light term set to zero):
   tau^2 <= 3200 (3+D) / ( W_B(g) * M (M-1) * S ),   M = ceil(g)-2,
valid whenever an exact certificate (cache/certI_*.json, ok=True) has Mb >= M, n <= M and theta_a <= theta(g).
Prints the bound at a list of gaps and compares with the paper's tau_*."""
import json, glob, math
from fractions import Fraction as Fr
from cert_lemmaM import WB, theta
certs = []
for f in glob.glob('cache/certI_*.json'):
    d = json.load(open(f))
    if d['ok']: certs.append((d['Mb'], d['n'], Fr(d['theta_a']), Fr(d['S']), f))
def tau_new(g):
    g = Fr(g); M = math.ceil(g) - 2; best = None
    for Mb, n, th, S, f in certs:
        if Mb >= M and n <= M and th <= theta(g):
            t2 = Fr(3200*3)/(WB(g)*M*(M-1)*S)
            if best is None or t2 < best[0]: best = (t2, f)
    return best
def tau_paper(g):
    K = 729/4800*(g-2)**2*(1-49/(2*g+9)**2); return 4*math.sqrt(3/K)
if __name__ == '__main__':
    for g in [32, 50, 100, 300, 1000, 3000, 10000]:
        r = tau_new(g)
        if r is None: print(g, 'no certificate'); continue
        t = math.sqrt(float(r[0]))
        print(f"g={g}: tau_new <= {t:.5e}   tau_new*g^2 = {t*g*g:.2f}   paper tau_* = {tau_paper(g):.5e}   ratio = {tau_paper(g)/t:.1f}   [{r[1]}]")
