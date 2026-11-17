"""s3 (referee M3): sgn(v'-u') on the frame domain of the moments.
Frame (x1,x2,x3,x4) = (-y, x/x^2, y, -x/x^2), x, y future timelike (the y>0, x>0 Laplace domain of the M-tables).
u' = x12^2 x34^2/(x13^2 x24^2), v' = x14^2 x23^2/(x13^2 x24^2) (CLS-D (31); 1<->3 swaps u' and v').
Claim [D]: v'-u' = eta(1+r^2)/(2r) with eta = -xhat.yhat >= 1, r = |x||y| -> v'-u' > 0 on the whole domain; it is a cross-ratio, hence independent of the
polarisations (celestial angle).  Under 1<->3 (y -> -y, past-directed: outside the domain) eta -> -eta.  Checked on random configurations."""
import numpy as np
rng = np.random.default_rng(7)
g = np.diag([-1.0, 1, 1]); dot = lambda a, b: a @ g @ b
def timelike():
    v = rng.normal(size=2)*rng.uniform(0, 3); t = np.sqrt(1 + v @ v)
    return np.array([t, v[0], v[1]])
worst = 0.0; npos = 0; N = 20000
for _ in range(N):
    xh, yh = timelike(), timelike()
    X, Y = rng.uniform(0.05, 3)*xh, rng.uniform(0.05, 3)*yh
    x1, x2, x3, x4 = -Y, X/dot(X, X), Y, -X/dot(X, X)
    s = lambda a, b: dot(a - b, a - b)
    up = s(x1, x2)*s(x3, x4)/(s(x1, x3)*s(x2, x4)); vp = s(x1, x4)*s(x2, x3)/(s(x1, x3)*s(x2, x4))
    r = np.sqrt(dot(X, X)*dot(Y, Y)); eta = -dot(X, Y)/r
    pred = eta*(1 + r*r)/(2*r)
    worst = max(worst, abs((vp - up) - pred)/abs(pred))
    npos += (vp - up) > 0
print(f"s3: {N} random frame configurations: max rel |(v'-u') - eta(1+r^2)/(2r)| = {worst:.2e};  v'-u' > 0 in {npos}/{N};  min eta = 1 by construction")
# 1<->3 swap: y -> -y
X, Y = 1.3*timelike(), 0.4*timelike()
x1, x2, x3, x4 = Y, X/dot(X, X), -Y, -X/dot(X, X)
s = lambda a, b: dot(a - b, a - b)
up = s(x1, x2)*s(x3, x4)/(s(x1, x3)*s(x2, x4)); vp = s(x1, x4)*s(x2, x3)/(s(x1, x3)*s(x2, x4))
print(f"    after 1<->3 (y -> -y): v'-u' = {vp - up:.4f} (negative: the image lies on the other sheet)")
