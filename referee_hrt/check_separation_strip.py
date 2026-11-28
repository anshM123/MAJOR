# (1) Separation Lemma S8.1: for an acausal circle (graph of increasing homeo phi of RP^1) and complementary arcs
#     alpha, beta cut at r,r', any positive combos X in K(alpha), Y in K(beta) satisfy N(X,Y) >= 2 sqrt(detX detY).
# (2) Explicit invariant strip at a rotating horizon, Eqs. (S27)-(S28): det X = 1, h X(s,r) = X(s+L,r),
#     induced metric dr^2 + [1+(1-q^2) sinh^2 r] ds^2.
import numpy as np, sympy as sp
rng=np.random.default_rng(3)
def N(X,Y):  # tr(adj(X) Y)
    adj=np.array([[X[1,1],-X[0,1]],[-X[1,0],X[0,0]]]); return np.trace(adj@Y)
worst=np.inf
for trial in range(4000):
    # random increasing homeo of circle: theta -> phi(theta), with r,r' at random positions
    k=rng.integers(1,4); c=rng.normal(size=k)*0.3/np.arange(1,k+1); ph0=rng.uniform(0,2*np.pi)
    phi=lambda t: t+sum(c[j]*np.sin((j+1)*t+ph0) for j in range(k))/ (1+0*t)
    # ensure monotone: derivative 1+sum c_j (j+1) cos >0
    if 1-sum(abs(c[j])*(j+1) for j in range(k))<=0.05: continue
    p=lambda t: np.outer([np.cos(t/2),np.sin(t/2)],[np.cos(phi(t)/2),np.sin(phi(t)/2)])  # RP^1 as angle/2
    r1,r2=np.sort(rng.uniform(0,2*np.pi,2))   # t in [0,2pi): theta=t/2 covers RP^1 once
    ta=rng.uniform(r1,r2,rng.integers(2,5)); tb=rng.uniform(r2,r1+2*np.pi,rng.integers(2,5))
    X=sum(rng.uniform(0.1,2)*p(t) for t in ta); Y=sum(rng.uniform(0.1,2)*p(t) for t in tb)
    dX,dY=np.linalg.det(X),np.linalg.det(Y)
    if dX<=0 or dY<=0: print('positivity failed',dX,dY); break
    ratio=N(X,Y)/(2*np.sqrt(dX*dY)); worst=min(worst,ratio)
print('Separation lemma: min N/(2 sqrt(detX detY)) over random trials =',worst)
s,r,a,b=sp.symbols('s r a b',positive=True); L=(a+b)/2; q=(a-b)/(a+b)
c=sp.Matrix([[0,sp.exp(s)],[-sp.exp(-s),0]]); n=sp.Matrix([[sp.exp(q*s),0],[0,-sp.exp(-q*s)]])
X=sp.cosh(r)*c+sp.sinh(r)*n
gL=sp.diag(sp.exp(a/2),sp.exp(-a/2)); gR=sp.diag(sp.exp(b/2),sp.exp(-b/2))
print('det X - 1 =',sp.simplify(X.det()-1))
R=(gL*X*gR.inv()-X.subs(s,s+L)).applyfunc(lambda e: sp.simplify(sp.powsimp(sp.expand(e),combine='exp').rewrite(sp.exp)))
print('equivariance residual =',R.applyfunc(lambda e: sp.simplify(sp.logcombine(e))))
import random
for _ in range(3):
    vals={a:random.uniform(.5,4),b:random.uniform(.5,4),s:random.uniform(-2,2),r:random.uniform(-2,2)}
    print('  numeric equivariance residual:',max(abs(complex(x)) for x in (gL*X*gR.inv()-X.subs(s,s+L)).subs(vals).evalf()))
def metric(u,v):  # polarisation of -det
    return sp.simplify(-( (u+v).det()-u.det()-v.det())/2)
Xs=X.diff(s); Xr=X.diff(r)
print('g_rr =',sp.simplify(metric(Xr,Xr)),' g_rs =',sp.simplify(metric(Xr,Xs)),
      ' g_ss - [1+(1-q^2)sinh^2 r] =',sp.simplify(sp.expand(metric(Xs,Xs)-(1+(1-q**2)*sp.sinh(r)**2)).rewrite(sp.exp)))
