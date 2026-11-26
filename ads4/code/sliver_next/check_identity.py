#!/usr/bin/env python3
"""Exact Gaussian-rational identity checks; algebraic regression, not continuum coverage."""
from fractions import Fraction as F
from dataclasses import dataclass
@dataclass(frozen=True)
class C:
    r:F=F(0); i:F=F(0)
    def __add__(self,o):o=cast(o);return C(self.r+o.r,self.i+o.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,o):return self+-cast(o)
    def __rsub__(self,o):return cast(o)+-self
    def __mul__(self,o):o=cast(o);return C(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
    __rmul__=__mul__
    def __truediv__(self,o):o=cast(o);den=o.r**2+o.i**2;return self*C(o.r/den,-o.i/den)
    def __rtruediv__(self,o):return cast(o)/self
    def __pow__(self,n):
        if n<0:return 1/(self**(-n))
        q=C(F(1));v=self
        while n:
            if n&1:q=q*v
            v=v*v;n//=2
        return q

def cast(x):return x if isinstance(x,C) else C(F(x))

def run():
    count=0
    for b in (0,4):
      s=2*b+1;q=b*(b+1)
      for n in (2,5,8,11):
       for a in (F(7,2),F(998),F(1997,2)):
        N=n+a;D=N*(N+s);h=1/D;d=1-F(n*(n+s),1)/D;tau=h/(1+q*h)
        U=[F(1)]
        for k in range(n):U.append(U[-1]*F(n*(n+s)-k*(k+s),1)/(D-k*(k+s)))
        def coeff(beta):
            f=[cast(1)]
            for k in range(n+1):f.append(f[-1]*(beta+k)/(k+1))
            return f
        def S(beta):return sum(u*v for u,v in zip(U,coeff(beta)))
        for p in (3,4):
         alpha=C(F(2*p-1,2),F(17,3));nu=C(F(0),F(17,3));cc=cast(1)
         # Z_j/d^alpha = c_j S(alpha-j); Z_-1/d^alpha = S(alpha+1)
         zm=S(alpha+1);z=S(alpha);rp=cast(1);P0=cast(1);Pprev=cast(0);Q0=cast(0);Qprev=cast(1);P=cast(0);Q=cast(0);direct=cast(0)
         for j in range(12):
          beta=alpha-j;rr=(F(2*j+1,2)-nu)**2/((j+1)*(j+1-2*nu))
          A=(d+h*beta*(s-2*beta-1))*rr;BB=h*beta*(beta+1)*rr*rp;den=1-h*beta*(beta-s)
          assert z==cc*S(beta)
          f=coeff(beta);A0=F(9,2)-b*(b-1)-beta*(beta+8-4*b)/2;A1=beta*(beta+9-4*b)/2
          dressed=sum(U[k]*(f[k]+tau*((k+1)**2*f[k+1]/2+F(-k*k+(8-4*b)*k+9-2*b*(b-1),2)*f[k])) for k in range(n+1))
          assert dressed==(1+tau*A0)*S(beta)+tau*A1*S(beta+1)
          if p==3:
            C0=1+q*h+h*A0;C1=h*A1*rp;P+=C0*P0+C1*Pprev;Q+=C0*Q0+C1*Qprev;direct+=cc*dressed*(1+q*h)
          else:P+=P0;Q+=Q0;direct+=cc*S(beta)
          znext=(A*z+BB*zm)/den;assert znext==cc*rr*S(beta-1)
          zm,z=z,znext;Pprev,P0=P0,(A*P0+BB*Pprev)/den;Qprev,Q0=Q0,(A*Q0+BB*Qprev)/den;cc*=rr;rp=rr;count+=1
         assert P*S(alpha)+Q*S(alpha+1)==direct
    print(f'PASS: {count} exact recurrence and dressing checks, 48 signed transfer sums.')
if __name__=='__main__':run()
