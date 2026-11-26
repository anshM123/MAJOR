# Non-rigorous float probe (mpmath) of the relaxed prefix sum S_0(c,h) on the cell rectangle.
import sys,mpmath as mp
sys.set_int_max_str_digits(0)
mp.mp.dps=int(__import__("os").environ.get("DPS","900"))
f=open('direct_coefficients.txt');N,D=map(int,f.readline().split())
A=[];G=[]
for l in f:
 a,r,g,t=l.split();A.append(mp.mpf(a)/mp.mpf(10)**D);G.append(mp.mpf(g)/mp.mpf(10)**D)
def S0(c,h,b,K):
 s=2*b+1;q=b*(b+1);sg=1 if b==0 else -1;S=0
 for k in range(K,-1,-1):
  k1=(k+1)**2;k2=-k*k+(8-4*b)*k+9-2*b*(b-1)
  B=(k1*A[k+1]+k2*A[k])/2
  Fk=(1+q*h)*A[k]+h*(B+sg*G[k])
  p=k*(k+s)
  S=Fk+(c-p*h)/(1-p*h)*S if k<K else Fk
 return S*(1-c)**3/(1+q*h)**4
if __name__=='__main__':
 cl,cu,b,K=mp.mpf(sys.argv[1]),mp.mpf(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
 hc=cu/((1001-b)*(1002+b))
 for c in (cl,cu):
  print('c',float(c),' '.join('%.2e'%float(S0(c,hc*t,b,K)) for t in [0,.25,.5,.75,.9,1.0]))
def caps(cl,cu,b):
 return cu/((1001-b)*(1002+b)), (1-mp.sqrt(cl))**2/998**2
