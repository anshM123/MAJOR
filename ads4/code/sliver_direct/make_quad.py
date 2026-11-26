from fractions import Fraction as F
from math import lcm
from check_norm import cosine_nodes,NODE
N=96;xs=[int(NODE*x) for x in cosine_nodes(N)]
# exact interpolatory quadrature on rational cosine nodes
p=[1]
for x in xs:
 q=[0]*(len(p)+1)
 for k,a in enumerate(p):q[k]-=x*a;q[k+1]+=a
 p=q
w=[];L=lcm(*range(1,N+2))
for x in xs:
 q=[0]*N+[p[-1]]
 for k in range(N-1,-1,-1):q[k]=p[k+1]+x*q[k+1]
 assert p[0]+x*q[0]==0
 den=1
 for y in xs:
  if y!=x:den*=x-y
 num=sum(2*q[k]*NODE**k*(L//(k+1))for k in range(0,N+1,2))
 w.append(F(num,L*den))
# Each weight is exact by integrated cardinal Lagrange polynomial.
# Bound total variation without adding huge unrelated rational denominators.
S=10**80
norm=sum((abs(a.numerator)*S+a.denominator-1)//a.denominator for a in w)
assert norm<3*S
print('PASS quadrature',N+1,'nodes, norm upper',norm,'scale',S)
with open('quad.txt','w') as f:
 f.write(str(N+1)+'\n')
 for x,a in zip(xs,w):f.write(f'{x} {NODE} {a.numerator} {a.denominator}\n')
