# Independent: which B_k (k=2..9) see contact P = w^i t^L (w=su) or t^i w^L; residue method, from (2.26) directly.
import sympy as sp
s,q=sp.symbols('s q',positive=True)
def mB(k,f):
    br_e=f(s,-q)+f(q-s,-q); br_o=-f(s,-q)+f(q-s,-q)
    e=(2*s-q)*br_e/(s*(s-q))**sp.Rational(k-2,2) if k%2==0 else br_o/(s*(s-q))**sp.Rational(k-3,2)
    return sp.factor(sum(sp.residue(e,s,zz) for zz in [0,q])/2)
for (i,L,kind) in [(1,6,'w^i t^L'),(2,6,'w^i t^L'),(3,6,'w^i t^L'),(1,6,'t^i w^L'),(2,6,'t^i w^L'),(3,6,'t^i w^L')]:
    if kind=='w^i t^L': f=lambda S,U,i=i,L=L:(S*U)**i*(-S-U)**L
    else: f=lambda S,U,i=i,L=L:(-S-U)**i*(S*U)**L
    print(kind,'i=%d L=%d'%(i,L),[ (k,mB(k,f)) for k in range(2,10)])
