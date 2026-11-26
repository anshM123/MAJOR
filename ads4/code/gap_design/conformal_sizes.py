# Sizes of the conformal coefficients Atilde_j = sum_{k=1}^j 4^k (-1)^(j-k) C(j+k-1,j-k) A_k (A~_0=A_0)
# computed exactly from table midpoints; reliability limited by midpoint error * amplification.
import sys,math
sys.set_int_max_str_digits(0)
J=int(sys.argv[1]);fn=sys.argv[2]
f=open(fn);N,D=map(int,f.readline().split())
A=[];R=[]
for i,l in enumerate(f):
    if i>J:break
    a,r,g,t=l.split();A.append(int(a));R.append(int(r))
def l10(x):
    x=abs(x)
    return -999 if x==0 else (x.bit_length()-1)*math.log10(2)+math.log10(x/(1<<(x.bit_length()-1))) if x.bit_length()<1000 else math.log10(x>>(x.bit_length()-60))+(x.bit_length()-60)*math.log10(2)
out=[]
for j in range(1,J+1):
    # C(j+k-1,j-k) for k=1..j ; iterate k downward from j: C(2j-1,0)=1
    s=0;amp=0;bn=1  # k=j: C(2j-1,0)=1
    for k in range(j,0,-1):
        term=bn*(4**k)
        s+=term*A[k] if (j-k)%2==0 else -term*A[k]
        amp+=term*R[k]
        # next k-1: C(j+k-2, j-k+1) = C(j+k-1,j-k)*(j-k+1)... use direct update
        if k>1:
            n_old=j+k-1;r_old=j-k   # C(n_old,r_old) -> C(n_old-1,r_old+1)
            bn=bn*(n_old-r_old-1)*(n_old-r_old)//((r_old+1)*n_old) if n_old>0 else 0
    out.append((j,l10(s)-D,l10(amp)-D))
    if j%100==0 or j<5:print(j,'log10|At_j| %.1f  log10 err %.1f'%(l10(s)-D,l10(amp)-D),flush=True)
import json;json.dump(out,open('conformal_sizes.json','w'))
