import sys,math
sys.set_int_max_str_digits(0)
def rd(fn):
 f=open(fn);N,D=map(int,f.readline().split());rows=[list(map(int,l.split())) for l in f];assert len(rows)==N+1;return N,D,rows
N,D,new=rd('direct_coefficients.txt');N0,D0,old=rd('direct_coefficients_1253_1600.txt')
print('new N',N,'D',D,'old N',N0,'D0',D0)
def l10(x):return -9999 if x==0 else math.log10(x) if x<10**300 else len(str(x))-1+math.log10(int(str(x)[:15])/10**14)
mx=-1e9;mg=-1e9
for k,(a,r,g,t) in enumerate(new):
 ra=l10(r)-D;rg=l10(t)-D;mx=max(mx,ra);mg=max(mg,rg)
 if k%300==0 or k==N:print(k,'log10|A|',round(l10(abs(a))-D,1),'log10 rA',round(ra,1),'log10 rG',round(rg,1))
print('max log10 radius A',round(mx,1),'G',round(mg,1))
sc=10**(D-D0);bad=0;agree=[]
for k,((a,r,g,t),(a0,r0,g0,t0)) in enumerate(zip(new,old)):
 for m,e,m0,e0 in ((a,r,a0,r0),(g,t,g0,t0)):
  if abs(m-m0*sc)>e+e0*sc+sc:bad+=1
 agree.append(l10(abs(a-a0*sc)+1)-D)
print('overlap 0..%d intervals intersect: bad=%d; max log10|midpoint diff|=%.1f'%(N0,bad,max(agree)))
