from racah_explore import *
import sys
g=int(sys.argv[1]); s=int(sys.argv[2]); lo=int(sys.argv[3]); hi=int(sys.argv[4])
ts=[0,0.25,0.5,1,2,4,8,16,32,64,128,256,1000,10000]
for n in range(lo,hi):
    w,wh,L=scan(g,n,s,Jmax=6*g,ts=ts)
    tau=4*mp.sqrt(3/L)
    print(g,s,n,'n/g=%.3f'%(n/g),'minmargin=%.4f'%float(w),wh,'L=%.4g'%float(L),'tau=%.4g'%float(tau),'tau*g^3=%.1f'%float(tau*g**3),flush=True)
