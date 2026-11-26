from contour_direct import *
from check_norm import parameters
c=json.load(open('boxes_23_p3.json'))[0];d,dd,*_=parameters(c);dl,du=d-dd,d+dd
for b in [4]:
 for y in map(F,['0','0.25','0.5','1','2','3','4','5','6','8','10','15','20','25']):
  u=y+F('.25');P,Q=abs_transfer(dl,du,RIGHT,RIGHT,y,u,4,b);v=hi((u-y)*osc_integrand(dl,du,4,b,RIGHT,RIGHT,y,u,P,Q));print(float(y),float(v),flush=True)
 for Y in map(F,['20.25','30.25','40.25']):
  t=time.perf_counter();print('top',float(Y),float(top_abs(dl,du,4,b,Y,step=F(25))[0]),time.perf_counter()-t,flush=True)
