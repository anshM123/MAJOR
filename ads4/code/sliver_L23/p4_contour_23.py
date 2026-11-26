from contour_direct import *
from check_norm import parameters
p=4;Y=F('30.25')
for row in json.load(open(sys.argv[1])):
 b=row['b']; assert row['status']=='ENCLOSED' and row['p']==p
 c=row['config'];d,dd,h,dh,v,dv,im=parameters(c);dl,du=d-dd,d+dd
 assert c['fullTwist'] and c.get('nuAxis')=='imag' and v==RIGHT and im-dv<=0 and im+dv>=25
 cap_check(c)
 P,Q=F(row['P_bound']),F(row['Q_bound']);t=time.perf_counter();bot,left=near_contours(dl,du,p,b,Y)
 top=top_abs(dl,du,p,b,Y,F(25))[0];right=F(0);y=F(0);worst=(F(0),None)
 while y<Y:
  u=min(Y,y+F('.25'));R,T=(P,Q) if u<=25 else abs_transfer(dl,du,RIGHT,RIGHT,y,u,p,b)
  e=hi((u-y)*osc_integrand(dl,du,p,b,RIGHT,RIGHT,y,u,R,T));right+=e
  if e>worst[0]:worst=(e,(str(y),str(u)))
  y=u
 x=extra(dl,du,p,b)
 r=dict(p=p,b=b,Y=str(Y),dlo=str(dl),dhi=str(du),bottom=str(bot),left=str(left),top=str(top),right=str(right),extra=str(x),total=str(bot+left+top+right+x),seconds=time.perf_counter()-t)
 print(json.dumps({**r,**{k+'_approx':float(F(r[k])) for k in ['bottom','left','top','right','extra','total']},'worst_approx':float(worst[0]),'worst_interval':worst[1]}),flush=True)
 Path(f'ct_23_4_{b}.json').write_text(json.dumps([r],indent=2)+'\n')
