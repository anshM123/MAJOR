import sys; sys.dont_write_bytecode=True
import ls_core as C, time
for J,mfl,h in [(100,C.mgap(1000,100),'+-'),(40,2.0,'++'),(300,1.3,'+-'),(5,C.mgap(1000,5),'+-')]:
    t=time.time(); a=C.hyb(1000,J,mfl,h,npn=512); t1=time.time(); b=C.hyb_fast(1000,J,mfl,h); t2=time.time()
    print(J,h,'Ns',a['Ns'],b['Ns'],'cv2',a['cv2'],b['cv2'],'cv3',a['cv3'],b['cv3'],'times',t1-t,t2-t1)
