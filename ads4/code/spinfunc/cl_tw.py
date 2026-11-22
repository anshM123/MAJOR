"""cl_tw: twist-edge mass m/Lam (rounded up at 1e-4) of spin J at Delta_g (state Delta = Delta_g + J).  usage: python3 cl_tw.py Dg J"""
import sys, math
Dg, J = int(sys.argv[1]), int(sys.argv[2])
L = math.sqrt((Dg - 2)*(Dg - 1))
print('%.4f' % (math.ceil(math.sqrt((Dg - 2)*(Dg - 1 + 2*J))/L*1e4)/1e4))
