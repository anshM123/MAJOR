import sys, math, os, json
sys.path.insert(0, '../analytic')
import rt_loc3 as R3
La = float(sys.argv[1]); Lb = float(sys.argv[2]); J = int(sys.argv[3])
lI2, lI3 = R3.integrals(La, J, False, Lb)
print('La', La, 'Lb', Lb, 'J', J, 'lI2', lI2, 'lI3', lI3, flush=True)
