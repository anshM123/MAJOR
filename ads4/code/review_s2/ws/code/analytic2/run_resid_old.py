import sys, math
sys.path.insert(0, '../analytic')
import rt_cert as RT
L = float(sys.argv[1])
Er, parts = RT.residual(L)
d = math.exp(-L)
F = RT.fcert_min(RT.beta_of_L(L), RT.beta_of_L(L + 0.25))
print('L', L, 'F', F, 'Er/F', Er/F)
for a, b in parts:
    print('  ', a, b*d**3*(2/RT.LAMf)/F)
