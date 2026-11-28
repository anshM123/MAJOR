import sys, n0_sector as m
for r in [float(x) for x in sys.argv[1:]]:
    print(r, repr(m.T_of_r(r, 60.0, 16)), flush=True)
