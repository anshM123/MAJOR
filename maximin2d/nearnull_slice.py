"""Crossing density d_s d_t D on boundary slices that are nearly null at the two endpoints
(slope +S at P, -S at Q), NEC vs NEC-violating Vaidya-AdS3; plus a step-size check.
usage: python3 -I nearnull_slice.py"""
import importlib.util, os
spec = importlib.util.spec_from_file_location("j", os.path.join(os.path.dirname(os.path.abspath(__file__)), "jacobi_sign_vaidya.py"))
j = importlib.util.module_from_spec(spec); spec.loader.exec_module(j)
P, Q = (0.0, 0.0), (0.4, 1.5)
for name, prm in [("vacuum", (1.0, 0.0, 0.0, 0.3)), ("NEC dm=+1.5", (1.0, 1.5, 0.0, 0.3)),
                  ("NEC-viol dm=-0.8", (1.0, -0.8, 0.0, 0.3))]:
    for S in (0.9, 0.99):
        nu = j.mixed(P, Q, (S, 1.0), (-S, 1.0), 0.02, prm)
        print(f"{name:18s} slope {S}: nu-density = {nu:+.4e}")
prm = (1.0, 1.5, 0.0, 0.3)
for h in (0.04, 0.02, 0.01):
    print("step check NEC dm=+1.5, D+- with h =", h, ":", f"{j.mixed(P, Q, (1.0, 1.0), (1.0, -1.0), h, prm):+.5e}")
