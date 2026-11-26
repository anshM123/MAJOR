"""rv15_rectangle.py -- the note's J<=32 certificate box is the full rectangle (q,u) in [0,1]^2 (q = 997002 t, u = L^2/m^2),
which includes non-physical masses Lambda <= m < m_gap.  Scan the closed form on that rectangle for J<=64 (numerical only),
to confirm J<=32 is positive on the whole box and to locate the non-physical failures the note reports for J>=33."""
from flint import arb, fmpq, ctx
from rv7_sweep import AGtab, S_value
qs = [fmpq(i, 8) for i in range(9)]
ts = [q/997002 if q != 0 else fmpq(0) for q in qs]
tab = AGtab(66, ts, 600)
neg = []; mins = {}
for J in range(0, 65):
    for b in (0, 4):
        if (b == 0 and J % 2) or J < b:
            continue
        mn = None
        for t in ts:
            for i in range(0, 201):
                u = arb(fmpq(i, 200))
                S = S_value(J, b, t, u, tab)
                v = float(S.mid())
                if mn is None or v < mn[0]:
                    mn = (v, float(t*997002), i/200)
                if not S > 0:
                    neg.append((J, b, float(t*997002), i/200, v))
        mins[(J, b)] = mn
bad32 = [x for x in neg if x[0] <= 32]
print("J<=32: non-positive points on the full [0,1]^2 rectangle:", len(bad32))
print("J<=32 smallest value on the rectangle:", min((v, k) for k, v in mins.items() if k[0] <= 32))
Js = sorted(set(x[0] for x in neg))
print("J with non-positive points on the rectangle (non-physical u > u_edge):", Js)
for J in Js[:6]:
    pts = [x for x in neg if x[0] == J]
    print(f"  J={J}: {len(pts)} points, e.g. {pts[:3]}; smallest u among them = {min(x[3] for x in pts)}")
# for those points check they are non-physical: u_edge(q) = L^2/m_gap^2 with a_g from q: q = 999 s^2/(998+s) etc.
