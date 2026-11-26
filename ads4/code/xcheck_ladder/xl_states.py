"""
xcheck_ladder/xl_states.py: the directly computed (light-ray / K2 / RQMC / Gauss-Kronrod) values of lambda_min/N_s for Omega^EE,
transcribed from prl/supplement.tex (Table tab:computed, its continuation, Table tab:Dstates, Table tab:gap2000) with their leverage L
(margin / error budget, margin = the value itself, so budget = value/L) or stated budget.  mass: 'tw' = twist edge Delta = J + Dg
(exact), or m/Lam as printed.  ch: 'min' = 2x2 (++,+-) block of even J (lambda_min), '+-' = 1x1 parity-odd odd-J block,
'OO' = parity-odd even-J block (b = 0 ladder).  For Dg = 2000 the 'img' mass is the chart-A image of the Dg = 1000 twist edge,
m/Lam = sqrt((998 + 2J + 1)/999) (the table prints it to 4 decimals).
Each row: (Dg, J, mass, ch, value, L, budget, stat, note)
"""
S = []


def add(Dg, J, mass, ch, val, L=None, budget=None, stat=None, note=''):
    if budget is None and L is not None:
        budget = val/L
    S.append(dict(Dg=Dg, J=J, mass=mass, ch=ch, val=val, L=L, budget=budget, stat=stat, note=note))


# ---------------- Delta_g = 1000, tab:computed (error-budgeted rows) ----------------
add(1000, 5, 'tw', '+-', 0.20088, budget=2.6e-3, note='binding state; own-mass K2 row, 9 direct nodes [NE+-]')
add(1000, 4, 'tw', 'min', 0.20413, L=22, note='own row; L>=22')
add(1000, 6, 'tw', 'min', 0.22953, L=5.0, note='own row; L>=5.0')
add(1000, 7, 'tw', '+-', 0.26832, budget=0.040, stat=0.00004, note='own direct row')
add(1000, 7, '1.01', '+-', 0.3115, L=7.2)
add(1000, 7, '1.02', '+-', 0.5832, L=10.3)
add(1000, 7, '1.05', '+-', 1.0000, L=25)
add(1000, 8, 'tw', 'min', 0.42973, L=1.2, note='sign-free L=1.2 (lambda-level 4.6)')
add(1000, 9, 'tw', '+-', 0.6757, L=38, stat=0.0006)
add(1000, 9, '1.03', '+-', 0.4543, L=23)
add(1000, 10, 'tw', 'min', 0.5892, L=5.6, stat=0.0006)
add(1000, 11, 'tw', '+-', 0.3485, L=16, stat=0.0006)
add(1000, 12, 'tw', 'min', 0.8528, L=7.4, stat=0.0004)
add(1000, 12, '1.02', 'min', 1.0006, L=7.5)
add(1000, 12, '1.05', 'min', 0.6872, L=3.6)
add(1000, 13, 'tw', '+-', 0.4076, L=6.8, stat=0.0006)
add(1000, 14, 'tw', 'min', 1.0009, L=2.07, stat=0.0006)
add(1000, 14, '1.03', 'min', 0.7264, L=1.52, note='tab:Dstates; L = minimum over the four budgeted masses')
add(1000, 15, 'tw', '+-', 0.6943, L=1.75, stat=0.0006)
add(1000, 16, 'tw', 'min', 0.9402, L=4.4, stat=0.0012)
add(1000, 16, '1.02', 'min', 0.7681, L=3.7)
add(1000, 16, '1.03', 'min', 0.7509, L=2.6)
add(1000, 16, '1.05', 'min', 1.0009, L=1.7)
add(1000, 17, 'tw', '+-', 1.001, L=1.53, stat=0.018)
add(1000, 17, '1.05', '+-', 0.848, L=18.4, note='held at the own row (H_hold [M])')
add(1000, 18, 'tw', 'min', 0.746, L=1.83, stat=0.018, note='structure-basis repaired; unmapped basis gave -0.12')
add(1000, 19, 'tw', '+-', 1.006, L=1.84, stat=0.018)
add(1000, 19, '1.05', '+-', 0.928, L=10.1, note='held at the own row (H_hold [M])')
add(1000, 20, 'tw', 'min', 0.835, L=1.41, stat=0.018, note='chain gave +0.965; direct edge +0.83+-0.14')
add(1000, 21, 'tw', '+-', 0.998, L=1.20, stat=0.020)
add(1000, 21, '1.03', '+-', 0.885, L=1.7)
add(1000, 21, '1.05', '+-', 1.005, L=2.9)
add(1000, 21, '1.075', '+-', 0.968, L=6.2)
add(1000, 21, '1.10', '+-', 0.959, L=6.2)
add(1000, 22, 'tw', 'min', 0.996, L=0.73, stat=0.018, note='[NE+M]')
add(1000, 22, '1.05', 'min', 0.932, L=4.07, note='held')
add(1000, 23, 'tw', '+-', 0.926, L=1.62, stat=0.012, note='lambda_c = 0.904+-0.012')
add(1000, 24, 'tw', 'min', 1.024, L=0.76, stat=0.018, note='[NE+M]')
add(1000, 24, '1.04', 'min', 0.953, L=1.78, note='lambda-level L')
add(1000, 25, 'tw', '+-', 0.960, L=1.65, stat=0.018)
add(1000, 26, 'tw', 'min', 0.985, L=0.95, stat=0.018, note='[NE+M]')
add(1000, 26, '1.03', 'min', 0.956, L=1.24, note='held')
add(1000, 27, 'tw', '+-', 0.976, L=1.81, stat=0.019)
add(1000, 28, 'tw', 'min', 0.950, L=0.66, stat=0.018, note='[NE+M]')
add(1000, 28, '1.03', 'min', 0.948, L=0.65, note='[NE+M]')
add(1000, 29, 'tw', '+-', 0.972, budget=0.784, stat=0.018, note='signed window L=1.24; chain variants -3.62 / +1.39 / +0.948')
add(1000, 30, 'tw', 'min', 0.900, L=0.47, stat=0.018, note='[NE+M]')
add(1000, 31, 'tw', '+-', 0.920, L=0.90, stat=0.032, note='lambda_c = 0.896+-0.032; chain +1.20 -> +1.02 without ref-zero node')
add(1000, 31, '1.04', '+-', 0.920, L=1.35)
add(1000, 31, '1.05', '+-', 0.954, L=1.64)
add(1000, 31, '1.075', '+-', 0.963, L=1.88)
add(1000, 31, '1.10', '+-', 0.995, L=2.13)
add(1000, 32, 'tw', 'min', 0.955, L=0.41, stat=0.018, note='[NE+M]; chain +0.904 -> +0.997 without ref-zero node')
# OO block (parity-odd even J): b = 0 ladder, normalised ++ diagonal
add(1000, 4, 'tw', 'OO', 0.2059, L=22.7)
add(1000, 6, 'tw', 'OO', 0.3551, L=8.0)
add(1000, 8, 'tw', 'OO', 0.5884, L=2.7)
add(1000, 10, 'tw', 'OO', 0.9378, L=1.4)
add(1000, 4, '1.01', 'OO', 0.2150, L=24)

# ---------------- Delta_g = 2000, tab:gap2000 ----------------
for J, ch, tw2, img, mn in [
        (4, 'min', (0.2021, 103), (0.2041, 105.5), None),
        (5, '+-', (0.2004, 320), (0.2009, 334), None),
        (6, 'min', (0.2144, 20.1), (0.2295, 22.0), None),
        (7, '+-', (0.2298, 26.0), (0.2684, 27.5), None),
        (8, 'min', (0.2849, 3.9), (0.4292, 5.2), None),
        (9, '+-', (0.9275, 8.0), (0.6763, 8.2), ('1.03', 0.4544, 9.3)),
        (10, 'min', (0.5820, 1.8), (0.5879, 2.0), ('1.0075', 0.5759, 1.9)),
        (11, '+-', (0.3349, 2.6), (0.3490, 2.6), None),
        (12, 'min', (0.4095, 15.9), (0.8499, 31.8), None),
        (13, '+-', (0.3027, 1.03), (0.4084, 1.03), None),
        (14, 'min', (0.6842, 0.5), (0.9991, 0.6), ('1.035', 0.6992, 0.7)),
        (15, '+-', (0.3659, 0.6), (0.6946, 0.6), None),
        (16, 'min', (0.9996, 4.7), (0.9448, 13.0), ('1.025', 0.7271, 27.7)),
        (17, '+-', (0.5329, 0.4), (1.0000, 0.4), None),
        (19, '+-', (0.8552, 0.2), (0.9989, 0.2), None),
        (20, 'min', (0.9651, 1.35), (0.8948, 3.3), ('1.015', 0.8179, 2.2)),
        (21, '+-', (1.0118, 1.31), (1.0101, 2.5), ('1.03', 0.8909, 4.1)),
        (23, '+-', (1.0446, 1.41), (0.9322, 4.2), ('1.025', 0.9237, 4.9)),
        (27, '+-', (1.0028, 1.20), (1.0365, 6.1), ('1.0175', 0.9443, 2.1)),
        (29, '+-', (0.9407, 0.71), (1.0357, 3.1), ('1.015', 0.9404, 0.77))]:
    add(2000, J, 'tw', ch, tw2[0], L=tw2[1])
    add(2000, J, 'img', ch, img[0], L=img[1])
    if mn:
        add(2000, J, mn[0], ch, mn[1], L=mn[2])
add(2000, 4, 'tw', 'OO', 0.2029, L=106)
add(2000, 4, 'img', 'OO', 0.2059, L=108)
add(2000, 6, 'img', 'OO', 0.3552, L=35.5)
add(2000, 8, 'img', 'OO', 0.5887, L=11.6)
add(2000, 10, 'img', 'OO', 0.9391, L=5.9)
# ---------------- Delta_g = 5000 (tab:gap2000) ----------------
add(5000, 5, 'tw', '+-', 0.2002, L=1922)
add(5000, 5, 'img', '+-', 0.2009, L=2066)
add(5000, 12, 'tw', 'min', 0.2856, L=3.4)
add(5000, 12, 'img', 'min', 0.8482, L=6.4)

# published chain values / artefacts (Dg = 1000), compared for sign only
ARTEFACTS = [
    (1000, 18, 'tw', 'min', -0.12, 'unmapped structure basis (S6, "basis artefact")'),
    (1000, 29, 'tw', '+-', -3.62, 'w2 odd chain with reference-zero node ("chain artefact")'),
    (1000, 29, 'tw', '+-', 1.39, 'spline through 5 direct edge nodes (aliasing)'),
    (1000, 29, 'tw', '+-', 0.948, 'published chain (node dropped)'),
    (1000, 31, 'tw', '+-', 1.20, 'published chain with reference-zero node'),
    (1000, 31, 'tw', '+-', 1.02, 'chain without that node'),
    (1000, 32, 'tw', 'min', 0.904, 'published chain with reference-zero node'),
    (1000, 32, 'tw', 'min', 0.997, 'chain without that node'),
    (1000, 20, 'tw', 'min', 0.965, 'chain value before the direct row'),
]
