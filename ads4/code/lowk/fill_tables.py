"""Refreshes the result tables of ads4/SPINNING_LOWK.md (sect. 5.1, 5.2 from out_t3_tables.txt; sect. 7 from out_t2_corrections.txt).
Usage: python3 fill_tables.py"""
import os
import re
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, '..', '..') + '/'
L = HERE + '/'
doc = open(D + 'SPINNING_LOWK.md').read()
tab = open(L + 'out_t3_tables.txt').read()
def section(l):
    parts = re.split(r'(?=^### l = )', tab, flags=re.M)
    return "\n\n".join(p.replace('### ', '#### ').rstrip() for p in parts if p.startswith(f'### l = {l}, k = '))
def between(text, a, b, new):
    i = text.index(a) + len(a); j = text.index(b, i)
    return text[:i] + "\n\n" + new + "\n\n" + text[j:]
doc = between(doc, "### 5.1 Currents, `l = 1`", "### 5.2 Stress tensors", section(1))
doc = between(doc, "### 5.2 Stress tensors, `l = 2`", "## 6. Non-analytic", section(2))
t2 = open(L + 'out_t2_corrections.txt').read().split('# ---- summary')[1].splitlines()
rows = ["| l | β | k | ν | comp | c_2: B2 [χ²] | B3 [χ²] | A [χ²] | N [χ²] | **adopted c_2** | c_1 (A) | c_{k+½} (N) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for line in t2:
    if line.startswith('#') or not line.strip():
        continue
    f = line.split()
    if f[2] not in ('2', '3'):
        continue
    # f: l beta k nu comp  B2 [chi] B3 [chi] A [chi] N [chi] adopted c1 cN
    rows.append(f"| {f[0]} | {f[1]} | {f[2]} | {f[3]} | {f[4]} | {f[5]} {f[6]} | {f[7]} {f[8]} | {f[9]} {f[10]} | {f[11]} {f[12]} | **{f[13]}** | {f[14]} | {f[15]} |")
doc = between(doc, "for `c_1` and `c_N`; they ignore the truncation of the ansatz.", "## 8. What remains", "\n".join(rows))
open(D + 'SPINNING_LOWK.md', 'w').write(doc)
print(re.findall(r'@@[A-Z0-9]+@@', doc))
